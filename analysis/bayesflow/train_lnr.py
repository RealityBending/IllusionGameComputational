import numpy as np
import pandas as pd
import bayesflow as bf
import multiprocessing
import joblib

# ===============================================================
# 1. Parameters
# ===============================================================
# Load your actual data
df: pd.DataFrame = pd.read_csv(
    "https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/illusion_part1.csv"
)

# Global constants mapping to your data reality
# N_OBS = 120  # Number of trials for a single participant
N_OBS = 1000
MIN_RT = df["RT"].min()


# --- Link Functions ---
def softplus(x):
    """Link function for sigma parameters (lb=0)"""
    # Clipped to prevent np.exp overflow from extreme prior draws
    return np.log1p(np.exp(np.clip(x, -50, 50)))


def expit(x):
    """Inverse-logit link function for tau (bounds 0-1)"""
    return 1 / (1 + np.exp(-np.clip(x, -50, 50)))


# ===============================================================
# 2. Generative model (prior + likelihood/simulator)
# ===============================================================
def prior_and_likelihood(seed, n_obs, min_rt):
    """
    Generates a single trial dataset (parameters + simulated data) on a local CPU core.
    """
    # 1. Initialize a mathematically isolated RNG for a specific core
    rng = np.random.default_rng(seed)

    # ---------------------------------------------------------------
    # 2. Prior: Draw all 10 target parameters
    # ---------------------------------------------------------------
    b_Intercept = rng.normal(0, 10)
    b_Illusion_Difference = rng.normal(0, 10)
    b_nuone_Intercept = rng.normal(0, 10)
    b_nuone_Illusion_Difference = rng.normal(0, 10)
    b_sigmazero_Intercept = rng.normal(0, 10)
    b_sigmazero_Illusion_Difference = rng.normal(0, 10)
    b_sigmaone_Intercept = rng.normal(0, 10)
    b_sigmaone_Illusion_Difference = rng.normal(0, 10)
    b_tau_Intercept = rng.normal(0, 10)
    b_tau_Illusion_Difference = rng.normal(0, 10)

    # ---------------------------------------------------------------
    # 3. Likelihood: Simulate the Lognormal Race
    # ---------------------------------------------------------------
    x = rng.uniform(0, 2, size=n_obs)

    # Linear Equations
    mu_val = b_Intercept + b_Illusion_Difference * x
    nuone_val = b_nuone_Intercept + b_nuone_Illusion_Difference * x

    # Inlined Link Functions (Softplus & Expit)
    sigmazero_val = np.log1p(
        np.exp(
            np.clip(
                b_sigmazero_Intercept + b_sigmazero_Illusion_Difference * x, -50, 50
            )
        )
    )
    sigmaone_val = np.log1p(
        np.exp(
            np.clip(b_sigmaone_Intercept + b_sigmaone_Illusion_Difference * x, -50, 50)
        )
    )
    tau_val = 1 / (
        1 + np.exp(-np.clip(b_tau_Intercept + b_tau_Illusion_Difference * x, -50, 50))
    )

    ndt = tau_val * min_rt

    # SAFEGUARDS: Clip parameters
    safe_mean_0 = np.clip(-mu_val, a_min=-10, a_max=5)
    safe_mean_1 = np.clip(-nuone_val, a_min=-10, a_max=5)
    safe_sig_0 = np.clip(sigmazero_val, a_min=1e-3, a_max=2)
    safe_sig_1 = np.clip(sigmaone_val, a_min=1e-3, a_max=2)

    # Simulated Races
    draws0 = rng.lognormal(mean=safe_mean_0, sigma=safe_sig_0)
    draws1 = rng.lognormal(mean=safe_mean_1, sigma=safe_sig_1)

    dec = (draws1 < draws0).astype(np.float32)
    rt_decision = np.minimum(draws0, draws1)

    rt = np.clip(rt_decision + ndt, a_min=0, a_max=60.0).astype(np.float32)

    # ---------------------------------------------------------------
    # 4. Return Flat Dictionary containing BOTH params and data
    # ---------------------------------------------------------------
    return {
        # Wrap scalars into 1D arrays of shape (1,)
        "b_Intercept": np.array([b_Intercept], dtype=np.float32),
        "b_Illusion_Difference": np.array([b_Illusion_Difference], dtype=np.float32),
        "b_nuone_Intercept": np.array([b_nuone_Intercept], dtype=np.float32),
        "b_nuone_Illusion_Difference": np.array(
            [b_nuone_Illusion_Difference], dtype=np.float32
        ),
        "b_sigmazero_Intercept": np.array([b_sigmazero_Intercept], dtype=np.float32),
        "b_sigmazero_Illusion_Difference": np.array(
            [b_sigmazero_Illusion_Difference], dtype=np.float32
        ),
        "b_sigmaone_Intercept": np.array([b_sigmaone_Intercept], dtype=np.float32),
        "b_sigmaone_Illusion_Difference": np.array(
            [b_sigmaone_Illusion_Difference], dtype=np.float32
        ),
        "b_tau_Intercept": np.array([b_tau_Intercept], dtype=np.float32),
        "b_tau_Illusion_Difference": np.array(
            [b_tau_Illusion_Difference], dtype=np.float32
        ),
        # Reshape 1D data arrays into 2D column vectors of shape (N_OBS, 1)
        "x": x.astype(np.float32).reshape(-1, 1),
        "rt": rt.astype(np.float32).reshape(-1, 1),
        "dec": dec.astype(np.float32).reshape(-1, 1),
    }


def parallel_batched_simulator(batch_shape):
    # CRITICAL FIX: Extract the integer if BayesFlow passes a tuple like (64,)
    batch_size = batch_shape[0] if isinstance(batch_shape, tuple) else batch_shape

    num_cores = 4

    ss = np.random.SeedSequence()
    # Now this is guaranteed to be a pure integer
    child_seeds = ss.spawn(batch_size)

    # (Note: Ensure this matches the name of your worker function!)
    results = joblib.Parallel(n_jobs=num_cores)(
        joblib.delayed(prior_and_likelihood)(seed, N_OBS, MIN_RT)
        for seed in child_seeds
    )

    batched_dict = {
        key: np.array([res[key] for res in results]) for key in results[0].keys()
    }
    return batched_dict


# Tell BayesFlow to use the batched simulator
# (Check your specific BayesFlow version's kwargs: it may be `is_batched=True` or `simulator_is_batched=True`)
simulator = bf.make_simulator(parallel_batched_simulator, is_batched=True)

# ===============================================================
# 3. Adapter + networks + approximator, then train
# ===============================================================
# 1. Generate a single dummy dataset to inspect the keys
dummy_sample = prior_and_likelihood(seed=42, n_obs=N_OBS, min_rt=MIN_RT)

# 2. Dynamically extract all keys EXCEPT the data arrays ("x", "rt", "dec")
target_params = [key for key in dummy_sample.keys() if key not in ["x", "rt", "dec"]]

adapter = (
    bf.Adapter()
    .as_set(["x", "rt", "dec"])
    .convert_dtype("float64", "float32")
    .concatenate(["x", "rt", "dec"], into="summary_variables")
    .concatenate(target_params, into="inference_variables")
)

# DeepSet is ideal for un-ordered trial data
summary_net = bf.networks.DeepSet(summary_dim=64)

# CouplingFlow handles the actual posterior mapping
inference_net = bf.networks.CouplingFlow()

workflow = bf.BasicWorkflow(
    simulator=simulator,
    adapter=adapter,
    inference_network=inference_net,
    summary_network=summary_net,
)

# Train the network
history = workflow.fit_online(
    epochs=10,
    num_batches_per_epoch=600,
    batch_size=64,
)

# Run Diagnostics -----------------------------------------------
# RMSE: Closer to 0 is perfect. A value of 1.0 means the network is just blindly guessing the mean of your prior distribution.
# Posterior Contraction: Measures how much the network "learned" from the data. A value of 1.0 means absolute certainty (the posterior shrank to a single point). A value of 0.0 means the posterior is identical to the prior (the data provided zero new information).
# Calibration Error: Measures whether the model's confidence matches reality. Values < 0.05 are generally considered excellent.
# Log Gamma: A measure of posterior sharpness. Positive values indicate a concentrated posterior; negative values indicate a highly dispersed one.
# metrics = workflow.compute_default_diagnostics(test_data=300)
# print(metrics)
fig = workflow.plot_default_diagnostics(test_data=300)


# ==============================================================
# 4. Fit to Data
# ==============================================================

# Extract the relevant columns, ensure float32, and ADD A BATCH DIMENSION (1, N)
# .reshape(1, -1) turns an array of shape (N,) into (1, N)
data = {
    "x": df["Illusion_Difference"].values.astype(np.float32).reshape(1, -1),
    "rt": df["RT"].values.astype(np.float32).reshape(1, -1),
    "dec": df["Error"].values.astype(np.float32).reshape(1, -1),
}

# Generate posterior samples using the trained workflow
posterior_samples = workflow.sample(conditions=data, num_samples=2000)

# ==============================================================
# 5. Process Posterior Samples and Save to CSV
# ==============================================================
post_means = []
post_sds = []
components = []

# Dynamically extract the parameter names directly from the BayesFlow output!
inferred_params = list(posterior_samples.keys())

for param in inferred_params:
    # Extract the array and squeeze out the batch dimension -> (2000,)
    samples_1d = posterior_samples[param].squeeze()

    # Calculate stats
    post_means.append(np.mean(samples_1d))
    post_sds.append(np.std(samples_1d))

    # Parse the component name based on brms naming conventions
    parts = param.split("_")
    if parts[1] in ["Intercept", "Illusion"]:
        components.append("nuzero")
    else:
        components.append(parts[1])

# Create the summary DataFrame using the dynamically extracted names
results_df = pd.DataFrame(
    {
        "Parameter": inferred_params,
        "Component": components,
        "Mean": post_means,
        "SD": post_sds,
    }
)

# Display the summary in the console
print("\n=== Posterior Summaries ===")
print(results_df.to_string(index=False))

# Save to CSV
output_filename = "posterior_summaries.csv"
results_df.to_csv(output_filename, index=False)
print(f"\nSuccessfully saved posterior summaries to '{output_filename}'.")
