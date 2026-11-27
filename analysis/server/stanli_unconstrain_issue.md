# `unconstrain()` rejects a length-1 R value for a `vector[1]` parameter (regression since 0.25.0)

Draft of an issue for [seantalts/stanli](https://github.com/seantalts/stanli/issues).
Everything below the line is the issue. Reproduced on Artemis on 2026-11-27
(job 11850756; script `stanli_vector1_reprex.R` in
`IGComputational/stanli_vector1/` on scratch). See `AGENT.md` §6.19.

---

**Title:** `unconstrain()` rejects a length-1 value for a `vector[1]` parameter ("dims declared=(1); dims found=()"), which 0.25.0 accepted

### What happens

In stanli 0.31.0, `unconstrain()` stops when a parameter declared with one
dimension of size 1 (`vector[1]`, `array[1] real`, ...) is given as a plain R
number:

```
mismatch in number dimensions declared and found in context; processing stage=parameter initialization; variable name=tau; dims declared=(1); dims found=()
```

The same value as `array(x, dim = 1)` works. stanli 0.25.0 accepted both.

R has no length-1 vector distinct from a scalar, so a plain number is how
these values come out of most R code. In particular, brms declares the SDs of
each group-level term as `vector<lower=0>[M_k] sd_k`, so any formula with a
lone varying intercept (`(1 | g)`) has a `vector[1]`. We hit this with brms
2.21.0 programs whose inits came from a generic init function. All four jobs
of a full-data fit stopped at initialisation, before sampling. With nutpie
(`sample_nutpie()`, `init` as a function) the same inits went through.

### Reprex

```r
library(stanli)
m <- stanli_model(code = "
parameters {
  real mu;
  vector<lower=0>[1] tau;
}
model {
  mu ~ normal(0, 1);
  tau ~ normal(0, 1);
}")
unconstrain(m, list(mu = 0.1, tau = 0.5))                          # error in 0.31.0
unconstrain(m, list(mu = 0.1, tau = 0.5), partial = TRUE)          # error in 0.31.0
unconstrain(m, list(mu = 0.1, tau = array(0.5, dim = 1)))          # works
```

Output, on the same node, with each version's release runtime:

```
stanli 0.25.0
tau = 0.5                           OK: 0.1 -0.6931
tau = array(0.5, dim = 1)           OK: 0.1 -0.6931
tau = 0.5, partial = TRUE           OK: 0.1 -0.6931
tau = array(0.5, 1), partial = TRUE OK: 0.1 -0.6931

stanli 0.31.0
tau = 0.5                           ERROR: mismatch in number dimensions declared and found in context; processing stage=parameter initialization; variable name=tau; dims declared=(1); dims found=()
tau = array(0.5, dim = 1)           OK: 0.1 -0.6931
tau = 0.5, partial = TRUE           ERROR: mismatch in number dimensions declared and found in context; processing stage=parameter initialization; variable name=tau; dims declared=(1); dims found=()
tau = array(0.5, 1), partial = TRUE OK: 0.1 -0.6931
```

### Expected

That a length-1 value is accepted for a parameter declared with one
dimension of size 1, as in 0.25.0. cmdstanr does this deliberately for inits.
`process_init.list()` (0.9.0) checks each parameter's declared dimensions and
wraps such a value before writing the JSON:

```r
if (model_variables$parameters[[par_name]]$dimensions == 1 && length(init[[i]][[par_name]]) == 1) {
  init[[i]][[par_name]] <- array(init[[i]][[par_name]], dim = 1)
}
```

`unconstrain()` already knows the declared dimensions (the error prints them),
so it could do the same. A length-1 value given for a declared scalar should
of course stay a scalar.

We have not checked 0.26.0 to 0.30.x, nor whether the other entry points
that take R values for parameters behave the same (`sample_model(init = )`
with a list, `brm(backend = "stanli", init = )`).

### Environment

- stanli 0.31.0 with the `stanli-runtime-linux-x86_64` release runtime
  (`STANLI_RUNTIME`), against 0.25.0 with its own release runtime
- R 4.3.2 (EasyBuild, GCCcore 12.3.0), Linux x86_64 (a SLURM cluster)
- stanli built from the release tarball with `-O2 -march=x86-64-v3`

### Workaround

Box the value before the call: `array(x, dim = 1)` for every parameter whose
declaration has one dimension and whose value has length 1.
