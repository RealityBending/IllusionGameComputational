"""Generate Muller-Lyer illustrations for the axes of the parameter heatmaps.

The heatmaps in `analysis/1_modelcomparison.qmd` are drawn over the raw units of
the Muller-Lyer stimulus space: illusion strength (x) from -49 to 49 degrees of
fin tilt, and task difficulty (y) from 0.04 to 0.46 of relative line length.
Those are the observed bounds of the data, and the ones `restore_units()` in the
.qmd files inverts back to (see `analysis/server/fit_model.R`).

This script renders the two ends and the middle of each of those axes, so the
figure can show what the axes actually mean:

- `mullerlyer_difference{1,2,3}.png`, at zero illusion strength: the objective
  length difference goes from barely visible (hardest) to obvious (easiest);
- `mullerlyer_strength{1,2,3}.png`, at mid difficulty: the fins go from
  facilitating (negative, they exaggerate the real difference) through absent to
  conflicting (positive, they fight it). The sign convention is pyllusion's and
  is the one the task data uses -- `fit_model.R` calls `Illusion_Strength >= 0`
  "Conflicting".

Both series are numbered in the direction of their axis, low value to high, so
the .qmd can load them in a loop (and reverse the difference one, the y axis
running bottom to top). The middle image of each series is the same stimulus
(mid difficulty, no illusion); it is written twice so each series can be read on
its own.

Usage
-----
    python analysis/illustrations/make_illustrations.py

pyllusion is not on PyPI in this environment, so a local checkout is used. Set
PYLLUSION_PATH if yours is elsewhere.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

try:
    import pyllusion
except ModuleNotFoundError:
    sys.path.insert(
        0,
        os.environ.get("PYLLUSION_PATH", r"C:\Users\domma\Dropbox\Software\Pyllusion"),
    )
    import pyllusion

from PIL import Image, ImageChops

# Observed bounds of the Muller-Lyer trials in data/illusion_part*.csv. Keep in
# step with restore_units() in the .qmd files if the stimulus set ever changes.
DIFFERENCE_MIN, DIFFERENCE_MAX = 0.04, 0.46
STRENGTH_MIN, STRENGTH_MAX = -49.0, 49.0
DIFFERENCE_MID = (DIFFERENCE_MIN + DIFFERENCE_MAX) / 2

DIFFERENCES = [DIFFERENCE_MIN, DIFFERENCE_MID, DIFFERENCE_MAX]
STRENGTHS = [STRENGTH_MIN, 0.0, STRENGTH_MAX]

WIDTH, HEIGHT = 800, 600
OUTLINE = 12  # line width in px, at the size above
PAD = 16  # px of white kept around the cropped stimuli


def render(illusion_strength: float, difference: float) -> Image.Image:
    illusion = pyllusion.MullerLyer(
        illusion_strength=illusion_strength, difference=difference
    )
    return illusion.to_image(width=WIDTH, height=HEIGHT, outline=OUTLINE)


def union_bbox(images: list[Image.Image]) -> tuple[int, int, int, int]:
    """Smallest box containing the ink of every image.

    One box for all six, never one per image: cropping each to its own ink would
    rescale them against each other, and the length difference between the lines
    is the very thing these images are meant to show. It also leaves every file
    the same size, so the two strips can be drawn at one scale in the figure.
    """
    boxes = [ImageChops.invert(im.convert("L")).getbbox() for im in images]
    x0 = max(min(b[0] for b in boxes) - PAD, 0)
    y0 = max(min(b[1] for b in boxes) - PAD, 0)
    x1 = min(max(b[2] for b in boxes) + PAD, WIDTH)
    y1 = min(max(b[3] for b in boxes) + PAD, HEIGHT)
    return x0, y0, x1, y1


def main() -> None:
    series = {
        "difference": [
            (value, render(illusion_strength=0, difference=value))
            for value in DIFFERENCES
        ],
        "strength": [
            (value, render(illusion_strength=value, difference=DIFFERENCE_MID))
            for value in STRENGTHS
        ],
    }

    box = union_bbox([image for levels in series.values() for _, image in levels])

    for name, levels in series.items():
        for i, (value, image) in enumerate(levels, start=1):
            path = HERE / f"mullerlyer_{name}{i}.png"
            cropped = image.crop(box)
            cropped.save(path)
            print(f"{path.name}  {value:>6.2f}  {cropped.size}")


if __name__ == "__main__":
    main()
