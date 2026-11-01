"""Draw the Illusion Game logo.

Two open chevrons and an eye, each part carrying a meaning:

- the red chevron points up: bottom-up sensory evidence;
- the blue chevron points down: top-down priors;
- the purple eye between them is the conscious, illusory percept, which in a
  predictive coding account is the integration of the two. Its purple is the
  midpoint of the red and the blue in OKLCh (a perceptual colour space), so it
  is their mix rather than a third colour picked by eye. The two points where
  the chevrons cross are drawn in it too: that is where they meet. The iris is
  a plain disc of a darker shade of the same purple.

The chevrons are the two triangles of a hexagram with their bases taken away.
The echo of Solomon's seal -- the union of opposites -- is wanted, but the star
is left open so that it does not read as the sigil. What remains of each base
is its middle third, the stretch that lies inside the other chevron, drawn in
black and tucked under that chevron's legs: two equal bars across converging
lines, which is the Ponzo illusion, once in each chevron. Each bar is level with
the cut ends of the chevron whose base it was.

The eye is a vesica piscis, the lens where two circles overlap and the other
classical figure for the meeting of two worlds. It is also the one almond that
fits the star: its tips have the same 120 degree angle as the crossings, so
with EYE_SCALE = 1 the eye's tips sit on the crossings and the upper lid leaves
each tip tangent to the red leg, the lower lid tangent to the blue one. Below 1
the eye nests in those corners with its lids parallel to the legs.

The GIF tells the same story in time. It opens on the priors alone: the blue
chevron and a blue eye, with the top bar longer than the bottom one -- an
eyebrow over the eye, and the two bars seen as unequal. The red chevron then
rises in from below, the eye turns purple, and the top bar shrinks until the two
are equal, the percept the evidence corrects it to. After a hold on the logo,
the red fades out still drifting up, and takes the correction with it, back to
where the loop began. The evidence only ever moves up, so the seam does not
show.

Usage
-----
    python img/make_logo.py

Writes logo.svg, logo.png and logo.gif (black bars, for light backgrounds) and
logo_dark.svg, logo_dark.png and logo_dark.gif (light bars, for dark
backgrounds) next to this script. The SVGs and PNGs are on a transparent
background. The GIFs are on an opaque one, since GIF transparency is all or
nothing and would leave the anti-aliased edges ragged.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Circle, Polygon
from PIL import Image

HERE = Path(__file__).resolve().parent

RED = "#E53935"
BLUE = "#1E88E5"
# Colour of the bars, and the background of the GIF
THEMES = {"logo": ("#1A1A1A", "#FFFFFF"), "logo_dark": ("#F2F2F2", "#121212")}

# Every length is in units of the hexagram's circumradius, R = 1.
STROKE = 0.075  # line width
EYE_SCALE = 1.0  # 1 puts the tips of the eye on the crossings
IRIS = 0.8  # iris radius, as a share of the eye's half-height
IRIS_SHADE = -0.2  # OKLab lightness of the iris relative to the purple
MARGIN = 0.06  # empty border around the mark

PNG_SIZE = 1024  # px; the canvas is square
FIG_INCHES = 4  # only sets the scale of the SVG

GIF_SIZE = 512  # px
# GIF frame delays are whole hundredths of a second, and browsers slow any delay
# under 2 of them down to 10, so 50 fps is as smooth as a GIF plays.
FPS = 50
# Seconds spent on the priors alone, the evidence rising in, the logo, and the
# evidence rising out
PHASES = (1.0, 1.4, 1.5, 1.2)
# How far below its place the red enters, and above it leaves. The way out is
# shorter: there is only MARGIN above the apex, so the edge of the frame starts
# cutting into it at once, and the red has to be faint before that shows.
RISE = (0.6, 0.25)
# Top bar under the priors alone, as a share of its length in the logo
TOP_BAR_PRIOR = 1.5

SQRT3 = np.sqrt(3)
# Each chevron runs end - apex - end. Red is the upward triangle without its
# base, blue the downward one.
RED_CHEVRON = np.array([[-SQRT3 / 2, -0.5], [0, 1], [SQRT3 / 2, -0.5]])
BLUE_CHEVRON = np.array([[-SQRT3 / 2, 0.5], [0, -1], [SQRT3 / 2, 0.5]])


# Colour ----------------------------------------------------------------------

# sRGB <-> OKLab, Ottosson (2020), https://bottosson.github.io/posts/oklab/
_LMS = np.array(
    [
        [0.4122214708, 0.5363325363, 0.0514459929],
        [0.2119034982, 0.6806995451, 0.1073969566],
        [0.0883024619, 0.2817188376, 0.6299787005],
    ]
)
_LAB = np.array(
    [
        [0.2104542553, 0.7936177850, -0.0040720468],
        [1.9779984951, -2.4285922050, 0.4505937099],
        [0.0259040371, 0.7827717662, -0.8086757660],
    ]
)


def to_oklch(colour: str) -> np.ndarray:
    rgb = np.array(to_rgb(colour))
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    lightness, a, b = _LAB @ np.cbrt(_LMS @ linear)
    return np.array([lightness, np.hypot(a, b), np.arctan2(b, a)])


def from_oklch(lch: np.ndarray) -> str:
    lightness, chroma, hue = lch
    lab = np.array([lightness, chroma * np.cos(hue), chroma * np.sin(hue)])
    linear = np.linalg.solve(_LMS, (np.linalg.solve(_LAB, lab)) ** 3)
    linear = np.clip(linear, 0, 1)
    rgb = np.where(
        linear <= 0.0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - 0.055
    )
    return to_hex(rgb)


def mix(colour1: str, colour2: str, weight: float = 0.5) -> str:
    """Interpolate two colours in OKLCh, round the shorter arc of hues.

    Straight through OKLab, red and blue meet in a greyish mauve; going round
    the hue circle keeps the chroma, and the shorter way from red to blue
    passes through purple.
    """
    lch1, lch2 = to_oklch(colour1), to_oklch(colour2)
    dhue = (lch2[2] - lch1[2] + np.pi) % (2 * np.pi) - np.pi
    lch = lch1 + weight * (lch2 - lch1)
    lch[2] = lch1[2] + weight * dhue
    return from_oklch(lch)


PURPLE = mix(RED, BLUE)


def darker(colour: str) -> str:
    """The shade of `colour` that fills the iris."""
    return from_oklch(to_oklch(colour) + [IRIS_SHADE, 0, 0])


# Geometry --------------------------------------------------------------------


def intersect(p, d, q, e) -> np.ndarray:
    """Where the line through `p` along `d` meets the line through `q` along `e`."""
    return p + np.linalg.solve(np.column_stack([d, -e]), q - p)[0] * d


def normal(d) -> np.ndarray:
    return np.array([-d[1], d[0]]) / np.linalg.norm(d)


def chevron(points: np.ndarray, width: float) -> np.ndarray:
    """Outline of a chevron stroked `width` wide, its open ends cut level.

    `points` runs end - apex - end. The apex is mitred; each end is cut along
    the horizontal through it rather than square to its leg.
    """
    end1, apex, end2 = points
    leg1, leg2 = apex - end1, end2 - apex
    level = np.array([1.0, 0.0])
    sides = []
    for offset in (width / 2, -width / 2):
        a, b = offset * normal(leg1), offset * normal(leg2)
        sides.append(
            [
                intersect(end1 + a, leg1, end1, level),
                intersect(end1 + a, leg1, apex + b, leg2),
                intersect(apex + b, leg2, end2, level),
            ]
        )
    return np.array(sides[0] + sides[1][::-1])


def stroke_overlap(p, d, q, e, width: float) -> np.ndarray:
    """Corners of the parallelogram where two straight strokes overlap.

    The strokes run through `p` along `d` and through `q` along `e`.
    """
    n, m = normal(d), normal(e)
    return np.array(
        [
            intersect(p + s * width / 2 * n, d, q + t * width / 2 * m, e)
            for s, t in [(1, 1), (1, -1), (-1, -1), (-1, 1)]
        ]
    )


def vesica(half_width: float, n: int = 200) -> np.ndarray:
    """Outline of a vesica piscis with its tips at (+-half_width, 0).

    Each lid is an arc of a circle of radius 2w/sqrt(3) centred at
    (0, -+w/sqrt(3)): the circles pass through each other's centres, which is
    what makes the tips 120 degrees.
    """
    radius = 2 * half_width / SQRT3
    upper = np.linspace(np.pi / 6, 5 * np.pi / 6, n)
    # Without the tips, which the upper lid already has: a repeated vertex
    # leaves the mitre at the tip without a direction.
    lower = upper[1:-1] + np.pi
    angles = np.concatenate([upper, lower])
    centres = np.where(angles < np.pi, -radius / 2, radius / 2)
    return np.column_stack([radius * np.cos(angles), radius * np.sin(angles) + centres])


# Drawing ---------------------------------------------------------------------


def canvas(size_px: int, background: str | None = None):
    """A square figure filled by its axes, and the stroke width in points."""
    # The apexes are mitred: at 60 degrees the point reaches one stroke width
    # past the apex, which sets the height of the mark.
    half_span = 1 + STROKE + MARGIN
    fig = plt.figure(
        figsize=(FIG_INCHES, FIG_INCHES),
        dpi=size_px / FIG_INCHES,
        facecolor=background or "none",
    )
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-half_span, half_span)
    ax.set_ylim(-half_span, half_span)
    ax.set_aspect("equal")
    ax.axis("off")
    lw = STROKE * FIG_INCHES * 72 / (2 * half_span)  # data units -> points
    return fig, ax, lw


def eye_colour(evidence: float) -> str:
    """Blue under the priors alone, purple once the evidence is all in."""
    return PURPLE if evidence == 1 else mix(BLUE, PURPLE, evidence)


def draw_mark(ax, lw: float, bar_colour: str, evidence=1.0, rise=0.0) -> None:
    """Draw the mark, with as much of the evidence in it as `evidence` says.

    At 0 there are only the priors: no red, a blue eye, and a top bar
    TOP_BAR_PRIOR times as long as in the logo. Up to 1, the logo, the red fades
    in while the eye turns purple and the top bar comes to its length. `rise`
    moves the red chevron up (or down, negative) from its place. The red's
    opacity is evidence squared, so it stays faint while it is far from there.
    """
    colour = eye_colour(evidence)
    opacity = evidence**2

    # Layers, bottom to top: the eye and bars, the chevrons, the crossings.
    # Under the chevrons, lids that are tangent to the legs near the tips
    # peel out from behind them instead of riding over them.
    eye = EYE_SCALE / SQRT3  # the crossings are at (+-1/sqrt(3), 0)
    ax.add_patch(
        Polygon(
            vesica(eye),
            closed=True,
            fill=False,
            edgecolor=colour,
            lw=lw,
            joinstyle="miter",
            zorder=1,
        )
    )
    iris = IRIS * eye / SQRT3
    ax.add_patch(Circle((0, 0), iris, color=darker(colour), lw=0, zorder=1))

    # Each bar runs to the centrelines of the legs it crosses, which are
    # 1/sqrt(3) apart at y = +-1/2 (the middle third of the missing base).
    # Cut square there (matplotlib's default cap projects), its ends lie wholly
    # under the legs; a round cap would poke out of their far side.
    bar = 1 / (2 * SQRT3)
    top = bar * (TOP_BAR_PRIOR + (1 - TOP_BAR_PRIOR) * evidence)
    for half, y in [(top, 0.5), (bar, -0.5)]:
        ax.plot(
            [-half, half],
            [y, y],
            color=bar_colour,
            lw=lw,
            solid_capstyle="butt",
            zorder=1,
        )

    ax.add_patch(Polygon(chevron(BLUE_CHEVRON, STROKE), color=BLUE, lw=0, zorder=2))
    if evidence == 0:
        return
    red = RED_CHEVRON + [0, rise]
    ax.add_patch(
        Polygon(chevron(red, STROKE), color=RED, alpha=opacity, lw=0, zorder=2)
    )
    for side in (-1, 1):
        red_leg = red[1] - red[1 + side]
        blue_leg = BLUE_CHEVRON[1] - BLUE_CHEVRON[1 + side]
        overlap = stroke_overlap(red[1], red_leg, BLUE_CHEVRON[1], blue_leg, STROKE)
        ax.add_patch(Polygon(overlap, color=PURPLE, alpha=opacity, lw=0, zorder=3))


def draw(path: Path, bar_colour: str, size_px: int = PNG_SIZE) -> None:
    fig, ax, lw = canvas(size_px)
    draw_mark(ax, lw, bar_colour)
    fig.savefig(path, transparent=True)
    plt.close(fig)


# Animation -------------------------------------------------------------------


def ease(t: float) -> float:
    return (1 - np.cos(np.pi * t)) / 2


def timeline() -> list[tuple[float, float]]:
    """(evidence, rise) of every frame of one loop; see draw_mark()."""
    priors, rise_in, logo, rise_out = (round(seconds * FPS) for seconds in PHASES)
    below, above = RISE
    frames = [(0.0, -below)] * priors
    for i in range(1, rise_in + 1):
        t = ease(i / rise_in)
        frames.append((t, -below * (1 - t)))
    frames += [(1.0, 0.0)] * logo
    for i in range(1, rise_out + 1):
        t = ease(i / rise_out)
        frames.append((1 - t, above * t))
    return frames


def to_palette(frame: np.ndarray, palette: np.ndarray) -> Image.Image:
    """`frame` as a paletted image, each colour mapped to its nearest entry.

    Pillow's own mapping goes through a coarse colour cache, which let a near
    white from the median cut stand in for a white background.
    """
    packed = frame.reshape(-1, 3).astype(np.int64) @ [1 << 16, 1 << 8, 1]
    colours, pixels = np.unique(packed, return_inverse=True)
    rgb = np.column_stack([colours >> 16, colours >> 8 & 255, colours & 255])
    distance = (rgb**2).sum(1)[:, None] - 2 * rgb @ palette.T + (palette**2).sum(1)
    indices = distance.argmin(1)[pixels].reshape(frame.shape[:2])
    height, width = indices.shape
    image = Image.frombytes("P", (width, height), indices.astype(np.uint8).tobytes())
    image.putpalette(palette.ravel().tolist())
    return image


def animate(path: Path, bar_colour: str, background: str) -> None:
    fig, ax, lw = canvas(GIF_SIZE, background)
    images = []
    for evidence, rise in timeline():
        for artist in [*ax.patches, *ax.lines]:
            artist.remove()
        draw_mark(ax, lw, bar_colour, evidence, rise)
        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[..., :3]

        # Each frame gets its own palette. Shared by the whole loop, a palette
        # had too few shades for the eye's change of colour and the red's fade,
        # which then moved in uneven steps of one to four frames. The flat
        # colours of the frame go in as they are, so that they hold still from
        # one frame to the next; median cut, which averages the colours in each
        # of its boxes (white came out #FAFAFA), fills the rest.
        colour = eye_colour(evidence)
        exact = [background, bar_colour, RED, BLUE, PURPLE, colour, darker(colour)]
        exact = np.round(np.array([to_rgb(c) for c in exact]) * 255).astype(np.int64)
        n_fill = 256 - len(exact)
        fill = Image.fromarray(frame).quantize(n_fill, method=Image.Quantize.MEDIANCUT)
        fill = np.array(fill.getpalette()[: 3 * n_fill], dtype=np.int64).reshape(-1, 3)
        images.append(to_palette(frame, np.concatenate([exact, fill])))
    plt.close(fig)

    # Pillow merges runs of identical frames into one, so the holds cost nothing
    images[0].save(
        path, save_all=True, append_images=images[1:], duration=1000 // FPS, loop=0
    )


def main() -> None:
    print(f"red {RED}  blue {BLUE}  purple {PURPLE}")
    for name, (bar_colour, background) in THEMES.items():
        for ext in ("svg", "png"):
            path = HERE / f"{name}.{ext}"
            draw(path, bar_colour)
            print(path.name)
        path = HERE / f"{name}.gif"
        animate(path, bar_colour, background)
        print(path.name)


if __name__ == "__main__":
    main()
