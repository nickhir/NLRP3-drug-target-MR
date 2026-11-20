# shared look for all the python figures, same as style.R
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg")
from matplotlib import rcParams, font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullFormatter

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT_DIR = ROOT / "figures_out"
OUT_DIR.mkdir(exist_ok=True)

IVW_C, WM_C = "#C0392B", "#2C6FB5"
GREY = "#8A8A8A"


def apply_style():
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    rcParams["font.family"] = "Open Sans"
    rcParams["font.weight"] = "semibold"
    # Bold and ExtraBold share a weight, the font list order decides which is used
    font_manager.fontManager.ttflist = [
        f for f in font_manager.fontManager.ttflist
        if "OpenSans-ExtraBold" not in f.fname
    ]


def minor_ticks(ax, step, length, width):
    # unlabelled minor ticks every step, call after set_xlim() and set_xticks()
    lo, hi = sorted(ax.get_xlim())
    majors = {round(float(t), 9) for t in ax.get_xticks()}
    ticks = [round(n * step, 9)
             for n in range(math.ceil(lo / step - 1e-9), math.floor(hi / step + 1e-9) + 1)]
    ax.xaxis.set_minor_locator(FixedLocator([t for t in ticks if t not in majors]))
    # a log axis spanning less than a decade would otherwise label its minors
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="minor", length=length, width=width)


def legend_entry(fig, cx, y, colour, text, half, pad, lw, ms, mew, fontsize):
    # one legend entry, a line through an open diamond and then the label
    fig.add_artist(Line2D([cx - half, cx + half], [y, y], color=colour,
                          linewidth=lw, transform=fig.transFigure))
    fig.add_artist(Line2D([cx], [y], marker="D", markersize=ms, color=colour,
                          markerfacecolor="white", markeredgewidth=mew,
                          linestyle="none", transform=fig.transFigure))
    fig.text(cx + half + pad, y, text, fontsize=fontsize, va="center")


def centred_legend(fig, items, y, centre, gap, half, pad, lw, ms, mew, fontsize):
    # a row of legend entries centred on x = centre
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    widths = []
    for _, text in items:
        probe = fig.text(0, 0, text, fontsize=fontsize)
        widths.append(probe.get_window_extent(renderer=renderer).width / fig.bbox.width)
        probe.remove()
    x = centre - (sum(2 * half + pad + w for w in widths) + gap * (len(items) - 1)) / 2
    for (colour, text), w in zip(items, widths):
        legend_entry(fig, x + half, y, colour, text, half, pad, lw, ms, mew, fontsize)
        x += 2 * half + pad + w + gap
