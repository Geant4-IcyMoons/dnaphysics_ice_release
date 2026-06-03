#!/usr/bin/env python3
"""
Plot Europa electron energy maps derived from e_bombardment files.

Leading hemisphere uses the per-coordinate minimum energy threshold.
Trailing hemisphere uses the per-coordinate maximum energy (max over energy axis).
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from constants import (
    FONT_COURIER,
    FONTSIZE_24,
    OUTPUT_DIR,
    RC_BASE_ELASTIC,
    TOP_ROOT,
    rcparams_with_fontsize,
)


plt.rcParams["font.family"] = FONT_COURIER
plt.rcParams["mathtext.rm"] = FONT_COURIER
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams.update(rcparams_with_fontsize(RC_BASE_ELASTIC, FONTSIZE_24))


def _load_pickle_array(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing data file: {path}")
    with open(path, "rb") as handle:
        data = pickle.load(handle)
    return np.asarray(data, dtype=float)


def _leading_min_map(path: Path) -> np.ndarray:
    data = _load_pickle_array(path)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D leading map, got shape {data.shape}")
    return data


def _trailing_max_map(path: Path) -> np.ndarray:
    data = _load_pickle_array(path)
    if data.ndim != 3:
        raise ValueError(f"Expected 3D trailing map, got shape {data.shape}")
    return np.max(data, axis=-1)


def _min_positive(arr: np.ndarray) -> float:
    positive = arr[arr > 0]
    if positive.size == 0:
        return 0.0
    return float(np.min(positive))


def _degree_formatter(value: float, _pos: int) -> str:
    if np.isclose(value, 0.0):
        value = 0.0
    if float(value).is_integer():
        label = f"{int(value)}"
    else:
        label = f"{value:g}"
    return rf"${label}\!^\circ$"


def _add_adaptive_grid(
    ax: plt.Axes,
    data: np.ndarray,
    im: plt.AxesImage,
    extent: tuple[float, float, float, float],
    linewidth: float = 0.35,
    alpha: float = 0.55,
) -> None:
    if data.ndim != 2 or data.size == 0:
        return

    n_rows, n_cols = data.shape
    rgba = im.to_rgba(data)
    luminance = (
        0.2126 * rgba[..., 0]
        + 0.7152 * rgba[..., 1]
        + 0.0722 * rgba[..., 2]
    )

    x_edges = np.linspace(extent[0], extent[1], n_cols + 1)
    y_edges = np.linspace(extent[3], extent[2], n_rows + 1)

    segments: list[list[tuple[float, float]]] = []
    colors: list[tuple[float, float, float, float]] = []

    for i in range(n_rows):
        y_top = float(y_edges[i])
        y_bottom = float(y_edges[i + 1])
        for j in range(n_cols):
            x_left = float(x_edges[j])
            x_right = float(x_edges[j + 1])
            is_bright = bool(np.isfinite(luminance[i, j]) and luminance[i, j] > 0.5)
            if is_bright:
                color = (0.0, 0.0, 0.0, alpha)
            else:
                color = (1.0, 1.0, 1.0, alpha)
            segments.extend(
                [
                    [(x_left, y_top), (x_right, y_top)],
                    [(x_left, y_bottom), (x_right, y_bottom)],
                    [(x_left, y_top), (x_left, y_bottom)],
                    [(x_right, y_top), (x_right, y_bottom)],
                ]
            )
            colors.extend([color, color, color, color])

    # grid = LineCollection(segments, colors=colors, linewidths=linewidth, zorder=3)
    # ax.add_collection(grid)


def _contrast_color_from_axes_fraction(
    luminance: np.ndarray,
    x_axes: float,
    y_axes: float,
) -> str:
    n_rows, n_cols = luminance.shape
    col = int(np.clip(np.floor(x_axes * n_cols), 0, n_cols - 1))
    row = int(np.clip(np.floor((1.0 - y_axes) * n_rows), 0, n_rows - 1))
    lum = luminance[row, col]
    if np.isfinite(lum) and lum > 0.5:
        return "black"
    return "white"


def _apply_adaptive_tick_mark_colors(
    ax: plt.Axes,
    data: np.ndarray,
    im: plt.AxesImage,
) -> None:
    if data.ndim != 2 or data.size == 0:
        return

    rgba = im.to_rgba(data)
    luminance = (
        0.2126 * rgba[..., 0]
        + 0.7152 * rgba[..., 1]
        + 0.0722 * rgba[..., 2]
    )

    for tick_value, tick in zip(ax.get_xticks().tolist(), ax.xaxis.get_major_ticks()):
        x_axes, _ = ax.transAxes.inverted().transform(
            ax.transData.transform((tick_value, np.mean(ax.get_ylim())))
        )
        color = _contrast_color_from_axes_fraction(luminance, float(x_axes), 0.02)
        tick.tick1line.set_color(color)
        tick.tick2line.set_color(color)
        tick.tick1line.set_markeredgecolor(color)
        tick.tick1line.set_markerfacecolor(color)
        tick.tick2line.set_markeredgecolor(color)
        tick.tick2line.set_markerfacecolor(color)

    for tick_value, tick in zip(ax.get_yticks().tolist(), ax.yaxis.get_major_ticks()):
        _, y_axes = ax.transAxes.inverted().transform(
            ax.transData.transform((np.mean(ax.get_xlim()), tick_value))
        )
        color = _contrast_color_from_axes_fraction(luminance, 0.02, float(y_axes))
        tick.tick1line.set_color(color)
        tick.tick2line.set_color(color)
        tick.tick1line.set_markeredgecolor(color)
        tick.tick1line.set_markerfacecolor(color)
        tick.tick2line.set_markeredgecolor(color)
        tick.tick2line.set_markerfacecolor(color)


def _plot_map(
    ax: plt.Axes,
    data: np.ndarray,
    title: str | None,
    cmap_name: str,
    vmin: float | None,
    vmax: float | None,
    under: str | None,
    over: str | None,
    extent: tuple[float, float, float, float],
    x_label: str,
    y_label: str | None,
    xticks: list[float] | None = None,
    norm: BoundaryNorm | None = None,
) -> plt.Axes:
    cmap = plt.get_cmap(cmap_name, 10).copy()
    if under is not None:
        cmap.set_under(under)
    if over is not None:
        cmap.set_over(over)
    imshow_kwargs = {
        "origin": "upper",
        "cmap": cmap,
        "extent": extent,
        "aspect": "equal",
    }
    if norm is not None:
        imshow_kwargs["norm"] = norm
    else:
        imshow_kwargs["vmin"] = vmin
        imshow_kwargs["vmax"] = vmax
    im = ax.imshow(data, **imshow_kwargs)
    _add_adaptive_grid(ax, data, im, extent)
    if title:
        ax.set_title(title)
    ax.set_xlabel(x_label)
    if y_label:
        ax.set_ylabel(y_label)
    if xticks is None:
        xticks = [0, 30, 60, 90, 120, 150, 180]
    ax.set_xticks(xticks)
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90])
    ax.xaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.yaxis.set_major_formatter(FuncFormatter(_degree_formatter))
    ax.tick_params(
        axis="both",
        which="major",
        length=6,
        width=1.2,
        direction="in",
        bottom=True,
        top=False,
        left=True,
        right=False,
        colors="black",
        labelcolor="black",
    )
    return im


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot electron energy maps for Europa."
    )
    parser.add_argument(
        "--leading-path",
        type=Path,
        default=TOP_ROOT / "e_bombardment_leading",
        help="Path to e_bombardment_leading pickle.",
    )
    parser.add_argument(
        "--trailing-path",
        type=Path,
        default=TOP_ROOT / "e_bombardment_trailing",
        help="Path to e_bombardment_trailing pickle.",
    )
    parser.add_argument(
        "--leading-cap",
        type=float,
        default=100.0,
        help="Cap for leading energies; values above use the 'over' color.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUTPUT_DIR / "europa_electron_energy_minmax.png",
        help="Output plot path.",
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="inferno",
        help="Matplotlib colormap.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Save the plot without opening a window.",
    )
    parser.add_argument(
        "--save-npz",
        type=Path,
        default=None,
        help="Optional NPZ output for the derived maps.",
    )
    args = parser.parse_args()

    leading_map = _leading_min_map(args.leading_path)
    trailing_map = _trailing_max_map(args.trailing_path)

    leading_valid = leading_map <= args.leading_cap
    if np.any(leading_valid):
        leading_vmin = float(np.min(leading_map[leading_valid]))
    else:
        leading_vmin = float(np.min(leading_map))
    leading_vmax = float(args.leading_cap)
    leading_bounds = np.linspace(leading_vmin, leading_vmax, 11)
    leading_norm = BoundaryNorm(leading_bounds, ncolors=10, clip=False)

    trailing_vmin = _min_positive(trailing_map)
    trailing_vmax = float(np.max(trailing_map))
    trailing_vmin_plot = trailing_vmin if trailing_vmin > 0 else 0.0
    trailing_bounds = np.linspace(trailing_vmin_plot, trailing_vmax, 11)
    trailing_norm = BoundaryNorm(trailing_bounds, ncolors=10, clip=False)

    print(
        "Leading map: min={:.3f} MeV, max={:.3f} MeV".format(
            leading_vmin, float(np.max(leading_map))
        )
    )
    print(
        "Trailing map: min>0={:.3f} MeV, max={:.3f} MeV".format(
            trailing_vmin, trailing_vmax
        )
    )

    fig, axes = plt.subplots(2, 1, figsize=(10, 16))
    fig.subplots_adjust(left=0.08, right=0.95, bottom=0.07, top=0.97, hspace=0.16)
    leading_extent = (0.0, 180.0, -90.0, 90.0)
    trailing_extent = (180.0, 360.0, -90.0, 90.0)

    im_leading = _plot_map(
        axes[0],
        leading_map,
        None,
        args.cmap,
        leading_vmin,
        leading_vmax,
        under=None,
        over=plt.get_cmap(args.cmap)(1.0),
        extent=leading_extent,
        x_label="",
        y_label=r"Latitude ($^\circ$)",
        xticks=[0, 30, 60, 90, 120, 150, 180],
        norm=leading_norm,
    )
    axes[0].invert_xaxis()
    axes[0].tick_params(labelbottom=False)
    _apply_adaptive_tick_mark_colors(axes[0], leading_map, im_leading)
    axes[0].text(
        0.0,
        1.04,
        "(a)",
        transform=axes[0].transAxes,
        ha="left",
        va="bottom",
        fontsize=float(FONTSIZE_24),
        clip_on=False,
    )
    divider1 = make_axes_locatable(axes[0])
    cax1 = divider1.append_axes("right", size="4%", pad=0.05)
    cbar1 = fig.colorbar(
        im_leading,
        cax=cax1,
        extend="max",
        boundaries=leading_bounds,
        ticks=leading_bounds,
        spacing="proportional",
    )
    cbar1.set_label("Min. energy ($T$; MeV)")

    extend_trailing = "min" if trailing_vmin > 0 else "neither"
    im_trailing = _plot_map(
        axes[1],
        trailing_map,
        None,
        args.cmap,
        trailing_vmin if trailing_vmin > 0 else None,
        trailing_vmax,
        under="black",
        over=None,
        extent=trailing_extent,
        x_label=r"Longitude ($^\circ$ W)",
        y_label=r"Latitude ($^\circ$)",
        xticks=[180, 210, 240, 270, 300, 330, 360],
        norm=trailing_norm,
    )
    axes[1].invert_xaxis()
    _apply_adaptive_tick_mark_colors(axes[1], trailing_map, im_trailing)
    axes[1].text(
        0.0,
        1.04,
        "(b)",
        transform=axes[1].transAxes,
        ha="left",
        va="bottom",
        fontsize=float(FONTSIZE_24),
        clip_on=False,
    )
    divider2 = make_axes_locatable(axes[1])
    cax2 = divider2.append_axes("right", size="4%", pad=0.05)
    cbar2 = fig.colorbar(
        im_trailing,
        cax=cax2,
        extend=extend_trailing,
        boundaries=trailing_bounds,
        ticks=trailing_bounds,
        spacing="proportional",
    )
    cbar2.set_label("Max. energy ($T$; MeV)")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out, bbox_inches="tight")
    print(f"Plot saved to: {args.out}")

    if args.save_npz is not None:
        args.save_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            args.save_npz,
            leading_min_mev=np.asarray(leading_map, dtype=float),
            trailing_max_mev=np.asarray(trailing_map, dtype=float),
        )
        print(f"NPZ saved to: {args.save_npz}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
