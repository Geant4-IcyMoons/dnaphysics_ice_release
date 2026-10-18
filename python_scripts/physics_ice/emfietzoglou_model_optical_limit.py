r"""
Dielectric model for water ice (amorphous, hexagonal) with q-dispersion.

What this implements
- 5 excitation bands with per-band dispersion parameters (a_j, b_j, c_j) applied to f_j(q)
- 4 ionization shells with global dispersion for E_i(q) and gamma_i(q)
- Optional O K-shell kept optical
- Returns \epsilon_1(E,q), \epsilon_2(E,q), and ELF(E,q) = Im[-1/\\epsilon(E,q)]

Units
- Energies in eV
- q in a0^{-1}
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from constants import (
    FONT_COURIER,
    FONTSIZE_18,
    HFONT_COURIER,
    ICE_DATA_XLSX_PATH,
    RC_BASE_STANDARD,
    rcparams_with_fontsize,
)

font = FONT_COURIER
hfont = HFONT_COURIER
plt.rcParams['font.family'] = font
plt.rcParams['mathtext.rm'] = font
plt.rcParams['mathtext.fontset'] = 'custom'

FONTSIZE = FONTSIZE_18
plt.rcParams.update(rcparams_with_fontsize(RC_BASE_STANDARD, FONTSIZE))


import emfietzoglou_model_finite_q as finite_q
Osc, IceOpticalSet = finite_q.Osc, finite_q.IceOpticalSet
DispersionCoeffs = finite_q.DispersionCoeffs
epsilon_optical = finite_q.epsilon_optical
epsilon1_valence_E0 = finite_q.epsilon1_valence_E0
epsilon2_Kshell_E0 = finite_q.epsilon2_Kshell_E0


def epsilon2_valence_E0(E, s, partitioned=True):
    """Use the production response; do not maintain a separate optical fit."""
    result = finite_q.epsilon2_valence_Eq(np.atleast_1d(E), 0., s,
        DispersionCoeffs(0., 0., 0.), partitioned=partitioned)
    return {k: ([a[0] for a in v] if isinstance(v, list) else v[0])
            for k, v in result.items()}


def elf_E0(E, s, include_kshell=True, partitioned=True):
    real = epsilon1_valence_E0(E, s)["total"]
    imag = epsilon2_valence_E0(E, s, partitioned)["total"]
    valence = imag/(real*real + imag*imag)
    return valence + (epsilon2_Kshell_E0(E, s) if include_kshell else 0.)

def plot_Im_epsilon_channel_resolved(
    E: np.ndarray,
    s: IceOpticalSet,
    C: DispersionCoeffs,  # unused; kept for API compatibility
    q: float = 0.0,       # unused; always optical (q=0)
    include_kshell: bool = False,  # unused here
    ax=None,
    alpha: float = 0.95,
    linewidth: float = 1.8,
    legend: bool = False,
    also_plot_composite: bool = False,
    partitioned: bool = False
):
    """Overlay per-channel ε₂ at q=0 (valence only)."""
    if ax is None:
        ax = plt.gca()
    E = np.asarray(E, float)

    res = epsilon2_valence_E0(E, s, partitioned=partitioned)

    exc_colors = ["#1f77b4", "#2ca02c", "#17becf", "#8c564b", "#9467bd"]
    ion_colors = ["#d62728", "#ff7f0e", "#bcbd22", "#e377c2"]

    for j, y in enumerate(res["excitations"]):
        ax.plot(E, y, color=exc_colors[j % len(exc_colors)], alpha=alpha, lw=linewidth, label=f"Exc {j+1}")
    for k, y in enumerate(res["ionizations"]):
        ax.plot(E, y, color=ion_colors[k % len(ion_colors)], alpha=alpha, lw=linewidth, label=f"Ion {k+1}")

    if also_plot_composite:
        ax.plot(E, res["total"], color="k", lw=2.2, label=r"Total $\epsilon_2$ (valence, q=0)")

    if legend:
        ax.legend(loc='upper right')

    ax.set_xlabel("Energy Transfer ($E$; eV)")
    ax.set_ylabel(r"$\epsilon_2(E, q{=}0)$")

def plot_Re_epsilon_channel_resolved(
    E: np.ndarray,
    s: IceOpticalSet,
    C: DispersionCoeffs,  # unused; kept for API compatibility
    q: float = 0.0,       # unused; always optical (q=0)
    include_kshell: bool = False,  # unused here
    ax=None,
    alpha: float = 0.95,
    linewidth: float = 1.8,
    legend: bool = False,
    also_plot_composite: bool = False,
    include_baseline_one: bool = False,
    partitioned: bool = False  # Added for compatibility, though Re is not partitioned per se
):
    """Overlay per-channel ε₁ at q=0 (valence only)."""
    if ax is None:
        ax = plt.gca()
    E = np.asarray(E, float)

    res = epsilon1_valence_E0(E, s)

    exc_colors = ["#1f77b4", "#2ca02c", "#17becf", "#8c564b", "#9467bd"]
    ion_colors = ["#d62728", "#ff7f0e", "#bcbd22", "#e377c2"]

    for j, y in enumerate(res["excitations"]):
        ax.plot(E, y, color=exc_colors[j % len(exc_colors)], alpha=alpha, lw=linewidth, label=f"Exc {j+1}")
    for k, y in enumerate(res["ionizations"]):
        ax.plot(E, y, color=ion_colors[k % len(ion_colors)], alpha=alpha, lw=linewidth, label=f"Ion {k+1}")

    if also_plot_composite:
        ax.plot(E, res["total"], color="k", lw=2.2, label=r"Total $\epsilon_1$ (valence, q=0)")

    if legend:
        ax.legend(loc='upper right')

    ax.set_xlabel("Energy Transfer ($E$; eV)")
    ax.set_ylabel(r"$\epsilon_1(E, q{=}0)$")

def plot_ELF_channel_resolved(
    E: np.ndarray,
    s: IceOpticalSet,
    C: DispersionCoeffs,  # unused; kept for API compatibility
    q: float = 0.0,       # unused; always optical (q=0)
    include_kshell: bool = True,
    ax=None,
    alpha: float = 0.95,
    linewidth: float = 1.8,
    also_plot_composite: bool = True,
    legend: bool = True,
    partitioned: bool = False
):
    """Channel-resolved ELF at q=0 from channel-resolved ε1, ε2."""
    if ax is None:
        ax = plt.gca()
    E = np.asarray(E, float)

    # Channel-resolved ε at q=0 (must already be implemented as dicts)
    e1 = epsilon1_valence_E0(E, s)     # {"excitations":[...], "ionizations":[...], "total":...}
    e2 = epsilon2_valence_E0(E, s, partitioned=partitioned)     # {"excitations":[...], "ionizations":[...], "total":...}

    denom = e1["total"]**2 + e2["total"]**2
    if np.any(denom <= 0):
        raise ValueError("Non-positive |ε|^2 encountered; check ε construction.")

    # Per-channel ELF contributions: each valence channel divided by |ε^(m)|^2
    exc_elf = [y / denom for y in e2["excitations"]]
    ion_elf = [y / denom for y in e2["ionizations"]]

    # K-shell (optical) added only to ELF, not to ε
    kshell = epsilon2_Kshell_E0(E, s) if include_kshell else None

    # Totals
    elf_valence = e2["total"] / denom
    elf_total = elf_valence + (kshell if kshell is not None else 0.0)

    # Plotting
    exc_colors = ["#1f77b4", "#2ca02c", "#17becf", "#8c564b", "#9467bd"]
    ion_colors = ["#d62728", "#ff7f0e", "#bcbd22", "#e377c2"]

    for j, y in enumerate(exc_elf):
        ax.plot(E, y, color=exc_colors[j % len(exc_colors)], alpha=alpha, lw=linewidth, label=f"Exc {j+1}")
    for k, y in enumerate(ion_elf):
        ax.plot(E, y, color=ion_colors[k % len(ion_colors)], alpha=alpha, lw=linewidth, label=f"Ion {k+1}")

    if kshell is not None:
        ax.plot(E, kshell, color="0.3", ls="--", lw=linewidth, alpha=alpha, label="O K-shell")

    if also_plot_composite:
        ax.plot(E, elf_total, color="k", lw=2.2, label=r"Total ELF (q=0)")

    if legend:
        ax.legend(loc='upper right')

    ax.set_xlabel("Energy Transfer ($E$; eV)")
    ax.set_ylabel(r"$\mathrm{ELF}=\mathrm{Im}\!\left[\frac{1}{\epsilon}\right]$ at $q=0$")

def plot_Kshell_channel_resolved(
    E: np.ndarray,
    s: IceOpticalSet,
    include_kshell: bool = True,
    ax=None,
    alpha: float = 0.95,
    linewidth: float = 1.8,
    legend: bool = False,
    also_plot_composite: bool = False,  # kept for API symmetry; no effect here
):
    """
    Plot only the optical K-shell ε₂ curve with the same interface and style
    as the other plot_*_channel_resolved functions.

    Notes:
      - Uses epsilon2_Kshell_E0(E, s) which must be HARD-gated at s.kshell.Bth.
      - This curve is NOT included in ε plots; it is a reference for ELF.
    """
    if ax is None:
        ax = plt.gca()
    E = np.asarray(E, float)

    if not include_kshell:
        # Keep interface behavior consistent: do nothing but keep axes labeled.
        ax.set_xlabel("Energy Transfer ($E$; eV)")
        ax.set_ylabel(r"$\epsilon_2(E, q{=}0)$")
        if legend:
            ax.legend(loc='upper right')
        return

    y = epsilon2_Kshell_E0(E, s)  # should be gated at Bth internally

    # Same visual language as others
    ax.plot(E, y, color="0.3", ls="--", lw=linewidth, alpha=alpha, label="O K-shell")

    if legend:
        ax.legend(loc='upper right')

    ax.set_xlabel("Energy Transfer ($E$; eV)")
    ax.set_ylabel(r"$\epsilon_2(E, q{=}0)$")

def plot_neff_and_I(E_min=0.1, E_max=1.0e6, npts=50000, savepath=None, partitioned=False):
    r"""
    Plot N_eff(E_max) and I(E_max) at q=0 for amorphous and hexagonal ice,
    reproducing Fig. 3: valence vs total (valence + K-shell).

    Definitions:
      N_eff(E_max) = (2 / (π Ep^2)) * ∫_0^{E_max} E * Im[1/ε](E,0) dE
      ln I(E_max)  = [∫_0^{E_max} E ln(E) Im[1/ε](E,0) dE] / [∫_0^{E_max} E Im[1/ε](E,0) dE]
    """

    def _cumtrapz(x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        out = np.zeros_like(x)
        dx  = np.diff(x)
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dx)
        return out

    # energy grid
    E = np.logspace(np.log10(E_min), np.log10(E_max), int(npts))

    # get parameter sets
    sets = {
        "amorphous": epsilon_optical("amorphous"),
        "hexagonal": epsilon_optical("hexagonal"),
    }

    results = {}
    for name, s in sets.items():
        # Valence ELF from the model (q=0), WITHOUT K shell
        elf_val = elf_E0(E, s, include_kshell=False, partitioned=partitioned)

        # K-shell ε2 added directly to Im[1/ε] per Eq. (7)
        e2k = epsilon2_Kshell_E0(E, s)

        # Total optical ELF
        elf_tot = elf_val + e2k

        # Running integrals
        S_val = _cumtrapz(E, E * elf_val)
        S_tot = _cumtrapz(E, E * elf_tot)

        # N_eff(E_max)
        neff_val = (2.0 / (np.pi * s.Ep**2)) * S_val
        neff_tot = (2.0 / (np.pi * s.Ep**2)) * S_tot

        # I(E_max)
        num_val = _cumtrapz(E, E * np.log(E) * elf_val)
        num_tot = _cumtrapz(E, E * np.log(E) * elf_tot)

        # Safe ratios; undefined at index 0 where S=0
        with np.errstate(divide="ignore", invalid="ignore"):
            logI_val = np.where(S_val > 0.0, num_val / S_val, np.nan)
            logI_tot = np.where(S_tot > 0.0, num_tot / S_tot, np.nan)
        I_val = np.exp(logI_val)
        I_tot = np.exp(logI_tot)

        results[name] = dict(E=E, neff_val=neff_val, neff_tot=neff_tot,
                             I_val=I_val, I_tot=I_tot, s=s)

    # Report maximum Neff values for each material
    for name, res in results.items():
        neff_val_max = float(np.nanmax(res["neff_val"]))
        neff_tot_max = float(np.nanmax(res["neff_tot"]))
        print(f"[Neff] {name}: max valence = {neff_val_max:.4f}, max total = {neff_tot_max:.4f}")

    # Plot with legend bars underneath each panel
    from matplotlib.gridspec import GridSpec
    # Use a taller figure and relatively larger legend row to avoid overlap
    fig = plt.figure(figsize=(10.5, 6.0))
    gs = GridSpec(2, 2, height_ratios=[3.5, 2.0], hspace=0.25, wspace=0.25)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Panel (a): N_eff
    ax1.set_xscale("log")
    ax1.plot(results["amorphous"]["E"], results["amorphous"]["neff_val"], lw=1.9, label="amorphous — valence")
    ax1.plot(results["amorphous"]["E"], results["amorphous"]["neff_tot"], lw=1.9, ls="--", label="amorphous — total")
    ax1.plot(results["hexagonal"]["E"], results["hexagonal"]["neff_val"], lw=1.9, label="hexagonal — valence")
    ax1.plot(results["hexagonal"]["E"], results["hexagonal"]["neff_tot"], lw=1.9, ls="--", label="hexagonal — total")
    # mark the K edge (use amorphous Bth)
    ks = getattr(sets["amorphous"], "kshell", None)
    ax1.set_ylim(0.0, 1.05)
    ax1.set_xlabel(r"$E_{\max}$ (eV)")
    ax1.set_ylabel(r"$N_{\mathrm{eff}}$")
    ax1.set_title(r"$N_{\mathrm{eff}}$ vs $E_{\max}$")

    # Legend bar for panel (a) — stack entries vertically to avoid any overlap
    ax1_leg = fig.add_subplot(gs[1, 0])
    ax1_leg.axis("off")
    handles1, labels1 = ax1.get_legend_handles_labels()
    ax1_leg.legend(handles1, labels1, loc="center left", ncol=1, frameon=False)

    # Panel (b): I(E_max)
    ax2.set_xscale("log")
    ax2.plot(results["amorphous"]["E"], results["amorphous"]["I_val"], lw=1.9, label="amorphous — valence")
    ax2.plot(results["amorphous"]["E"], results["amorphous"]["I_tot"], lw=1.9, ls="--", label="amorphous — total")
    ax2.plot(results["hexagonal"]["E"], results["hexagonal"]["I_val"], lw=1.9, label="hexagonal — valence")
    ax2.plot(results["hexagonal"]["E"], results["hexagonal"]["I_tot"], lw=1.9, ls="--", label="hexagonal — total")
    ax2.set_xlabel(r"$E_{\max}$ (eV)")
    ax2.set_ylabel("I-value (eV)")
    ax2.set_title(r"$I(E_{\max})$")

    # Legend bar for panel (b) — stack entries vertically to avoid any overlap
    ax2_leg = fig.add_subplot(gs[1, 1])
    ax2_leg.axis("off")
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax2_leg.legend(handles2, labels2, loc="center left", ncol=1, frameon=False)

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight")
    plt.show()

    return results

def plot_model_vs_experiment_two_panel(use_partitioning: bool = True, savepath: str | Path | None = "output/Model_vs_Experiment_both.pdf") -> Path:
    """
    Create a two-panel horizontal comparison (amorphous, hexagonal) of model vs
    experimental dielectric properties (Im eps, Re eps, ELF) at q=0, with a
    matching residual panel beneath each main panel. Legend is centered below
    both columns. Formatting matches the single-panel plots; grid is omitted.
    """
    from matplotlib.gridspec import GridSpec

    ice_data_file = ICE_DATA_XLSX_PATH
    E = np.linspace(1.0, 60.0, 20000)
    vib_colors = plt.get_cmap("plasma", 4)(np.arange(3))

    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(
        3,
        2,
        height_ratios=[3.0, 1.2, 0.8],
        width_ratios=[1.0, 1.0],
        hspace=0.22,
        wspace=0.25,
    )
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    residual_axes = [
        fig.add_subplot(gs[1, 0], sharex=axes[0]),
        fig.add_subplot(gs[1, 1], sharex=axes[1]),
    ]
    legend_ax = fig.add_subplot(gs[2, :])
    legend_ax.axis("off")

    handles_all = []
    labels_all = []

    for ax, rax, ice, panel_label in zip(
        axes,
        residual_axes,
        ["amorphous", "hexagonal"],
        ["(a)", "(b)"],
    ):
        s = epsilon_optical(ice)

        sheet_name = "Hexagonal" if ice == "hexagonal" else "Amorphous"
        df_exp = pd.read_excel(ice_data_file, sheet_name=sheet_name)
        mask_e2 = df_exp["eV"].notna() & df_exp["e2"].notna()
        mask_e1 = df_exp["eV.1"].notna() & df_exp["e1"].notna()
        mask_elf = df_exp["eV.2"].notna() & df_exp["ELF"].notna()
        exp_e2_E = df_exp.loc[mask_e2, "eV"].values
        exp_e2 = df_exp.loc[mask_e2, "e2"].values
        exp_e1_E = df_exp.loc[mask_e1, "eV.1"].values
        exp_e1 = df_exp.loc[mask_e1, "e1"].values
        exp_elf_E = df_exp.loc[mask_elf, "eV.2"].values
        exp_elf = df_exp.loc[mask_elf, "ELF"].values

        e1_dict = epsilon1_valence_E0(E, s)
        e2_dict = epsilon2_valence_E0(E, s, partitioned=use_partitioning)
        e1_model = e1_dict["total"]
        e2_model = e2_dict["total"] + epsilon2_Kshell_E0(E, s)
        elf_model = elf_E0(E, s, include_kshell=True, partitioned=use_partitioning)

        # Model predictions
        ln1 = ax.plot(E, e2_model, "k-", linewidth=2, label=r"Model: Im($\epsilon$)", zorder=5)[0]
        ln2 = ax.plot(E, e1_model, "k--", linewidth=2, label=r"Model: Re($\epsilon$)", zorder=5)[0]
        ln3 = ax.plot(E, elf_model, "k-.", linewidth=2, label="Model: ELF", zorder=5)[0]

        # Experimental
        ln4 = ax.plot(exp_e2_E, exp_e2, "d-", color=vib_colors[0], linewidth=1.5, markersize=5,
                      markerfacecolor=vib_colors[0], markeredgecolor=vib_colors[0],
                      label=r"Exp: Im($\epsilon$)", zorder=10)[0]
        ln5 = ax.plot(exp_e1_E, exp_e1, "s-", color=vib_colors[1], linewidth=1.5, markersize=5,
                      markerfacecolor=vib_colors[1], markeredgecolor=vib_colors[1],
                      label=r"Exp: Re($\epsilon$)", zorder=10)[0]
        ln6 = ax.plot(exp_elf_E, exp_elf, "^-", color=vib_colors[2], linewidth=1.5, markersize=5,
                      markerfacecolor=vib_colors[2], markeredgecolor=vib_colors[2],
                      label="Exp: ELF", zorder=10)[0]

        # Residuals: model - experiment at experimental sampling points.
        e2_resid = np.interp(exp_e2_E, E, e2_model) - exp_e2
        e1_resid = np.interp(exp_e1_E, E, e1_model) - exp_e1
        elf_resid = np.interp(exp_elf_E, E, elf_model) - exp_elf

        rax.plot(
            exp_e2_E,
            e2_resid,
            "d-",
            color=vib_colors[0],
            linewidth=1.2,
            markersize=4,
            markerfacecolor=vib_colors[0],
            markeredgecolor=vib_colors[0],
            zorder=5,
        )
        rax.plot(
            exp_e1_E,
            e1_resid,
            "s-",
            color=vib_colors[1],
            linewidth=1.2,
            markersize=4,
            markerfacecolor=vib_colors[1],
            markeredgecolor=vib_colors[1],
            zorder=5,
        )
        rax.plot(
            exp_elf_E,
            elf_resid,
            "^-",
            color=vib_colors[2],
            linewidth=1.2,
            markersize=4,
            markerfacecolor=vib_colors[2],
            markeredgecolor=vib_colors[2],
            zorder=5,
        )
        rax.axhline(0.0, color="black", linewidth=1.0, zorder=1)

        ax.set_xlim(0, 30)
        ax.set_ylim(0, 2.9)
        ax.set_xlabel("")
        if ice == "amorphous":
            ax.set_ylabel("Dielectric properties")
        else:
            ax.set_ylabel("")
        ax.tick_params(
            axis="both",
            which="major",
            direction="in",
            length=6,
            width=1.0,
            bottom=True,
            top=False,
            left=True,
            right=False,
        )
        ax.minorticks_off()
        # No grid as requested
        ax.tick_params(labelbottom=False)
        ax.text(
            0.0,
            1.04,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=float(FONTSIZE),
            clip_on=False,
        )

        rax.set_xlim(0, 30)
        rax.set_ylim(-0.5, 0.5)
        rax.set_yticks([-0.3, 0.0, 0.3])
        rax.set_xlabel("Energy transfer ($E$; eV)")
        if ice == "amorphous":
            rax.set_ylabel("Residuals")
        else:
            rax.set_ylabel("")
        rax.tick_params(
            axis="both",
            which="major",
            direction="in",
            length=6,
            width=1.0,
            bottom=True,
            top=False,
            left=True,
            right=False,
        )
        rax.minorticks_off()

        handles_all.extend([ln1, ln2, ln3, ln4, ln5, ln6])
        labels_all.extend([
            r"Model: Im($\epsilon$)", r"Model: Re($\epsilon$)", "Model: ELF",
            r"Exp: Im($\epsilon$)", r"Exp: Re($\epsilon$)", "Exp: ELF"
        ])

    # Legend centered below panels
    unique = []
    seen = set()
    for h, lbl in zip(handles_all, labels_all):
        if lbl not in seen:
            unique.append((h, lbl))
            seen.add(lbl)
    handle_by_label = {lbl: h for h, lbl in unique}
    labels_u = (
        r"Model: Im($\epsilon$)",
        r"Exp: Im($\epsilon$)",
        r"Model: Re($\epsilon$)",
        r"Exp: Re($\epsilon$)",
        "Model: ELF",
        "Exp: ELF",
    )
    handles_u = tuple(handle_by_label[lbl] for lbl in labels_u)
    legend_ax.legend(
        handles_u,
        labels_u,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.8),
        ncol=3,
        frameon=False,
        columnspacing=1.5,
        handlelength=2.5,
        borderpad=1.5,
        labelspacing=0.8,
    )

    fig.tight_layout()
    out_path = Path(savepath) if savepath is not None else None
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")
        print(f"[saved] {out_path}")
        plt.show()
    return out_path

if __name__ == "__main__":

    # ice = "hexagonal"
    ice = "amorphous"
    s = epsilon_optical(ice)
    a_vec = np.array([3.82, 2.47, 2.47, 3.01, 2.44])
    b_vec = np.array([0.0272, 0.0295, 0.0311, 0.0111, 0.0633])
    c_vec = np.array([0.098, 0.075, 0.074, 0.765, 0.425])
    C = DispersionCoeffs(a_fj=a_vec, b_fj=b_vec, c_fj=c_vec)  # RR2017 defaults for c_disp,d_disp,b1,b2
    E = np.linspace(1.0, 60.0, 20000)
    qvals = np.array([0.0, 0.1, 0.3, 0.6, 0.9, 2.0])  # a0^{-1}

    # Enable Partitioning for corrected plots
    use_partitioning = True

    # Valence ε and ELF(+K) at multiple q
    e1q_val = epsilon1_valence_E0(E, s)
    e2q_val = epsilon2_valence_E0(E, s, partitioned=use_partitioning)

    q = 0.0

    # Load experimental data from ice data.xlsx

    ice_data_file = Path(__file__).parent.parent / 'tabular' / 'ice data.xlsx'
    sheet_name = 'Hexagonal' if ice == "hexagonal" else 'Amorphous'
    df_exp = pd.read_excel(ice_data_file, sheet_name=sheet_name)
    # Extract experimental data with valid values
    mask_e2 = df_exp['eV'].notna() & df_exp['e2'].notna()
    mask_e1 = df_exp['eV.1'].notna() & df_exp['e1'].notna()
    mask_elf = df_exp['eV.2'].notna() & df_exp['ELF'].notna()
    exp_e2_E = df_exp.loc[mask_e2, 'eV'].values
    exp_e2 = df_exp.loc[mask_e2, 'e2'].values
    exp_e1_E = df_exp.loc[mask_e1, 'eV.1'].values
    exp_e1 = df_exp.loc[mask_e1, 'e1'].values
    exp_elf_E = df_exp.loc[mask_elf, 'eV.2'].values
    exp_elf = df_exp.loc[mask_elf, 'ELF'].values

    from matplotlib.gridspec import GridSpec

    # Separate figure: channel-resolved Im(epsilon) at desired q, with legend bar underneath
    fig3 = plt.figure(figsize=(8.0, 7.5))
    gs3 = GridSpec(2, 1, height_ratios=[3.0, 2.0], hspace=0.30)
    ax3 = fig3.add_subplot(gs3[0, 0])
    plot_Im_epsilon_channel_resolved(E, s, C, q=q, include_kshell=True, legend=False,
                                     also_plot_composite=True, ax=ax3, partitioned=use_partitioning)
    # Add experimental data
    ax3.plot(exp_e2_E, exp_e2, 'd-', color='lightgray', linewidth=1.5, markersize=5,
             label='Experimental data', zorder=10)
    ax3.set_xlabel("Energy Transfer ($E$; eV)")
    ax3.set_ylabel(r"$\operatorname{Im}(\epsilon)$")
    ax3.set_title(f"Channel-resolved Im($\\epsilon$) at q = {q}; {ice} ice")
    # Legend bar below main plot (multi-row/column, fully separated from panel)
    ax3_leg = fig3.add_subplot(gs3[1, 0])
    ax3_leg.axis('off')
    handles3, labels3 = ax3.get_legend_handles_labels()
    ax3_leg.legend(handles3, labels3, loc='center', ncol=3, frameon=False)
    plt.savefig(f"output/Channel_resolved_Im_optical_{ice}.pdf", bbox_inches="tight")
    plt.show()

    # Separate figure: channel-resolved Re(epsilon) at desired q, with legend bar underneath
    fig4 = plt.figure(figsize=(8.0, 7.5))
    gs4 = GridSpec(2, 1, height_ratios=[3.0, 2.0], hspace=0.30)
    ax4 = fig4.add_subplot(gs4[0, 0])
    plot_Re_epsilon_channel_resolved(E, s, C, q=q, include_kshell=True, legend=False,
                                     also_plot_composite=True, ax=ax4, partitioned=use_partitioning)
    # Add experimental data
    ax4.plot(exp_e1_E, exp_e1, 'd-', color='lightgray', linewidth=1.5, markersize=5,
             label='Experimental data', zorder=10)
    ax4.set_xlabel("Energy Transfer ($E$; eV)")
    ax4.set_ylabel(r"$\operatorname{Re}(\epsilon)$")
    ax4.set_title(f"Channel-resolved Re($\\epsilon$) at q = {q}; {ice} ice")
    # Legend bar below main plot (multi-row/column, fully separated from panel)
    ax4_leg = fig4.add_subplot(gs4[1, 0])
    ax4_leg.axis('off')
    handles4, labels4 = ax4.get_legend_handles_labels()
    ax4_leg.legend(handles4, labels4, loc='center', ncol=3, frameon=False)
    plt.savefig(f"output/Channel_resolved_Re_optical_{ice}.pdf", bbox_inches="tight")
    plt.show()

    # Separate figure: channel-resolved ELF at desired q (optical limit), with legend bar underneath
    fig5 = plt.figure(figsize=(8.0, 7.5))
    gs5 = GridSpec(2, 1, height_ratios=[3.0, 2.0], hspace=0.30)
    ax5 = fig5.add_subplot(gs5[0, 0])
    plot_ELF_channel_resolved(E, s, C, q=q, include_kshell=True,
                              legend=False, also_plot_composite=True, ax=ax5, partitioned=use_partitioning)
    # Add experimental data
    ax5.plot(exp_elf_E, exp_elf, 'd-', color='lightgray', linewidth=1.5, markersize=5,
             label='Experimental data', zorder=10)
    ax5.set_title(f"Channel-resolved ELF at q = {q}; {ice} ice")
    # Legend bar below main plot (multi-row/column, fully separated from panel)
    ax5_leg = fig5.add_subplot(gs5[1, 0])
    ax5_leg.axis('off')
    handles5, labels5 = ax5.get_legend_handles_labels()
    ax5_leg.legend(handles5, labels5, loc='center', ncol=3, frameon=False)
    plt.savefig(f"output/Channel_resolved_ELF_optical_{ice}.pdf", bbox_inches="tight")
    plt.show()

    # K-shell only, same look-and-feel; window 400–600 eV
    E_k = np.linspace(400.0, 600.0, 4000)
    figK, axK = plt.subplots()
    plot_Kshell_channel_resolved(E_k, s, include_kshell=True,
                                 legend=True, also_plot_composite=False, ax=axK)
    axK.set_title(f"O K-shell $\\epsilon_2$ (optical, gated) 400–600 eV; {ice} ice")
    plt.savefig(f"output/Kshell_channel_resolved_optical_{ice}.pdf", bbox_inches="tight")
    plt.show()

    # New figure: Model vs Experimental comparison (single panel, all three quantities)
    fig_comp, ax_comp = plt.subplots(figsize=(10, 7))

    # Get total model predictions at q=0 (optical limit)
    e1_dict = epsilon1_valence_E0(E, s)
    e2_dict = epsilon2_valence_E0(E, s, partitioned=use_partitioning)
    e1_model = e1_dict["total"]
    e2_model = e2_dict["total"]
    # Add K-shell contribution to epsilon2
    e2_model = e2_model + epsilon2_Kshell_E0(E, s)
    # Compute ELF with K-shell
    elf_model = elf_E0(E, s, include_kshell=True, partitioned=use_partitioning)

    # Plot all quantities on the same panel
    # Model predictions (black lines)
    ax_comp.plot(E, e2_model, 'k-', linewidth=2, label=r'Model: Im($\epsilon$)', zorder=5)
    ax_comp.plot(E, e1_model, 'k--', linewidth=2, label=r'Model: Re($\epsilon$)', zorder=5)
    ax_comp.plot(E, elf_model, 'k-.', linewidth=2, label='Model: ELF', zorder=5)

    # Experimental data (dark gray with filled shapes)
    ax_comp.plot(exp_e2_E, exp_e2, 'd-', color='darkgray', linewidth=1.5, markersize=5,
                markerfacecolor='darkgray', markeredgecolor='darkgray',
                label=r'Exp: Im($\epsilon$)', zorder=10)
    ax_comp.plot(exp_e1_E, exp_e1, 's-', color='darkgray', linewidth=1.5, markersize=5,
                markerfacecolor='darkgray', markeredgecolor='darkgray',
                label=r'Exp: Re($\epsilon$)', zorder=10)
    ax_comp.plot(exp_elf_E, exp_elf, '^-', color='darkgray', linewidth=1.5, markersize=5,
                markerfacecolor='darkgray', markeredgecolor='darkgray',
                label='Exp: ELF', zorder=10)

    ax_comp.set_xlabel("Energy Transfer ($E$; eV)")
    ax_comp.set_ylabel("Dielectric Properties")
    ax_comp.set_title(f"{ice.capitalize()} Ice")
    ax_comp.legend(ncol=2)
    ax_comp.grid(True, alpha=0.3)
    ax_comp.set_xlim(0, 30)

    fig_comp.tight_layout()
    plt.savefig(f"output/Model_vs_Experiment_{ice}.pdf", bbox_inches="tight")
    plt.show()

    plot_neff_and_I(savepath="output/Neff_I_Fig3_like.pdf", partitioned=use_partitioning)

    plot_model_vs_experiment_two_panel(savepath="output/dielectric_both.pdf", use_partitioning=use_partitioning)
