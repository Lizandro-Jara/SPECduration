import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# FIG 1: FULL DRS -> show + save
# ============================================================

def plot_drs_full(T_values, RESP, title="Displacement Response Spectrum (DRS)"):
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    plt.figure(figsize=(10, 5))
    plt.plot(T_values, RESP, color="navy", lw=2)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period Tn [s]", fontsize=12)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=12)
    plt.title(title, fontsize=14, fontweight="bold")

    plt.xlim(float(np.min(T_values)), float(np.max(T_values)))
    plt.ylim(0, float(np.nanmax(RESP)) * 1.10)
    plt.xticks(np.arange(0, np.ceil(np.max(T_values)) + 1, 1))

    plt.tight_layout()
    plt.show()


def save_drs_full_png(T_values, RESP, outpath, title="Displacement Response Spectrum (DRS)", dpi=300):
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    plt.figure(figsize=(10, 5))
    plt.plot(T_values, RESP, color="navy", lw=2)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period Tn [s]", fontsize=12)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=12)
    plt.title(title, fontsize=14, fontweight="bold")

    plt.xlim(float(np.min(T_values)), float(np.max(T_values)))
    plt.ylim(0, float(np.nanmax(RESP)) * 1.10)
    plt.xticks(np.arange(0, np.ceil(np.max(T_values)) + 1, 1))

    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ============================================================
# FIG 2: AUTO-ZOOM + DRS ZOOM -> show + save
# ============================================================

def drs_auto_zoom(T_values, RESP, zoom_frac=0.20, min_points=15):
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    idx_peak = int(np.nanargmax(RESP))
    Tp = float(T_values[idx_peak])
    Dmax = float(RESP[idx_peak])

    Tmin_full = float(T_values[0])
    Tmax_full = float(T_values[-1])

    Tmin_zoom = max(Tmin_full, Tp * (1 - zoom_frac))
    Tmax_zoom = min(Tmax_full, Tp * (1 + zoom_frac))

    mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)
    T_fit = T_values[mask_zoom]
    RESP_fit = RESP[mask_zoom]

    zoom_used = float(zoom_frac)

    if T_fit.size < min_points:
        zoom_used = 0.35
        Tmin_zoom = max(Tmin_full, Tp * (1 - zoom_used))
        Tmax_zoom = min(Tmax_full, Tp * (1 + zoom_used))
        mask_zoom = (T_values >= Tmin_zoom) & (T_values <= Tmax_zoom)
        T_fit = T_values[mask_zoom]
        RESP_fit = RESP[mask_zoom]

    info = {
        "Tp": Tp,
        "Dmax": Dmax,
        "Tmin_zoom": float(Tmin_zoom),
        "Tmax_zoom": float(Tmax_zoom),
        "n_points": int(T_fit.size),
        "zoom_used": float(zoom_used),
    }
    return T_fit, RESP_fit, info


def plot_drs_zoom(T_fit, RESP_fit, Tp, Dmax, title=None):
    T_fit = np.asarray(T_fit, dtype=float)
    RESP_fit = np.asarray(RESP_fit, dtype=float)

    plt.figure(figsize=(10, 5))
    plt.plot(T_fit, RESP_fit, color="navy", lw=2, label="DRS (predominant-window zoom)")
    plt.scatter([Tp], [Dmax], color="red", s=60, zorder=5,
                label=f"Tp = {Tp:.3f} s\nDmax = {Dmax:.3f} m")
    plt.axvline(Tp, color="black", ls="--", lw=1)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period Tn [s]", fontsize=12)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=12)

    if title is None:
        title = f"DRS (Predominant Window Zoom) — Tp = {Tp:.3f} s"
    plt.title(title, fontsize=14, fontweight="bold")

    if RESP_fit.size > 0:
        plt.ylim(float(np.min(RESP_fit)) * 0.98, float(np.max(RESP_fit)) * 1.02)

    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.legend(fontsize=11, loc="best")

    plt.tight_layout()
    plt.show()


def save_drs_zoom_png(T_fit, RESP_fit, Tp, Dmax, outpath, title=None, dpi=300):
    T_fit = np.asarray(T_fit, dtype=float)
    RESP_fit = np.asarray(RESP_fit, dtype=float)

    plt.figure(figsize=(10, 5))
    plt.plot(T_fit, RESP_fit, color="navy", lw=2, label="DRS (predominant-window zoom)")
    plt.scatter([Tp], [Dmax], color="red", s=60, zorder=5,
                label=f"Tp = {Tp:.3f} s\nDmax = {Dmax:.3f} m")
    plt.axvline(Tp, color="black", ls="--", lw=1)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period Tn [s]", fontsize=12)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=12)

    if title is None:
        title = f"DRS (Predominant Window Zoom) — Tp = {Tp:.3f} s"
    plt.title(title, fontsize=14, fontweight="bold")

    if RESP_fit.size > 0:
        plt.ylim(float(np.min(RESP_fit)) * 0.98, float(np.max(RESP_fit)) * 1.02)

    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.legend(fontsize=11, loc="best")

    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ============================================================
# FIG 3: R2 vs ALPHA -> save
# ============================================================

def save_r2_vs_alpha_png(df_alpha, best_alpha, outpath, dpi=300):
    df_alpha = df_alpha.copy()

    alpha_vals = np.asarray(df_alpha["alpha"], dtype=float)
    R2_vals = np.asarray(df_alpha["R2"], dtype=float)

    plt.figure(figsize=(11, 5))

    plt.plot(
        alpha_vals,
        R2_vals,
        color="crimson",
        lw=2.0,
        marker="o",
        ms=5,
        label="R² of sinusoidal fit"
    )

    if np.isfinite(best_alpha):
        plt.axvline(
            float(best_alpha),
            color="black",
            ls="--",
            lw=1.5,
            label=f"Optimal α = {float(best_alpha):.3f}"
        )

    plt.title(
        "Variation of the coefficient of determination R² with α",
        fontsize=15,
        fontweight="bold"
    )
    plt.xlabel("α", fontsize=13)
    plt.ylabel("R²", fontsize=13)

    plt.grid(True, ls="--", lw=0.5, alpha=0.65)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ============================================================
# FIG 4: FINAL SINUSOIDAL FIT -> save
# ============================================================

def save_final_sine_fit_png(T_win, R_win, fit_result, outpath, dpi=300):
    T_win = np.asarray(T_win, dtype=float)
    R_win = np.asarray(R_win, dtype=float)

    mask = np.isfinite(T_win) & np.isfinite(R_win)
    T_win = T_win[mask]
    R_win = R_win[mask]

    plt.figure(figsize=(10.5, 6.8))

    if T_win.size < 2:
        plt.scatter(T_win, R_win, s=28, color="royalblue", label="DRS (window)")
        plt.title("Sinusoidal fit within the detected window", fontsize=18, fontweight="bold")
        plt.xlabel("Period T [s]", fontsize=14)
        plt.ylabel("DRS", fontsize=14)
        plt.grid(True, ls="--", lw=0.6, alpha=0.65)
        plt.legend(fontsize=12, loc="best")
        plt.tight_layout()
        plt.savefig(outpath, dpi=dpi)
        plt.close()
        return

    # sort
    order = np.argsort(T_win)
    T_win = T_win[order]
    R_win = R_win[order]

    # remove duplicated T values
    T_unique, idx_unique = np.unique(T_win, return_index=True)
    R_unique = R_win[idx_unique]

    # smooth visual curve
    T_dense = np.linspace(float(T_unique.min()), float(T_unique.max()), 900)

    try:
        from scipy.interpolate import PchipInterpolator
        pchip = PchipInterpolator(T_unique, R_unique, extrapolate=False)
        R_smooth = pchip(T_dense)
    except Exception:
        R_smooth = np.interp(T_dense, T_unique, R_unique)

    # DRS visual
    plt.plot(
        T_dense,
        R_smooth,
        color="royalblue",
        lw=2.4,
        alpha=0.85,
        label="DRS (visual curve)"
    )

    from .fit import sine_fit  # local import to avoid cycles

    if fit_result is not None:
        plt.plot(
            T_dense,
            sine_fit(T_dense, *fit_result["popt"]),
            color="crimson",
            lw=2.2,
            label=f"Sinusoidal fit (R²={fit_result['R2']:.3f})"
        )

    plt.title("Sinusoidal fit within the detected window", fontsize=18, fontweight="bold")
    plt.xlabel("Period T [s]", fontsize=14)
    plt.ylabel("DRS", fontsize=14)
    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()
