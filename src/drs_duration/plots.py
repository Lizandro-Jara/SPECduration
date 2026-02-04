import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# FIG 1: DRS COMPLETO (PRO) -> show + save
# ============================================================

def plot_drs_full(T_values, RESP, title="Espectro de Respuesta de Desplazamiento (DRS)"):
    """Versión para mostrar en pantalla (ejemplo)."""
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    plt.figure(figsize=(10, 5))
    plt.plot(T_values, RESP, color="navy", lw=2)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Periodo Tn [s]", fontsize=12)
    plt.ylabel("Desplazamiento relativo pico |u| [m]", fontsize=12)
    plt.title(title, fontsize=14, fontweight="bold")

    plt.xlim(float(np.min(T_values)), float(np.max(T_values)))
    plt.ylim(0, float(np.nanmax(RESP)) * 1.10)
    plt.xticks(np.arange(0, np.ceil(np.max(T_values)) + 1, 1))

    plt.tight_layout()
    plt.show()


def save_drs_full_png(T_values, RESP, outpath, title="Espectro de Respuesta de Desplazamiento (DRS)", dpi=300):
    """Versión para guardar PNG (pipeline)."""
    T_values = np.asarray(T_values, dtype=float)
    RESP = np.asarray(RESP, dtype=float)

    mask = np.isfinite(T_values) & np.isfinite(RESP)
    T_values = T_values[mask]
    RESP = RESP[mask]

    plt.figure(figsize=(10, 5))
    plt.plot(T_values, RESP, color="navy", lw=2)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Periodo Tn [s]", fontsize=12)
    plt.ylabel("Desplazamiento relativo pico |u| [m]", fontsize=12)
    plt.title(title, fontsize=14, fontweight="bold")

    plt.xlim(float(np.min(T_values)), float(np.max(T_values)))
    plt.ylim(0, float(np.nanmax(RESP)) * 1.10)
    plt.xticks(np.arange(0, np.ceil(np.max(T_values)) + 1, 1))

    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ============================================================
# FIG 2: AUTO-ZOOM + DRS ZOOM (PRO) -> show + save
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

    # Si zoom queda muy angosto, ampliar
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
    """Versión para mostrar en pantalla (ejemplo)."""
    T_fit = np.asarray(T_fit, dtype=float)
    RESP_fit = np.asarray(RESP_fit, dtype=float)

    plt.figure(figsize=(10, 5))
    plt.plot(T_fit, RESP_fit, color="navy", lw=2, label="DRS (Zoom Vent. Pred)")
    plt.scatter([Tp], [Dmax], color="red", s=60, zorder=5,
                label=f"Tp = {Tp:.3f} s\nDmax = {Dmax:.3f} m")
    plt.axvline(Tp, color="black", ls="--", lw=1)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Periodo Tn [s]", fontsize=12)
    plt.ylabel("Desplazamiento relativo pico |u| [m]", fontsize=12)

    if title is None:
        title = f"DRS (Zoom Ventana Predominante) — Tp = {Tp:.3f} s"
    plt.title(title, fontsize=14, fontweight="bold")

    if RESP_fit.size > 0:
        plt.ylim(float(np.min(RESP_fit)) * 0.98, float(np.max(RESP_fit)) * 1.02)

    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.legend(fontsize=11, loc="best")

    plt.tight_layout()
    plt.show()


def save_drs_zoom_png(T_fit, RESP_fit, Tp, Dmax, outpath, title=None, dpi=300):
    """Versión para guardar PNG (pipeline)."""
    T_fit = np.asarray(T_fit, dtype=float)
    RESP_fit = np.asarray(RESP_fit, dtype=float)

    plt.figure(figsize=(10, 5))
    plt.plot(T_fit, RESP_fit, color="navy", lw=2, label="DRS (Zoom Vent. Pred)")
    plt.scatter([Tp], [Dmax], color="red", s=60, zorder=5,
                label=f"Tp = {Tp:.3f} s\nDmax = {Dmax:.3f} m")
    plt.axvline(Tp, color="black", ls="--", lw=1)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Periodo Tn [s]", fontsize=12)
    plt.ylabel("Desplazamiento relativo pico |u| [m]", fontsize=12)

    if title is None:
        title = f"DRS (Zoom Ventana Predominante) — Tp = {Tp:.3f} s"
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
# FIG 3: R2 vs ALPHA (PRO) -> save
# ============================================================

def save_r2_vs_alpha_png(df_alpha, best_alpha, outpath, dpi=300):
    df_alpha = df_alpha.copy()

    alpha_vals = np.asarray(df_alpha["alpha"], dtype=float)
    R2_sin = np.asarray(df_alpha["R2_noC"], dtype=float)
    R2_con = np.asarray(df_alpha["R2_C"], dtype=float)

    plt.figure(figsize=(11, 5))

    plt.plot(alpha_vals, R2_sin, color="red", lw=2.0, marker="o", ms=5,
             label="R² sin C (modelo físico)")
    plt.plot(alpha_vals, R2_con, color="green", lw=1.5, ls="--", marker="o", ms=4,
             label="R² con C (modelo extendido)")

    if np.isfinite(best_alpha):
        plt.axvline(float(best_alpha), color="black", ls="--", lw=1.5,
                    label=f"α óptimo (sin C) = {float(best_alpha):.3f}")

    plt.title("Variación del coeficiente de determinación R² con α",
              fontsize=15, fontweight="bold")
    plt.xlabel("α", fontsize=13)
    plt.ylabel("R²", fontsize=13)

    plt.grid(True, ls="--", lw=0.5, alpha=0.65)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ============================================================
# FIG 4: AJUSTE SENOIDAL FINAL (PRO + curva visual PCHIP) -> save
# ============================================================

def save_final_sine_fit_png(T_win, R_win, result_noC, result_C, outpath, dpi=300):
    T_win = np.asarray(T_win, dtype=float)
    R_win = np.asarray(R_win, dtype=float)

    mask = np.isfinite(T_win) & np.isfinite(R_win)
    T_win = T_win[mask]
    R_win = R_win[mask]

    plt.figure(figsize=(10.5, 6.8))

    if T_win.size < 2:
        plt.scatter(T_win, R_win, s=28, color="royalblue", label="DRS (ventana)")
        plt.title("Ajustes senoidales dentro de la ventana detectada", fontsize=18, fontweight="bold")
        plt.xlabel("Periodo T [s]", fontsize=14)
        plt.ylabel("DRS", fontsize=14)
        plt.grid(True, ls="--", lw=0.6, alpha=0.65)
        plt.legend(fontsize=12, loc="best")
        plt.tight_layout()
        plt.savefig(outpath, dpi=dpi)
        plt.close()
        return

    # Ordenar
    order = np.argsort(T_win)
    T_win = T_win[order]
    R_win = R_win[order]

    # Quitar duplicados en T
    T_unique, idx_unique = np.unique(T_win, return_index=True)
    R_unique = R_win[idx_unique]

    # Curva visual suave
    T_dense = np.linspace(float(T_unique.min()), float(T_unique.max()), 900)

    try:
        from scipy.interpolate import PchipInterpolator
        pchip = PchipInterpolator(T_unique, R_unique, extrapolate=False)
        R_smooth = pchip(T_dense)
    except Exception:
        R_smooth = np.interp(T_dense, T_unique, R_unique)

    # DRS puntos + curva visual
    
    plt.plot(T_dense, R_smooth, color="royalblue", lw=2.4, alpha=0.85,
             label="DRS (curva visual)")
    
    # Ajustes (usamos sine_fit / sine_fit_C del dict de resultados)
    # result_noC["popt"] = [A, td, phi]
    # result_C["popt"]   = [A, td, phi, C]
    from .fit import sine_fit, sine_fit_C  # import local para evitar ciclos

    if result_noC is not None:
        plt.plot(T_dense, sine_fit(T_dense, *result_noC["popt"]),
                 color="crimson", lw=2.2,
                 label=f"Sin C (R²={result_noC['R2']:.3f})")

    if result_C is not None:
        plt.plot(T_dense, sine_fit_C(T_dense, *result_C["popt"]),
                 color="green", lw=2.0, ls="--",
                 label=f"Con C (R²={result_C['R2']:.3f})")

    plt.title("Ajustes senoidales dentro de la ventana detectada", fontsize=18, fontweight="bold")
    plt.xlabel("Periodo T [s]", fontsize=14)
    plt.ylabel("DRS", fontsize=14)
    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()
