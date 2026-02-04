from pathlib import Path

from drs_duration.io import load_record_txt
from drs_duration.spectrum import compute_drs
from drs_duration.plots import plot_drs_full, drs_auto_zoom, plot_drs_zoom

from drs_duration.window import detect_window_by_derivative
from drs_duration.fit import fit_sinusoidal_in_window

# Si tu alpha_scan está como función, importa esa función.
# Si lo tienes como script, dime y lo acomodamos.
from drs_duration.alpha_scan import scan_alpha  # <- debe existir


def main():
    # 1) Elegir archivo (usa el primero .txt en examples)
    path = str(next(Path("examples").glob("*.txt")))
    print("Usando archivo:", path)

    # 2) Cargar registro
    t, ug, dt = load_record_txt(path)

    # 3) Calcular DRS completo
    T_values, RESP = compute_drs(ug, dt, zeta=0.05, Tmin=None, Tmax=15.0, dT=0.01)

    # 4) Plot DRS completo
    plot_drs_full(T_values, RESP)

    # 5) Auto-zoom + prints + plot zoom
    T_fit, RESP_fit, info = drs_auto_zoom(T_values, RESP, zoom_frac=0.20, min_points=15)

    print(f"Tp detectado automáticamente = {info['Tp']:.6f} s, Dmax = {info['Dmax']:.6e}")
    print(
        f"Zoom automático final: [{info['Tmin_zoom']:.6f}, {info['Tmax_zoom']:.6f}] s — "
        f"{info['n_points']} puntos usados. (±{info['zoom_used']*100:.0f}%)"
    )

    plot_drs_zoom(T_fit, RESP_fit, info["Tp"], info["Dmax"])

    # =========================================================
    # 6) Barrido de alpha para elegir alpha óptimo (sin C)
    # =========================================================
    print("\n--- Barrido de α para escoger α óptimo ---")
    df_resultados, mejor_alpha_sinC = scan_alpha(
        T_fit,
        RESP_fit,
        alphas=None,     # si None, tu función debe usar 0.01..0.60
        plot=True        # si quieres el gráfico R2 vs alpha
    )

    print(f"\n🎯 α óptimo (sin C) = {mejor_alpha_sinC:.3f}")

    # =========================================================
    # 7) Detección final de ventana con α óptimo + plots PRO
    # =========================================================
    print("\n--- Detección final de ventana ---")
    det_res = detect_window_by_derivative(
        T_fit, RESP_fit,
        alpha=mejor_alpha_sinC,
        plot=True
    )

    # =========================================================
    # 8) Ajuste senoidal final (sin C y con C) + plot
    # =========================================================
    print("\n--- Ajuste senoidal final ---")
    res_fit = fit_sinusoidal_in_window(
        T_fit,
        RESP_fit,
        detect_res=det_res,
        use_detect_function=False,
        plot=True,
        verbose=True
    )

    noC = res_fit.get("result_noC")
    conC = res_fit.get("result_C")

    print("\n=== RESUMEN FINAL ===")
    if noC:
        A, td, phi = noC["popt"]
        print(f"[SIN C] A={A:.6f}, td={td:.6f}, phi={phi:.6f}, R2={noC['R2']:.6f}, cycles={noC['cycles']:.6f}")
    else:
        print("[SIN C] No convergió")

    if conC:
        A, td, phi, C = conC["popt"]
        print(f"[CON C] A={A:.6f}, td={td:.6f}, phi={phi:.6f}, C={C:.6f}, R2={conC['R2']:.6f}, cycles={conC['cycles']:.6f}")
    else:
        print("[CON C] No convergió")


if __name__ == "__main__":
    main()



