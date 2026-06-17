import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------
# Basic plotting helpers
# ----------------------------------------------------------------------

def _clean_xy(x, y):
    """
    Remove non-finite values and sort paired x-y data by the x coordinate.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if x.size == 0:
        raise RuntimeError("No valid data available for plotting.")

    order = np.argsort(x)
    return x[order], y[order]


def _safe_ylim(y, margin=0.08):
    """
    Return a stable y-axis range with a small visual margin.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]

    if y.size == 0:
        return None

    ymin = float(np.nanmin(y))
    ymax = float(np.nanmax(y))

    if np.isclose(ymin, ymax):
        delta = abs(ymax) * 0.10 + 1e-12
        return ymin - delta, ymax + delta

    delta = (ymax - ymin) * margin
    return ymin - delta, ymax + delta


def _get_window_limits(T_values, T_left, T_right, left_frac=0.25, right_frac=0.25):
    """
    Define a local period range around the detected spectral window.
    """
    T_values = np.asarray(T_values, dtype=float)

    width = max(float(T_right) - float(T_left), 1e-12)

    x_min = max(float(np.min(T_values)), float(T_left) - left_frac * width)
    x_max = min(float(np.max(T_values)), float(T_right) + right_frac * width)

    if x_max <= x_min:
        x_min = float(np.min(T_values))
        x_max = float(np.max(T_values))

    return x_min, x_max


def _get_second_derivative_array(detect_result):
    """
    Return the second spectral derivative.

    The preferred key is d2S. The d2R key is kept only for backward
    compatibility with older output dictionaries.
    """
    if "d2S" in detect_result:
        return np.asarray(detect_result["d2S"], dtype=float)

    if "d2R" in detect_result:
        return np.asarray(detect_result["d2R"], dtype=float)

    raise KeyError("detect_result does not contain 'd2S' or 'd2R'.")


def _get_model_function():
    """
    Return the harmonic spectral model used for fitting.
    """
    try:
        from .fit import spectral_sine_model
        return spectral_sine_model
    except Exception:
        from .fit import sine_fit
        return sine_fit


# ----------------------------------------------------------------------
# Full displacement response spectrum
# ----------------------------------------------------------------------

def save_drs_full_png(
    T_values,
    RESP,
    outpath,
    title="Displacement Response Spectrum (DRS)",
    dpi=300,
):
    """
    Save the full displacement response spectrum.
    """
    T_values, RESP = _clean_xy(T_values, RESP)

    plt.figure(figsize=(10.5, 5.8))

    plt.plot(
        T_values,
        RESP,
        lw=2.2,
        label="DRS",
    )

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period T [s]", fontsize=13)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=13)
    plt.title(title, fontsize=17, fontweight="bold")

    plt.xlim(float(np.min(T_values)), float(np.max(T_values)))

    ylim = _safe_ylim(RESP)
    if ylim is not None:
        plt.ylim(max(0.0, ylim[0]), ylim[1])

    plt.legend(fontsize=12, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ----------------------------------------------------------------------
# Detected DRS window with harmonic fit
# ----------------------------------------------------------------------

def save_drs_window_png(
    T_values,
    RESP,
    detect_result,
    outpath,
    fit_result=None,
    title="DRS window detection and harmonic fit",
    dpi=300,
):
    """
    Save the main local DRS figure.

    The blue curve shows the original DRS in a local range around the
    detected window. The orange curve shows the harmonic fit evaluated
    inside the detected window. The vertical lines indicate T_left,
    T_peak, and T_right.
    """
    T_values, RESP = _clean_xy(T_values, RESP)

    T_peak = float(detect_result["T_peak"])
    T_left = float(detect_result["T_left"])
    T_right = float(detect_result["T_right"])

    idx_peak = int(np.argmin(np.abs(T_values - T_peak)))
    S_peak = float(RESP[idx_peak])

    mask_window = (T_values >= T_left) & (T_values <= T_right)

    x_min, x_max = _get_window_limits(
        T_values=T_values,
        T_left=T_left,
        T_right=T_right,
        left_frac=0.25,
        right_frac=0.25,
    )

    local_mask = (T_values >= x_min) & (T_values <= x_max)

    plt.figure(figsize=(10.8, 5.9))

    # Original DRS in the local plotting range.
    plt.plot(
        T_values[local_mask],
        RESP[local_mask],
        lw=2.2,
        label="DRS original",
    )

    # Harmonic fit inside the detected window.
    if fit_result is not None and "popt" in fit_result:
        model_fun = _get_model_function()

        T_win = T_values[mask_window]

        if T_win.size >= 2:
            T_dense = np.linspace(
                float(np.min(T_win)),
                float(np.max(T_win)),
                900,
            )

            S_fit = model_fun(T_dense, *fit_result["popt"])

            if "R2" in fit_result:
                fit_label = f"Harmonic fit (R² = {fit_result['R2']:.3f})"
            else:
                fit_label = "Harmonic fit"

            plt.plot(
                T_dense,
                S_fit,
                lw=2.0,
                label=fit_label,
            )

    # Peak period marker.
    plt.scatter(
        [T_peak],
        [S_peak],
        s=80,
        zorder=5,
        label=f"T_peak = {T_peak:.3f} s",
    )

    # Window limits and peak period.
    plt.axvline(
        T_left,
        ls="--",
        lw=1.6,
        label=f"T_left = {T_left:.3f} s",
    )

    plt.axvline(
        T_right,
        ls="--",
        lw=1.6,
        label=f"T_right = {T_right:.3f} s",
    )

    plt.axvline(
        T_peak,
        ls=":",
        lw=1.8,
    )

    td_text = ""
    if fit_result is not None and "td" in fit_result:
        td_text = f" | td = {fit_result['td']:.3f} s"

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period T [s]", fontsize=13)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=13)
    plt.title(f"{title}{td_text}", fontsize=17, fontweight="bold")

    plt.xlim(x_min, x_max)

    ylim = _safe_ylim(RESP[local_mask])
    if ylim is not None:
        plt.ylim(max(0.0, ylim[0]), ylim[1])

    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ----------------------------------------------------------------------
# Separate harmonic fit figure
# ----------------------------------------------------------------------

def save_drs_fit_overlay_png(
    T_values,
    RESP,
    detect_result,
    fit_result,
    outpath,
    title="Harmonic fitting inside detected window",
    dpi=300,
):
    """
    Save the harmonic fit as a separate figure.

    This function is kept for backward compatibility. The recommended
    report figure is save_drs_window_png, which already combines the
    detected DRS window and the harmonic fit.
    """
    T_values, RESP = _clean_xy(T_values, RESP)

    T_left = float(detect_result["T_left"])
    T_right = float(detect_result["T_right"])

    mask_window = (T_values >= T_left) & (T_values <= T_right)

    T_win = T_values[mask_window]
    S_win = RESP[mask_window]

    if T_win.size < 2:
        raise RuntimeError(
            "Not enough points inside the detected window to plot the fit."
        )

    plt.figure(figsize=(10.8, 5.9))

    plt.plot(
        T_win,
        S_win,
        lw=2.2,
        label="DRS window",
    )

    if fit_result is not None and "popt" in fit_result:
        model_fun = _get_model_function()

        T_dense = np.linspace(
            float(np.min(T_win)),
            float(np.max(T_win)),
            900,
        )

        S_fit = model_fun(T_dense, *fit_result["popt"])

        plt.plot(
            T_dense,
            S_fit,
            lw=2.0,
            label=f"Harmonic fit (R² = {fit_result['R2']:.3f})",
        )

    td_text = ""
    if fit_result is not None and "td" in fit_result:
        td_text = f" | td = {fit_result['td']:.3f} s"

    plt.title(
        f"{title}{td_text}",
        fontsize=17,
        fontweight="bold",
    )

    plt.xlabel("Period T [s]", fontsize=13)
    plt.ylabel("DRS", fontsize=13)

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)

    ylim = _safe_ylim(S_win)
    if ylim is not None:
        plt.ylim(max(0.0, ylim[0]), ylim[1])

    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ----------------------------------------------------------------------
# Second spectral derivative detector
# ----------------------------------------------------------------------

def save_second_derivative_window_png(
    detect_result,
    outpath,
    title="Second spectral derivative detector",
    dpi=300,
):
    """
    Save the S''(T) diagnostic figure around the detected window.

    The full local S''(T) curve is shown together with the portion
    inside the detected spectral window. The zero line is included to
    verify the sign changes used by the detector.
    """
    T = np.asarray(detect_result["T_u"], dtype=float)
    S2 = _get_second_derivative_array(detect_result)

    T, S2 = _clean_xy(T, S2)

    T_left = float(detect_result["T_left"])
    T_right = float(detect_result["T_right"])
    T_peak = float(detect_result["T_peak"])

    width = max(T_right - T_left, 1e-12)

    x_min = max(float(np.min(T)), T_left - 0.35 * width)
    x_max = min(float(np.max(T)), T_right + 0.35 * width)

    local_mask = (T >= x_min) & (T <= x_max)
    window_mask = (T >= T_left) & (T <= T_right)

    T_local = T[local_mask]
    S2_local = S2[local_mask]

    T_win = T[window_mask]
    S2_win = S2[window_mask]

    plt.figure(figsize=(10.8, 5.8))

    # Local S''(T) curve.
    plt.plot(
        T_local,
        S2_local,
        lw=2.0,
        label=r"$S''(T)$",
    )

    # S''(T) values inside the detected window.
    if T_win.size > 0:
        plt.plot(
            T_win,
            S2_win,
            lw=2.2,
            label=r"$S''(T)$ inside detected window",
        )

    # Zero-curvature reference.
    plt.axhline(
        0.0,
        lw=1.2,
        label=r"$S''(T)=0$",
    )

    # Window limits and peak period.
    plt.axvline(
        T_left,
        ls="--",
        lw=1.6,
        label=f"T_left = {T_left:.3f} s",
    )

    plt.axvline(
        T_peak,
        ls=":",
        lw=1.8,
        label=f"T_peak = {T_peak:.3f} s",
    )

    plt.axvline(
        T_right,
        ls="--",
        lw=1.6,
        label=f"T_right = {T_right:.3f} s",
    )

    # Controlled y-axis scaling for readable diagnostic plots.
    S2_local_finite = S2_local[np.isfinite(S2_local)]

    if S2_local_finite.size > 0:
        y_low = float(np.nanpercentile(S2_local_finite, 1))
        y_high = float(np.nanpercentile(S2_local_finite, 99))

        if T_win.size > 0:
            S2_win_finite = S2_win[np.isfinite(S2_win)]

            if S2_win_finite.size > 0:
                y_low = min(y_low, float(np.nanmin(S2_win_finite)))
                y_high = max(y_high, float(np.nanmax(S2_win_finite)))

        y_low = min(y_low, 0.0)
        y_high = max(y_high, 0.0)

        min_range = 0.50
        current_range = y_high - y_low

        if current_range < min_range:
            center = 0.5 * (y_low + y_high)
            y_low = center - min_range / 2.0
            y_high = center + min_range / 2.0

            y_low = min(y_low, 0.0)
            y_high = max(y_high, 0.0)

        margin = 0.12 * (y_high - y_low)
        y_low_plot = y_low - margin
        y_high_plot = y_high + margin

        plt.ylim(y_low_plot, y_high_plot)

        real_min = float(np.nanmin(S2_local_finite))
        real_max = float(np.nanmax(S2_local_finite))

        if real_min < y_low_plot or real_max > y_high_plot:
            plt.text(
                0.02,
                0.03,
                "Y-axis clipped for visualization",
                transform=plt.gca().transAxes,
                fontsize=9,
                alpha=0.70,
            )

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period T [s]", fontsize=13)
    plt.ylabel(r"$S''(T)$", fontsize=15)
    plt.title(title, fontsize=17, fontweight="bold")

    plt.xlim(x_min, x_max)

    plt.legend(fontsize=10, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


# ----------------------------------------------------------------------
# Backward-compatible plotting aliases
# ----------------------------------------------------------------------

def save_spectrum_fit_overlay_png(
    T_values,
    RESP,
    detect_result,
    fit_result,
    outpath,
    title="Harmonic fitting inside detected window",
    dpi=300,
):
    """
    Backward-compatible alias for the separate harmonic fit figure.
    """
    save_drs_fit_overlay_png(
        T_values=T_values,
        RESP=RESP,
        detect_result=detect_result,
        fit_result=fit_result,
        outpath=outpath,
        title=title,
        dpi=dpi,
    )


def save_drs_window_with_fit_inset_png(
    T_values,
    RESP,
    detect_result,
    fit_result,
    outpath,
    title="DRS window detection and harmonic fit",
    dpi=300,
):
    """
    Backward-compatible alias for the combined DRS window and fit figure.
    """
    save_drs_window_png(
        T_values=T_values,
        RESP=RESP,
        detect_result=detect_result,
        fit_result=fit_result,
        outpath=outpath,
        title=title,
        dpi=dpi,
    )


def save_drs_zoom_png(
    T_zoom,
    RESP_zoom,
    T_peak,
    Dmax,
    outpath,
    title=None,
    dpi=300,
):
    """
    Save a simple DRS zoom around the peak period.

    This function is kept for older workflows and quick internal checks.
    """
    T_zoom, RESP_zoom = _clean_xy(T_zoom, RESP_zoom)

    T_peak = float(T_peak)
    Dmax = float(Dmax)

    plt.figure(figsize=(10, 5))

    plt.plot(
        T_zoom,
        RESP_zoom,
        lw=2,
        label="DRS zoom",
    )

    plt.scatter(
        [T_peak],
        [Dmax],
        s=65,
        zorder=5,
        label=f"T_peak = {T_peak:.3f} s",
    )

    plt.axvline(
        T_peak,
        ls="--",
        lw=1.2,
    )

    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.xlabel("Period T [s]", fontsize=12)
    plt.ylabel("Peak relative displacement |u| [m]", fontsize=12)

    if title is None:
        title = f"DRS zoom around T_peak = {T_peak:.3f} s"

    plt.title(title, fontsize=14, fontweight="bold")

    ylim = _safe_ylim(RESP_zoom)
    if ylim is not None:
        plt.ylim(max(0.0, ylim[0]), ylim[1])

    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()


def save_final_sine_fit_png(
    T_win,
    R_win,
    fit_result,
    outpath,
    dpi=300,
):
    """
    Save the legacy final sine-fit figure.

    This function remains available for older workflows. New reports
    should use save_drs_window_png, which shows the harmonic fit together
    with the detected DRS window.
    """
    T_win, R_win = _clean_xy(T_win, R_win)

    plt.figure(figsize=(10.5, 6.2))

    if T_win.size < 2:
        plt.scatter(
            T_win,
            R_win,
            s=30,
            label="DRS window",
        )

        plt.title(
            "Harmonic fitting inside detected window",
            fontsize=15,
            fontweight="bold",
        )

        plt.xlabel("Period T [s]", fontsize=12)
        plt.ylabel("DRS", fontsize=12)
        plt.grid(True, ls="--", lw=0.6, alpha=0.65)
        plt.legend(fontsize=11, loc="best")
        plt.tight_layout()
        plt.savefig(outpath, dpi=dpi)
        plt.close()
        return

    order = np.argsort(T_win)
    T_win = T_win[order]
    R_win = R_win[order]

    T_unique, idx_unique = np.unique(T_win, return_index=True)
    R_unique = R_win[idx_unique]

    T_dense = np.linspace(
        float(T_unique.min()),
        float(T_unique.max()),
        900,
    )

    try:
        from scipy.interpolate import PchipInterpolator

        pchip = PchipInterpolator(
            T_unique,
            R_unique,
            extrapolate=False,
        )
        R_smooth = pchip(T_dense)

    except Exception:
        R_smooth = np.interp(T_dense, T_unique, R_unique)

    plt.plot(
        T_dense,
        R_smooth,
        lw=2.2,
        alpha=0.85,
        label="DRS window",
    )

    if fit_result is not None:
        model_fun = _get_model_function()

        plt.plot(
            T_dense,
            model_fun(T_dense, *fit_result["popt"]),
            lw=2.0,
            label=f"Harmonic fit (R² = {fit_result['R2']:.3f})",
        )

    td_text = ""
    if fit_result is not None and "td" in fit_result:
        td_text = f" | td = {fit_result['td']:.3f} s"

    plt.title(
        f"Harmonic fitting inside detected window{td_text}",
        fontsize=15,
        fontweight="bold",
    )

    plt.xlabel("Period T [s]", fontsize=12)
    plt.ylabel("DRS", fontsize=12)
    plt.grid(True, ls="--", lw=0.6, alpha=0.65)
    plt.legend(fontsize=11, loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.close()