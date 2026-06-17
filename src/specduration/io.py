import numpy as np


def load_record_txt(path: str, skip_header: int = 1):
    """
    Load a two-column ground-motion record from a text file.

    Expected format
    ---------------
    The file must contain at least two columns:

        time [s]    acceleration [m/s²]

    Returns
    -------
    t : ndarray
        Time vector, in seconds.
    ug : ndarray
        Ground acceleration vector, in m/s².
    dt : float
        Time step, in seconds.
    """
    data = np.genfromtxt(path, skip_header=skip_header)

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"Unexpected file format in {path}. "
            "At least two columns are required: time and acceleration."
        )

    t = data[:, 0].astype(float)
    ug = data[:, 1].astype(float)

    if len(t) < 2:
        raise ValueError(
            "The input record has fewer than two points; "
            "the time step cannot be computed."
        )

    dt = float(t[1] - t[0])

    if not np.isfinite(dt) or dt <= 0:
        raise ValueError(
            "Invalid time step detected. The time column must be increasing."
        )

    return t, ug, dt
