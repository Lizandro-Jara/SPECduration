import numpy as np

def load_record_txt(path: str, skip_header: int = 1):
    """
    Lee un registro tipo:
    t(s)   ag(m/s^2)
    0.00   -0.0004219
    0.01   -0.000437
    ...

    Retorna: t, ug, dt
    """
    A = np.genfromtxt(path, skip_header=skip_header)

    if A.ndim != 2 or A.shape[1] < 2:
        raise ValueError(f"Formato inesperado en {path}. Se esperaba al menos 2 columnas (t, ug).")

    t = A[:, 0].astype(float)
    ug = A[:, 1].astype(float)

    if len(t) < 2:
        raise ValueError("El registro tiene menos de 2 puntos; no se puede calcular dt.")

    dt = float(t[1] - t[0])
    return t, ug, dt

