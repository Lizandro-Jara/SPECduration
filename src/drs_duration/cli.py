import argparse
from .pipeline import run_pipeline, PipelineConfig


def main():
    p = argparse.ArgumentParser(prog="drs-duration")
    p.add_argument("--input", required=True, help="Ruta del acelerograma .txt")
    p.add_argument("--out", required=True, help="Carpeta de salida (outputs)")
    p.add_argument("--no-plots", action="store_true", help="No generar figuras")

    args = p.parse_args()

    cfg = PipelineConfig()
    result = run_pipeline(args.input, args.out, cfg=cfg, make_plots=(not args.no_plots))

    print("\n=== RESULTADOS FINALES ===")
    for k, v in result.items():
        print(f"{k}: {v}")
