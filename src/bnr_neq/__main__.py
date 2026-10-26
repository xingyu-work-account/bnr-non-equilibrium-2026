"""Command-line interface: ``python -m bnr_neq <command>``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import io, pipeline
from .config import CONFIG_PATH, load_config

NOT_IMPLEMENTED: dict[str, str] = {}
IMPLEMENTED = ("ingest", "validate", "segments", "thermo", "dynamics", "model", "states", "pulses", "report")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m bnr_neq")
    parser.add_argument("command", choices=[*IMPLEMENTED, *NOT_IMPLEMENTED, "run-all"])
    parser.add_argument("--config", type=Path, default=CONFIG_PATH)
    parser.add_argument("--data-dir", type=Path, default=io.DATA_DIR)
    parser.add_argument("--out-dir", type=Path, default=pipeline.OUTPUT_DIR)
    args = parser.parse_args(argv)

    if args.command in NOT_IMPLEMENTED:
        print(f"not implemented: {args.command} belongs to milestone {NOT_IMPLEMENTED[args.command]}",
              file=sys.stderr)
        return 2

    cfg = load_config(args.config)
    steps = IMPLEMENTED if args.command == "run-all" else (args.command,)
    for step in steps:
        if step == "ingest":
            ing = pipeline.ingest(cfg, args.data_dir)
            pipeline.write_ingest(ing, args.out_dir)
            print(f"ingest: {len(ing.dropped)} row(s) dropped; wrote {args.out_dir}/processed and tables")
        elif step == "validate":
            from . import validation

            result = validation.run(cfg, args.out_dir)
            print(f"validate: wrote {args.out_dir}/tables/sensor_validation.csv; verdicts: {result.verdicts}")
        elif step == "segments":
            from . import segments

            seg = segments.run(cfg, args.out_dir)
            n_ok = int(seg.references["pass"].sum())
            print(f"segments: wrote {args.out_dir}/tables/segments.csv; references {n_ok}/{len(seg.references)} reproduced")
        elif step == "thermo":
            from . import thermo

            th = thermo.run(cfg, args.out_dir)
            n_ok = int(th["references"]["pass"].sum())
            print(f"thermo: wrote {args.out_dir}/tables/thermo_*.csv and balance_*.csv; "
                  f"references {n_ok}/{len(th['references'])} reproduced")
        elif step == "dynamics":
            from . import dynamics

            dy = dynamics.run(cfg, args.out_dir)
            req = dy.validation[dy.validation["required"].astype(bool)]
            print(f"dynamics: estimator validation {int(req['pass'].sum())}/{len(req)} required checks passed; "
                  f"wrote {args.out_dir}/tables/dynamics.csv")
        elif step == "model":
            from . import model

            mo = model.run(cfg, args.out_dir)
            done = mo.checks[mo.checks["pass"].astype(str) != "not computed"]
            ok = bool(done["pass"].astype(bool).all())
            n_id = int(mo.identifiability["identifiable"].sum())
            print(f"model: selected {mo.selected}; consistency checks {'pass' if ok else 'FAIL'}; "
                  f"{n_id}/{len(mo.identifiability)} parameters identifiable; wrote {args.out_dir}/tables/model_fit.csv")
        elif step == "states":
            from . import states

            st = states.run(cfg, args.out_dir)
            whole = st.fractions[st.fractions["interval"] == "whole record"].iloc[0]
            print(f"states: {100 * whole['classified_fraction']:.0f} % of the record classified (coverage); idle "
                  f"(inactive, over-supplied) between {100 * whole['idle_fraction_lower_bound']:.0f} and "
                  f"{100 * whole['idle_fraction_upper_bound']:.0f} % of stirred time; "
                  f"wrote {args.out_dir}/tables/states_*.csv")
        elif step == "pulses":
            from . import pulses

            pu = pulses.run(cfg, args.out_dir)
            n_c = int((pu.metrics["status"] == "computed").sum())
            print(f"pulses: {n_c} metrics computed, {len(pu.metrics) - n_c} not computed; "
                  f"wrote {args.out_dir}/tables/pulses.csv")
        elif step == "report":
            from .report import assemble_report

            print(f"report: wrote {assemble_report(cfg, args.out_dir)}")
    if args.command == "run-all":
        print(f"run-all: ran {', '.join(IMPLEMENTED)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
