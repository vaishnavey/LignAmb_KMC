from __future__ import annotations

import argparse
from pathlib import Path

from .core import CHAIN_COMPATIBLE_MODES, DIMER_ONLY_MODES, MIXED_MODE, SUPPORTED_MODES, TopologyRequest, write_mixed_topology_ensemble, write_topology_bundle, write_topology_ensemble


def resolve_modes(mode_spec: str, chain_length: int) -> list[str]:
    normalized = mode_spec.strip().lower()
    if not normalized:
        raise ValueError("At least one topology mode is required.")
    if normalized == "all":
        modes = list(CHAIN_COMPATIBLE_MODES)
        if chain_length == 2:
            modes.extend(DIMER_ONLY_MODES)
        return modes
    if normalized == "all-supported":
        modes = list(CHAIN_COMPATIBLE_MODES)
        if chain_length == 2:
            modes.extend(DIMER_ONLY_MODES)
        return modes

    modes = [mode.strip() for mode in mode_spec.split(",") if mode.strip()]
    invalid = [mode for mode in modes if mode.lower() not in SUPPORTED_MODES and mode.lower() != MIXED_MODE]
    if invalid:
        raise ValueError(f"Unsupported topology mode(s): {', '.join(invalid)}")
    return modes


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a LignAmb25-compatible lignin topology from an S/G/H ratio.")
    parser.add_argument("--g", type=float, required=True, help="G ratio")
    parser.add_argument("--s", type=float, required=True, help="S ratio")
    parser.add_argument("--h", type=float, required=True, help="H ratio")
    parser.add_argument("--length", type=int, required=True, help="Number of monomers in the chain")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for monomer ordering")
    parser.add_argument(
        "--mode",
        default="bo4",
        help="Single topology mode or alias: bo4, ao4, b1, bb, b5, c5c5, mixed, all, or all-supported",
    )
    parser.add_argument("--ensemble-size", type=int, default=1, help="Number of topologies to generate")
    parser.add_argument(
        "--modes",
        default="bo4,ao4,b1,bb",
        help="Comma-separated topology modes to cycle through when building an ensemble, or use all/all-supported",
    )
    parser.add_argument("--name", default="lignin_chain", help="Output basename")
    parser.add_argument("--out-dir", default=".", help="Output directory")
    parser.add_argument("--no-run-tleap", action="store_true", help="Only write the tleap file")
    parser.add_argument("--tleap-cmd", default="tleap", help="tleap executable")
    args = parser.parse_args()

    request = TopologyRequest(
        g_ratio=args.g,
        s_ratio=args.s,
        h_ratio=args.h,
        chain_length=args.length,
        seed=args.seed,
        mode=args.mode,
        name=args.name,
        output_dir=Path(args.out_dir),
    )
    mode_alias_requested = args.mode.strip().lower() in {"all", "all-supported"}
    if args.ensemble_size == 1:
        if mode_alias_requested:
            modes = resolve_modes(args.mode, request.chain_length)
            tleap_paths = write_topology_ensemble(
                request,
                modes=modes,
                ensemble_size=len(modes),
                run_tleap=not args.no_run_tleap,
                tleap_cmd=args.tleap_cmd,
            )
        else:
            tleap_paths = [write_topology_bundle(request, run_tleap=not args.no_run_tleap, tleap_cmd=args.tleap_cmd)]
    else:
        if args.mode.strip().lower() == MIXED_MODE:
            tleap_paths = write_mixed_topology_ensemble(
                request,
                ensemble_size=args.ensemble_size,
                run_tleap=not args.no_run_tleap,
                tleap_cmd=args.tleap_cmd,
            )
        else:
            modes = resolve_modes(args.modes, request.chain_length)
            tleap_paths = write_topology_ensemble(
                request,
                modes=modes,
                ensemble_size=args.ensemble_size,
                run_tleap=not args.no_run_tleap,
                tleap_cmd=args.tleap_cmd,
            )

    for tleap_path in tleap_paths:
        print(tleap_path)


if __name__ == "__main__":
    main()
