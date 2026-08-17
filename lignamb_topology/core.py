from __future__ import annotations

import json
import sys
import shutil
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from random import Random
from typing import Any

from build_lignamb_chain import build_component_ops, chain_spec_from_kmc_result, render_tleap, resolve_tleap_executable


DEFAULT_LEAPRC = "leaprc.LignAmb25_HF"
MONOMER_TYPES = ("G", "S", "H")
SUPPORTED_MODES = ("bo4", "ao4", "b1", "bb", "b5", "c5c5")
CHAIN_COMPATIBLE_MODES = ("bo4", "ao4", "b1", "bb")
DIMER_ONLY_MODES = ("b5", "c5c5")
MIXED_MODE = "mixed"
PREDICTED_LINKAGE_LABELS = ("bo4", "ao4", "b1", "bb", "b5", "c5o4", "c5c5")
RENDERABLE_LINKAGE_LABELS = ("bo4", "ao4", "b1")


@dataclass(frozen=True)
class TopologyRequest:
    """User input for a single LignAmb-compatible lignin topology."""

    g_ratio: float
    s_ratio: float
    h_ratio: float
    chain_length: int
    seed: int = 1
    mode: str = "bo4"
    name: str = "lignin_chain"
    output_dir: Path | str = Path(".")


def _normalize_distribution(scores: dict[str, float]) -> dict[str, float]:
    filtered = {label: value for label, value in scores.items() if value > 0}
    total = sum(filtered.values())
    if total <= 0:
        raise ValueError("Could not derive a positive linkage distribution.")
    return {label: value / total for label, value in filtered.items()}


@lru_cache(maxsize=1)
def _load_kmc_rate_data() -> tuple[dict[str, Any], dict[str, str]] | None:
    try:
        from ligninkmc.kmc_common import BO4, AO4, B1, BB, B5, C5O4, C5C5, DEF_RXN_RATES, G, H, MON_MON, S
    except ImportError:
        workspace_root = Path(__file__).resolve().parents[2]
        candidate = workspace_root / "lignin_kmc_work" / "lignin-kmc"
        if candidate.is_dir() and str(candidate) not in sys.path:
            sys.path.append(str(candidate))
        try:
            from ligninkmc.kmc_common import BO4, AO4, B1, BB, B5, C5O4, C5C5, DEF_RXN_RATES, G, H, MON_MON, S
        except ImportError:
            return None

    return DEF_RXN_RATES, {"G": G, "S": S, "H": H}, {"bo4": BO4, "ao4": AO4, "b1": B1, "bb": BB, "b5": B5, "c5o4": C5O4, "c5c5": C5C5}, MON_MON


def normalize_ratio(request: TopologyRequest) -> dict[str, float]:
    total = request.g_ratio + request.s_ratio + request.h_ratio
    if total <= 0:
        raise ValueError("At least one of G, S, or H must be positive.")
    return {
        "G": request.g_ratio / total,
        "S": request.s_ratio / total,
        "H": request.h_ratio / total,
    }


def allocate_counts(request: TopologyRequest) -> dict[str, int]:
    if request.chain_length <= 0:
        raise ValueError("chain_length must be positive.")

    ratios = normalize_ratio(request)
    exact = {mono: ratios[mono] * request.chain_length for mono in MONOMER_TYPES}
    counts = {mono: int(value) for mono, value in exact.items()}
    remainder = request.chain_length - sum(counts.values())

    ranking = sorted(((exact[mono] - counts[mono], mono) for mono in MONOMER_TYPES), reverse=True)
    for _, mono in ranking[:remainder]:
        counts[mono] += 1
    return counts


def sample_sequence(counts: dict[str, int], seed: int) -> list[str]:
    pool: list[str] = []
    for mono, count in counts.items():
        pool.extend([mono] * count)
    rng = Random(seed)
    rng.shuffle(pool)
    return pool


def _predict_linkage_distribution(request: TopologyRequest) -> dict[str, float]:
    ratios = normalize_ratio(request)
    kmc_data = _load_kmc_rate_data()
    if kmc_data is None:
        scores = {
            "bo4": 0.50 + 0.20 * ratios["G"] + 0.10 * ratios["S"],
            "ao4": 0.03 + 0.01 * ratios["S"],
            "b1": 0.12 + 0.12 * ratios["H"] + 0.04 * ratios["G"],
            "bb": 0.10 + 0.08 * ratios["G"],
            "b5": 0.05 + 0.06 * ratios["G"],
            "c5o4": 0.04 + 0.05 * ratios["H"],
            "c5c5": 0.04 + 0.03 * ratios["G"],
        }
        return _normalize_distribution(scores)

    rates, labels, kmc_labels, mon_mon = kmc_data
    scores: dict[str, float] = {label: 0.0 for label in PREDICTED_LINKAGE_LABELS}
    for bond_label, kmc_label in kmc_labels.items():
        bond_rates = rates[kmc_label]
        for left_symbol, left_fraction in ratios.items():
            for right_symbol, right_fraction in ratios.items():
                left_label = labels[left_symbol]
                right_label = labels[right_symbol]
                pair_rates = bond_rates.get((left_label, right_label)) or bond_rates.get((right_label, left_label))
                if pair_rates is None:
                    pair_rates = next(iter(bond_rates.values()))
                score = pair_rates.get(mon_mon)
                if score is None:
                    score = next(iter(pair_rates.values()))
                scores[bond_label] += float(score) * left_fraction * right_fraction
    return _normalize_distribution(scores)


def _renderable_linkage_distribution(linkage_distribution: dict[str, float]) -> dict[str, float]:
    renderable = {label: linkage_distribution[label] for label in RENDERABLE_LINKAGE_LABELS if label in linkage_distribution}
    return _normalize_distribution(renderable)


def _sample_linkage_sequence(distribution: dict[str, float], count: int, seed: int) -> list[str]:
    if count <= 0:
        return []
    labels = list(distribution)
    weights = [distribution[label] for label in labels]
    rng = Random(seed)
    return rng.choices(labels, weights=weights, k=count)


def _sample_linkage_counts(linkages: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for linkage in linkages:
        counts[linkage] = counts.get(linkage, 0) + 1
    return counts


def _composition_metadata(request: TopologyRequest, counts: dict[str, int], sequence: list[str]) -> dict[str, Any]:
    return {
        "input_composition": {
            "G": request.g_ratio,
            "S": request.s_ratio,
            "H": request.h_ratio,
        },
        "normalized_composition": normalize_ratio(request),
        "chain_length": request.chain_length,
        "counts": counts,
        "sequence": sequence,
        "seed": request.seed,
        "mode": request.mode,
    }


def _terminal_template(mono: str) -> str:
    return {"G": "1G1", "S": "1S1", "H": "1H1"}[mono]


def _open_template(mono: str) -> str:
    return {"G": "2G4", "S": "2S4", "H": "2H4"}[mono]


def _mode_to_bond_type(mode: str) -> str:
    normalized = mode.strip().lower()
    if normalized not in SUPPORTED_MODES:
        raise ValueError(f"Unsupported topology mode: {mode!r}")
    return {
        "bo4": "BO4",
        "ao4": "AO4",
        "b1": "B1",
        "bb": "BB",
        "b5": "B5",
        "c5c5": "C5C5",
    }[normalized]


def _sequence_to_path_spec(
    sequence: list[str],
    mode: str,
    spec_name: str,
    ratio_label: str,
) -> dict[str, Any]:
    if not sequence:
        raise ValueError("Could not derive a monomer sequence from the requested ratio.")

    bond_type = _mode_to_bond_type(mode)
    monomers = {index: {"type": mono} for index, mono in enumerate(sequence)}
    order = list(range(len(sequence)))
    edges = [(left, right, bond_type) for left, right in zip(order, order[1:])]
    ops = build_component_ops(order, monomers, edges, spec_name)

    return {
        "leaprc": DEFAULT_LEAPRC,
        "name": spec_name,
        "mol_name": "chain",
        "chains": [
            {
                "name": spec_name,
                "title": f"{mode.upper()} lignin chain for G/S/H = {ratio_label}",
                "ops": ops,
            }
        ],
    }


def _mixed_path_spec(request: TopologyRequest, sequence: list[str]) -> dict[str, Any]:
    linkage_distribution = _predict_linkage_distribution(request)
    renderable_distribution = _renderable_linkage_distribution(linkage_distribution)
    sampled_bonds = _sample_linkage_sequence(renderable_distribution, max(0, request.chain_length - 1), request.seed + 1)

    monomers = [{"type": mono, "identity": index} for index, mono in enumerate(sequence)]
    edges = [
        {"left": left, "right": right, "bond_type": bond_type.upper()}
        for (left, right), bond_type in zip(zip(range(len(sequence)), range(1, len(sequence))), sampled_bonds)
    ]
    spec = chain_spec_from_kmc_result({"monomers": monomers, "edges": edges}, spec_name=request.name)
    spec.setdefault("metadata", {})
    spec["metadata"].update(
        {
            "predicted_linkage_distribution": linkage_distribution,
            "renderable_linkage_distribution": renderable_distribution,
            "sampled_linkages": sampled_bonds,
            "sampled_linkage_counts": _sample_linkage_counts(sampled_bonds),
        }
    )
    return spec


def generate_linear_bo4_spec(request: TopologyRequest) -> dict[str, Any]:
    return generate_topology_spec(request)


def generate_topology_spec(request: TopologyRequest) -> dict[str, Any]:
    counts = allocate_counts(request)
    sequence = sample_sequence(counts, request.seed)
    mode = request.mode.strip().lower()
    if mode == MIXED_MODE:
        spec = _mixed_path_spec(request, sequence)
    else:
        spec = _sequence_to_path_spec(sequence, request.mode, request.name, f"{request.g_ratio}:{request.s_ratio}:{request.h_ratio}")

    metadata = dict(spec.get("metadata", {}))
    metadata.update(_composition_metadata(request, counts, sequence))
    spec.update(
        {
            "output": {
                "parm7": f"{request.name}.parm7",
                "rst7": f"{request.name}.rst7",
                "pdb": f"{request.name}.pdb",
            },
            "metadata": metadata,
        }
    )
    return spec


def _rewrite_pdb_from_topology(output_dir: Path, name: str) -> None:
    import parmed as pmd

    parm7_path = output_dir / f"{name}.parm7"
    rst7_path = output_dir / f"{name}.rst7"
    pdb_path = output_dir / f"{name}.pdb"

    if not (parm7_path.exists() and rst7_path.exists() and pdb_path.exists()):
        return

    struct = pmd.load_file(str(parm7_path), str(rst7_path))
    scratch_path = output_dir / f".{name}.parmed.pdb"
    struct.save(str(scratch_path), overwrite=True)

    lines = [line for line in scratch_path.read_text().splitlines() if not line.startswith("TER") and not line.startswith("END")]

    conect_map: dict[int, set[int]] = {}
    for bond in struct.bonds:
        left = bond.atom1.idx + 1
        right = bond.atom2.idx + 1
        conect_map.setdefault(left, set()).add(right)
        conect_map.setdefault(right, set()).add(left)

    for atom_index in sorted(conect_map):
        partners = sorted(conect_map[atom_index])
        for start in range(0, len(partners), 4):
            chunk = partners[start:start + 4]
            lines.append(f"CONECT{atom_index:5d}" + "".join(f"{partner:5d}" for partner in chunk))

    lines.append("END")
    pdb_path.write_text("\n".join(lines) + "\n")
    scratch_path.unlink(missing_ok=True)


def _write_topology_bundle_from_spec(
    request: TopologyRequest,
    spec: dict[str, Any],
    *,
    run_tleap: bool = True,
    tleap_cmd: str = "tleap",
) -> Path:
    output_dir = Path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tleap_path = output_dir / f"{request.name}.tleap"
    tleap_path.write_text(render_tleap(spec))

    if request.mode.strip().lower() == MIXED_MODE:
        summary_path = output_dir / f"{request.name}.linkages.json"
        summary_path.write_text(json.dumps(spec.get("metadata", {}), indent=2, sort_keys=True) + "\n")

    if run_tleap:
        import subprocess

        tleap_executable = resolve_tleap_executable(tleap_cmd)
        completed = subprocess.run(
            [tleap_executable, "-f", tleap_path.name],
            cwd=str(output_dir),
            check=False,
            capture_output=True,
            text=True,
        )
        tleap_output = (completed.stdout or "") + "\n" + (completed.stderr or "")
        if "Exiting LEaP: Errors = 0" not in tleap_output:
            raise RuntimeError(tleap_output.strip() or f"tleap failed for {request.name}")
        _rewrite_pdb_from_topology(output_dir, request.name)

    return tleap_path


def write_topology_bundle(request: TopologyRequest, *, run_tleap: bool = True, tleap_cmd: str = "tleap") -> Path:
    spec = generate_topology_spec(request)
    return _write_topology_bundle_from_spec(request, spec, run_tleap=run_tleap, tleap_cmd=tleap_cmd)


def write_topology_ensemble(
    request: TopologyRequest,
    *,
    modes: list[str],
    ensemble_size: int,
    run_tleap: bool = True,
    tleap_cmd: str = "tleap",
) -> list[Path]:
    if ensemble_size <= 0:
        raise ValueError("ensemble_size must be positive.")
    if not modes:
        raise ValueError("At least one topology mode is required for an ensemble.")

    output_dir = Path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tleap_paths: list[Path] = []
    for index in range(ensemble_size):
        mode = modes[index % len(modes)]
        member_name = f"{request.name}_{index + 1:02d}_{mode.lower()}"
        member_dir = output_dir / member_name
        member_request = TopologyRequest(
            g_ratio=request.g_ratio,
            s_ratio=request.s_ratio,
            h_ratio=request.h_ratio,
            chain_length=request.chain_length,
            seed=request.seed + index,
            mode=mode,
            name=member_name,
            output_dir=member_dir,
        )
        tleap_paths.append(write_topology_bundle(member_request, run_tleap=run_tleap, tleap_cmd=tleap_cmd))

    return tleap_paths


def write_mixed_topology_ensemble(
    request: TopologyRequest,
    *,
    ensemble_size: int,
    run_tleap: bool = True,
    tleap_cmd: str = "tleap",
) -> list[Path]:
    if ensemble_size <= 0:
        raise ValueError("ensemble_size must be positive.")

    output_dir = Path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    count_seed_sequence = sample_sequence(allocate_counts(request), request.seed)
    manifest: dict[str, Any] = {
        "name": request.name,
        "ensemble_size": ensemble_size,
        "composition": _composition_metadata(request, allocate_counts(request), count_seed_sequence),
        "members": [],
        "rejected_members": [],
    }

    tleap_paths: list[Path] = []
    accepted_index = 0
    attempt_index = 0
    max_attempts = max(ensemble_size * 10, ensemble_size)
    while accepted_index < ensemble_size and attempt_index < max_attempts:
        index = attempt_index
        member_name = f"{request.name}_{index + 1:02d}_mixed"
        member_dir = output_dir / member_name
        member_request = TopologyRequest(
            g_ratio=request.g_ratio,
            s_ratio=request.s_ratio,
            h_ratio=request.h_ratio,
            chain_length=request.chain_length,
            seed=request.seed + attempt_index,
            mode=MIXED_MODE,
            name=member_name,
            output_dir=member_dir,
        )
        spec = generate_topology_spec(member_request)
        try:
            tleap_path = _write_topology_bundle_from_spec(member_request, spec, run_tleap=run_tleap, tleap_cmd=tleap_cmd)
        except RuntimeError as exc:
            manifest["rejected_members"].append(
                {
                    "name": member_name,
                    "seed": member_request.seed,
                    "reason": str(exc).splitlines()[-1] if str(exc).strip() else "tleap failed",
                }
            )
            shutil.rmtree(member_dir, ignore_errors=True)
            attempt_index += 1
            continue

        tleap_paths.append(tleap_path)
        manifest["members"].append(
            {
                "name": member_name,
                "seed": member_request.seed,
                "tleap": str(tleap_path),
                "summary": str(member_dir / f"{member_name}.linkages.json"),
                "metadata": spec.get("metadata", {}),
            }
        )
        accepted_index += 1
        attempt_index += 1

    if accepted_index < ensemble_size:
        raise RuntimeError(
            f"Only accepted {accepted_index} of {ensemble_size} mixed candidates after {attempt_index} attempts."
        )

    (output_dir / f"{request.name}.ensemble.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return tleap_paths
