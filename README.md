# lignamb_topology

This package generates LignAmb25-compatible lignin topology inputs from an S/G/H ratio.

Current scope:

- Takes a target G/S/H ratio plus a chain length and seed.
- Generates a monomer sequence from that ratio.
- Emits a tleap input using the LignAmb25 residue templates and linker conventions.
- Can generate an ensemble by cycling through supported topology modes.
- Supports a `mixed` mode that samples a linkage distribution from the composition and writes a linkage summary JSON alongside the tleap input.
- Supports mixed-linkage ensembles with one summary JSON per candidate plus an aggregate `*.ensemble.json` manifest.
- Runs tleap to produce `.pdb`, `.parm7`, and `.rst7` outputs.
- Supports `all` / `all-supported` mode aliases for composition-preserving mode coverage and 2-mer full-mode coverage, respectively.

Current topology model:

- Supported modes: `bo4`, `ao4`, `b1`, `bb`, `b5`, `c5c5`.
- Uses the matching LignAmb residue templates for each mode, including the non-linear residue families already present in the library.
- Designed to be a deprecation path away from lignin-kmc for topology generation, while still reusing the old algorithmic idea of ratio-driven chain construction.
- The mixed-linkage renderer currently covers BO4, AO4, and B1 paths; BB, 4-O-5 / C5O4, and 5-5 / B5 remain prediction-side classes until the mixed-chain templates are validated locally.

Recommended workflow guide:

- See [COMPOSITION_WORKFLOW.md](COMPOSITION_WORKFLOW.md) for the detailed composition-to-topology procedure, including ratio normalization, chain-length selection, single-structure generation, ensemble generation, and validation guidance.

Example:

```bash
python -m lignamb_topology --g 2.02 --s 0.19 --h 0.38 --length 6 --seed 7 --name lignin_6mer --out-dir out
```

This will write the tleap script and, unless `--no-run-tleap` is given, the parameterizable LignAmb outputs in the requested directory.

For a mixed-linkage prediction, use:

```bash
python -m lignamb_topology --g 2.02 --s 0.19 --h 0.38 --length 6 --seed 7 --mode mixed --name lignin_mixed --out-dir out
```

That writes `lignin_mixed.linkages.json` with the inferred linkage distribution and sampled edge plan.

For several mixed-linkage candidates, add `--ensemble-size`:

```bash
python -m lignamb_topology --g 2.02 --s 0.19 --h 0.38 --length 6 --seed 7 --mode mixed --ensemble-size 4 --name lignin_mixed --out-dir out
```

That writes one `*.linkages.json` per candidate and a top-level `lignin_mixed.ensemble.json` manifest.

For an ensemble, provide a count and a mode cycle:

```bash
python -m lignamb_topology --g 2.02 --s 0.19 --h 0.38 --length 6 --seed 7 --name lignin_6mer --out-dir out --ensemble-size 4 --modes bo4,ao4,b1,bb
```