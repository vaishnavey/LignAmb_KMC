# Composition-Driven Topology Workflow

This guide describes the recommended path for turning a user-provided lignin composition into a LignAmb-compatible topology bundle.

The short version is:

1. Convert the requested composition into normalized G/S/H ratios.
2. Choose a chain length that matches the granularity you want to model.
3. Generate a topology bundle with `lignamb_topology` rather than the older `lignin-kmc` output path.
4. Prefer the `.parm7` and `.rst7` files as the authoritative bonded structure.
5. Use the `.pdb` as a visualization/export convenience, not as the source of truth.

## Why This Is The Recommended Path

The topology-first workflow is the most stable way to build LignAmb-ready structures because it separates three concerns that were previously mixed together:

- composition handling
- topology assembly
- file export for simulation or viewing

That separation matters because the user may ask for a composition in multiple forms, such as:

- a ratio like `G/S/H = 2.02/0.19/0.38`
- absolute composition values in mmol/g
- a desired chain length with a ratio target
- a specific bonding mode such as `bo4`, `ao4`, `b1`, `bb`, `b5`, or `c5c5`

The recommended workflow always reduces those inputs to a clean, reproducible topology request, then lets the package handle the actual chain construction.

## When To Use This Workflow

Use this guide whenever the user gives you a composition and wants one of the following:

- a single representative lignin topology
- a small ensemble of plausible topologies
- LignAmb-compatible `.parm7`, `.rst7`, and `.pdb` outputs
- a reproducible structure that can be regenerated later from the same composition and seed

This is the preferred route for future work. It replaces the old instinct of building a raw lignin-kmc artifact first and then trying to clean it up afterward.

## Input Interpretation

The first job is to interpret the user’s composition correctly.

### If The User Gives Ratios Directly

If the user already gives something like `G = 2.02, S = 0.19, H = 0.38`, treat those as relative weights unless the user explicitly says they are absolute amounts that should be preserved as-is.

The topology builder normalizes them internally, so the exact numerical scale does not matter as long as the relative composition is correct.

### If The User Gives mmol/g Or Another Absolute Basis

If the user gives values in mmol/g, do not try to preserve the units in the topology request itself.

Instead:

1. Treat the values as a composition vector.
2. Normalize them to proportions.
3. Select a chain length that gives you a usable integer monomer count.

Example:

- `G: 2.02`
- `S: 0.19`
- `H: 0.38`

Normalized proportions are based on the sum `2.02 + 0.19 + 0.38 = 2.59`.

That means the approximate fractions are:

- `G ≈ 0.780`
- `S ≈ 0.073`
- `H ≈ 0.147`

The builder then converts those proportions into monomer counts for the requested chain length.

### If The User Gives Only A Target S/G/H Distribution

If the user only gives a qualitative target like “high G, low S, some H”, ask for one of the following:

- an explicit numeric ratio
- a chain length
- both

The package needs a numeric request to produce a reproducible topology.

## Choosing Chain Length

Chain length is the main control knob after composition.

### Recommended Default

For most use cases, choose a short-to-moderate chain length first, then increase it once the composition and bonding mode behave as expected.

Recommended starting points:

- `6` to `10` monomers for quick validation
- `12` to `20` monomers for a more realistic local composition sample
- larger lengths only when you need a broader topology distribution and are comfortable with longer generation and validation time

### How To Think About Length

Longer chains improve composition resolution because integer rounding has less effect.

For example, if the normalized fraction of `H` is small, a very short chain may contain no `H` units at all simply because the rounding drops them.

As a result:

- short chains are good for debugging and visual validation
- longer chains are better when the user cares about preserving the requested composition more faithfully

### Practical Rule

If the user’s request is meant for later simulation, choose the smallest chain length that still captures the intended composition without rounding away one of the monomer types.

## Choosing The Output Mode

The topology mode controls the bond family used to connect the chain.

Supported modes in the current package are:

- `bo4`
- `ao4`
- `b1`
- `bb`
- `b5`
- `c5c5`

### Recommended Default Mode Strategy

If the user does not specify a mode, start with `bo4`.

That is usually the safest default because:

- it is common in lignin topologies
- it exercises the bridge-residue path that is most useful to validate early
- it is a good sanity check for the LignAmb residue library and connectivity handling

### When To Use Ensembles

Use an ensemble when the user wants a family of plausible topologies rather than a single canonical one.

This is the right choice when:

- the user wants several candidates for downstream screening
- the exact bond arrangement is not meant to be unique
- you want to sample several modes from the same composition

The current recommended ensemble cycle is:

- `bo4`
- `ao4`
- `b1`
- `bb`

If the request specifically calls for `b5` or `c5c5`, generate those as targeted single-mode outputs or extend the mode cycle intentionally.

## Recommended User-Facing Workflow

This is the sequence I recommend using when a user provides a composition.

### Step 1: Normalize The Composition

Convert the user’s values to normalized ratios.

If the user gives:

- `G = 2.02`
- `S = 0.19`
- `H = 0.38`

then normalize those values and keep the proportions, not the original units, in the topology request.

### Step 2: Select A Chain Length

Pick a length that matches the purpose:

- debugging or visualization: `6` to `8`
- representative example: `8` to `12`
- broader downstream use: `12+`

If the user already asked for a specific number of monomers, keep that value unless it is clearly too small to preserve the composition.

### Step 3: Decide Between One Structure Or An Ensemble

Use a single structure when the user wants:

- one representative topology
- one reproducible test case
- a simple start for simulation or inspection

Use an ensemble when the user wants:

- multiple plausible topologies for the same composition
- a small panel of modes for screening
- better coverage of the chemistry before choosing a final candidate

### Step 4: Generate The Topology Bundle

Use `lignamb_topology` to write:

- the tleap input
- the `.parm7` topology
- the `.rst7` coordinates
- the `.pdb` visualization file

The package is designed so the `.parm7` and `.rst7` files are the authoritative bonded output.

### Step 5: Validate The Result

Always validate the generated bundle before treating it as final.

Minimum checks:

- tleap completes without errors
- the residue sequence matches the intended composition
- the final residue count matches the chosen chain length
- the `.parm7` bond graph includes the intended inter-residue links
- the `.pdb` opens cleanly in the intended viewer

If the viewer looks odd but the parm topology is correct, trust the parm topology first.

## Recommended CLI Patterns

### Single Topology

Use this when you want one representative chain.

```bash
python -m lignamb_topology \
  --g 2.02 \
  --s 0.19 \
  --h 0.38 \
  --length 6 \
  --seed 7 \
  --mode bo4 \
  --name lignin_6mer \
  --out-dir out
```

Recommended interpretation:

- `--g`, `--s`, `--h` define the composition
- `--length` defines the size of the chain
- `--seed` makes the monomer ordering reproducible
- `--mode` chooses the bond family
- `--out-dir` points to a clean output folder

### Ensemble Of Several Modes

Use this when you want a small panel of related candidates.

If you want the best composition-preserving coverage without listing modes manually, use `--mode all` or `--modes all`. That expands to the chain-compatible linkage families only:

- `bo4`
- `ao4`
- `b1`
- `bb`

If you specifically need the dimer-only modes too, use `--mode all-supported` with `--length 2`.

```bash
python -m lignamb_topology \
  --g 2.02 \
  --s 0.19 \
  --h 0.38 \
  --length 6 \
  --seed 7 \
  --name lignin_6mer \
  --out-dir out \
  --ensemble-size 4 \
  --modes bo4,ao4,b1,bb
```

Recommended interpretation:

- each member is generated from the same composition
- the seed is offset across members so the ordering changes
- the mode cycle gives you a compact cross-section of topology families

## Mixed-Linkage Design For This Composition

If the user wants a design that mixes linkage families for the same composition, the recommended implementation today is an ensemble that cycles through the linear, chain-compatible linkages:

- `bo4`
- `ao4`
- `b1`
- `bb`

That gives a practical mixed-linkage design without forcing the builder to invent unsupported mixed chemistry inside a single linear path.

Example command:

```bash
python -m lignamb_topology \
  --g 2.02 \
  --s 0.19 \
  --h 0.38 \
  --length 6 \
  --seed 7 \
  --name lignin_6mer_mixed \
  --out-dir out \
  --ensemble-size 4 \
  --modes bo4,ao4,b1,bb
```

What this means in practice:

- you get one topology per linkage family
- all members preserve the same target composition
- the resulting set is the best current stand-in for a single mixed-linkage design

Important limitation:

- `b5` and `c5c5` remain dimer-only in the current builder
- if you want those included, generate them as separate targeted examples rather than forcing them into the same linear chain workflow

## Script Order

When a user gives you a composition, run the scripts in this order.

### 1. Use `lignamb_topology` To Build The Primary Bundle

This is the default entry point for composition-driven work.

Run it first when you want a single structure:

```bash
python -m lignamb_topology \
  --g 2.02 \
  --s 0.19 \
  --h 0.38 \
  --length 6 \
  --seed 7 \
  --mode bo4 \
  --name lignin_6mer \
  --out-dir out
```

What this produces:

- `lignin_6mer.tleap`
- `lignin_6mer.parm7`
- `lignin_6mer.rst7`
- `lignin_6mer.pdb`

Use this output as the primary candidate unless the user explicitly asked for more than one topology.

### 2. Use `lignamb_topology` In Ensemble Mode If The User Wants Multiple Candidates

Run the ensemble form second if the request is for a small set of plausible structures rather than one canonical chain.

```bash
python -m lignamb_topology \
  --g 2.02 \
  --s 0.19 \
  --h 0.38 \
  --length 6 \
  --seed 7 \
  --name lignin_6mer \
  --out-dir out \
  --ensemble-size 4 \
  --modes bo4,ao4,b1,bb
```

What this produces:

- one subdirectory per ensemble member
- a full tleap script for each member
- matching `.parm7`, `.rst7`, and `.pdb` files for each topology mode

Recommended use:

- compare the members visually
- inspect which bond family best fits the requested composition
- keep the `.parm7` file that matches the final downstream intent

### 3. Inspect The Topology With The Parm Files Before Trusting The PDB

After generation, inspect the topology files if the viewer looks suspicious.

Example check:

```bash
python - <<'PY'
from pathlib import Path
import parmed as pmd

base = Path('out/lignin_6mer')
struct = pmd.load_file(str(base / 'lignin_6mer.parm7'), str(base / 'lignin_6mer.rst7'))
print('residues:', len(struct.residues))
print('last residue:', struct.residues[-1].name, struct.residues[-1].number)
PY
```

Use this kind of check when the viewer is not obviously matching the intended connectivity.

### 4. Use `build_lignamb_chain.py` Only For Legacy Or Graph-Driven Input

Run this only when you are starting from an existing graph export or a KMC-derived structure description.

It is not the preferred starting point for a fresh composition request.

Example legacy flow:

```bash
python build_lignamb_chain.py \
  --spec legacy_graph.json \
  --out-dir out_legacy
```

Use this path when:

- the input is already a monomer/edge graph
- you are reproducing an older chain construction result
- you need to translate a graph export into tleap operations

### 5. Use `make_pdbs_from_smiles.py` Only For The Older SMILES-Based Route

This script is for the older, chemistry-to-coordinate path, not the recommended composition-first topology path.

Use it only if you are intentionally working with a SMILES-derived workflow.

Example:

```bash
cd /data/srinivab/lignin_kmc_work/lignin-kmc
conda run -n ambertools python make_pdbs_from_smiles.py
```

Why it is last in the order:

- it is a legacy utility rather than the new user-facing path
- it produces coordinates from a different starting point than the topology builder
- it should not be used to replace a topology-first workflow for new composition requests

### Recommended Overall Order

For a new composition request, the practical order is:

1. run `python -m lignamb_topology ...` for a single topology
2. run `python -m lignamb_topology ... --ensemble-size ...` if the user wants multiple candidates
3. inspect the `.parm7` and `.rst7` files
4. only use `build_lignamb_chain.py` if the input is already a graph export
5. only use `make_pdbs_from_smiles.py` if the user explicitly wants the older SMILES route

## How To Handle Viewer Issues

For later use, remember the following rule:

The topology graph in `.parm7` is the source of truth.

The `.pdb` is only for visualization and handoff.

That matters because PDB viewers can misread bridge residues, especially around terminal linkages. If the viewer looks broken but the parm topology is correct, the right response is usually to regenerate or rewrite the PDB representation, not to change the underlying chemistry.

Practical checks:

- inspect the terminal residue count in `.parm7`
- confirm the intended residue is the actual chain end
- compare the bond list near the terminal bridge
- do not rely on visual continuity alone in the PDB viewer

## What To Do When The User Gives A New Composition

Use this decision sequence:

1. Convert the composition to normalized G/S/H proportions.
2. Select a chain length that preserves the composition after integer rounding.
3. Decide whether the user wants one topology or several.
4. Choose `bo4` as the default single-mode starting point unless the user requests otherwise.
5. Generate the topology bundle with `lignamb_topology`.
6. Validate the `.parm7` bond graph and the final residue count.
7. Only then hand off the `.pdb` for viewing or the `.parm7`/`.rst7` for simulation.

## Recommended Output Contract

When you build a final output set for a user, keep the following contract:

- `*.tleap` for reproducibility
- `*.parm7` for the bonded topology
- `*.rst7` for coordinates
- `*.pdb` for convenience and visualization

If you need to choose one file to trust for chemistry, choose `*.parm7`.

## Why We Are Deprecating The Old Path

The old `lignin-kmc` path is still useful as algorithmic reference material, but it should not be the primary user-facing output path anymore.

The reason is simple:

- the topology package is easier to validate
- it produces direct LignAmb-compatible outputs
- it is clearer for users who start from composition rather than from a low-level graph
- it avoids the extra cleanup step that was previously needed to make the output usable

In other words, keep the old code as a backend idea when needed, but prefer the new topology package for real work.

## Final Recommendation

For a future user request with a composition, the default response path should be:

1. normalize the composition
2. choose a modest chain length
3. generate a single `bo4` example first unless the user wants an ensemble
4. validate the resulting `.parm7`
5. expand to other modes only if the request needs them

That is the most reliable and maintainable workflow for producing LignAmb-compatible lignin topologies from a composition request.