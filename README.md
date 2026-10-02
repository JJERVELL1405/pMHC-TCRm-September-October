 # haddock_py

A Python layer over HADDOCK3, for the pMHC–TCRm work.

## What this is and isn't

HADDOCK3 is a command-line program. `haddock.py` doesn't replace that — it shells out to the same commands you'd type. What it buys you:

- **Repeatability.** The prep pipeline is written down once, not retyped each time with a different typo.
- **Config generation.** Restraints are per-module in HADDOCK3; omit them from one stage and that stage silently runs unrestrained. Generating configs in code removes that class of error, and makes parameter sweeps trivial.
- **Analysis in pandas.** This is where it genuinely beats the terminal. `capri_ss.tsv` has thirty-odd columns and wraps unreadably in a shell. Comparing five runs from the command line is miserable.

The docking itself still runs as a subprocess. Nothing changes about the compute.

## Install

```bash
# put haddock.py somewhere permanent
mkdir -p ~/Documents/haddock_py
# copy haddock.py, the notebook, and this file there

conda activate haddock3
pip install pandas matplotlib jupyter
```

Then from anywhere:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path.home() / "Documents" / "haddock_py"))
import haddock

haddock.check_install()
```

**Jupyter must be launched from the `haddock3` environment**, or none of the binaries are on PATH:

```bash
conda activate haddock3
cd ~/Documents/haddock-tutorial/HADDOCK3-antibody-antigen
jupyter notebook
```

## Reference

### Setup
| Function | Does |
|---|---|
| `check_install()` | confirm the binaries are reachable |
| `module_params(module)` | every parameter a module takes, with defaults — do this rather than trusting a copied config, since parameters shift between releases |

### Structure preparation
| Function | Does |
|---|---|
| `prepare_chain(pdb_id, chain, out, resrange=)` | fetch from the PDB and clean |
| `prepare_local_chain(file, chain, out)` | same from a local file — you'll use this more |
| `merge_chains([files], out, chain_id=)` | merge and renumber into one chain |
| `assign_chain(file, out, chain_id)` | clean a single-chain molecule, give it an ID |
| `make_ensemble([models], out)` | combine conformers for ensemble docking |
| `chain_summary(pdb)` | residue ranges per chain — sanity-check after merging |
| `numbering_map(orig, renumbered)` | map old residue numbers to new |

### Restraints
| Function | Does |
|---|---|
| `passive_from_active(pdb, active)` | surface neighbours, filtered at 15% accessibility |
| `write_actpass(path, active, passive)` | write the two-line .act-pass file |
| `make_ambig(ap1, ap2, out, segid_one=, segid_two=)` | generate the CNS TBL |
| `validate_tbl(tbl)` | syntax check |
| `check_restraint_residues(tbl, {segid: pdb})` | verify the residues actually exist |
| `restrain_bodies(pdb, out)` | CA–CA restraints for merged chains |

### Configuration and running
| Function | Does |
|---|---|
| `write_config(...)` | generate a .cfg, injecting restraints into every module that needs them |
| `setup_only(cfg)` | parse without executing |
| `run(cfg, log=, background=True)` | launch |

`CLASSIC_WORKFLOW` and `ENSEMBLE_WORKFLOW` are provided as starting points.

### Analysis
| Function | Does |
|---|---|
| `load_clusters(run_dir)` | capri_clt.tsv → DataFrame — **this is what you report** |
| `load_models(run_dir)` | capri_ss.tsv → DataFrame |
| `stage_progression(run_dir)` | how quality and ranking evolve across modules |
| `clusters_distinguishable(clusters)` | are the top clusters actually separated? |
| `compare_runs({label: dir})` | top clusters side by side |
| `add_quality(df)` | CAPRI quality classification |
| `traceback(run_dir)` | which input conformer produced which model |
| `contacts(run_dir, cluster=)` | intermolecular contacts |
| `score_model(pdb)` | score a single complex |

## Three things worth remembering

**Report clusters, not single models.** The best model is not reliably ranked first — the scoring function is imperfect. `stage_progression()` shows you this directly via `rank_of_best`.

**Check whether your top clusters are distinguishable.** If the top two overlap within their standard deviations, you haven't discriminated between them. Saying so is the correct conclusion, not a failed run. `clusters_distinguishable()` does the arithmetic.

**`validate_tbl` only checks syntax.** A restraint naming residue 500 of a 250-residue protein passes validation and then does nothing. `check_restraint_residues()` catches it. This matters most right after merging and renumbering, which is exactly when numbering mistakes happen.

## For the pMHC–TCRm work

Concretely, what carries over:

- `merge_chains` on both sides — pMHC is heavy + β2m + peptide, TCR is α + β
- `restrain_bodies` on both, non-optional
- `numbering_map` saved alongside your restraint files, because post-merge numbering is where silent errors live
- The active-residue set **is** your specificity hypothesis. Peptide-only vs peptide-plus-groove vs groove-only are different scientific claims. `compare_runs` across those is a real experiment.
- `make_ensemble` + the ensemble workflow when you have several plausible CDR3 conformations and no principled way to choose. The tutorial's Bonus 3 finding — AlphaFold underperforming because of its H3 loop — is your CDR3 problem in a different costume.

## Version note

Written against HADDOCK3 2026.8.0. Module parameters change between releases. If a config is rejected, check with `haddock3-cfg -m <module> -d` rather than trusting any written config, including the templates here.
