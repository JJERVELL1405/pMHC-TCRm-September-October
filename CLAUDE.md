# TCRm structure preparation — project notes

Preparing the TCR-mimic antibody set for the pMHC-TCRm thesis (Bonvin lab).
Supervisor sent `prep_structures.py` plus a worked example
(`renumbered_dataset_pMHC_TCR.pse`, the already-prepared TCR-pMHC set).
Task: run the same script over the TCRm antibody structures.

## Environment

conda env `tcrprep` (macOS). Activate before any run:

```
conda activate tcrprep
```

Dependencies: `anarci` (bioconda, noarch), `hmmer`, `numpy`, `scipy`.
ANARCI installed is 2024.05.21. `scipy` is optional in the code but the
greedy fallback is documented as wrong on 2AK4 — keep it installed.

## The ANARCI patch — IMPORTANT

`prep_structures.py` line 198 unpacks 3 values from `run_anarci`:

```python
numbering, details, _hits = run_anarci(
```

ANARCI 2024.05.21 returns **four** (`sequences, numbered, alignment_details,
hit_tables`). Every structure fails with
`ValueError: too many values to unpack (expected 3)`.

Fix lives in `prep_structures_patched.py`, which differs from the original by
that one line only (`_seqs, ` prepended). **Always run the patched copy.**
Original is kept unmodified for diffing. Not yet raised with supervisor.

## Commands

```
conda activate tcrprep
cd ~/Documents/Tcrm_Prep30Sep
python3 prep_structures_patched.py --ids TCR3d_data_2026-10-01.csv \
    --out sets/tcrm_ab --verbose 2>&1 | tee prep.log
```

No `--from-dir` — the script downloads from RCSB into `pdb_cache/`
(already populated, so re-runs are fast).

`manifest.tsv` is opened in write mode each run, so a partial run
**truncates it**. Never process a subset into `sets/tcrm_ab`; use a
separate `--out` for one-off tests.

## Output format

Two chains per complex, band-offset so residue numbers don't collide:

| File | Chain | Contents |
|---|---|---|
| `<id>_tcr.pdb` | A | rec1 +0, rec2 +1000 (IMGT, via ANARCI) |
| `<id>_pmhc.pdb` | B | MHC +0, b2m +1000, peptide +2000 |
| `<id>_complex.pdb` | A+B | both — the caprieval reference |

So CDR3 of rec1 is chain A 105-117, CDR3 of rec2 is chain A 1105-1117,
peptide is chain B 2001+.

## Input data

- `TCR3d_data.csv` — downloaded 2026-09-20, 48 rows, **missing 33CI**
- `TCR3d_data_2026-10-01.csv` — 49 rows, use this one

Source: https://tcr3d.ibbr.umd.edu/tcrm_ab ("Download table data as CSV").
33CG and 33CI have identical docking (147.5) and incident (7.6) angles —
verify whether that is a real near-identity or a duplicated row before
using either in an angle distribution.

## Known issues in the prepared set

**9L1L — exclude.** Class II (I-Ak). Both MHC chains (179 and 172 res)
classify as `mhc` because the rule is `n >= 150`; `pick_unit` keeps one and
drops the other. It passes validation because its peptide is unfused, so
all four of `want` are present. `9l1l_pmhc.pdb` contains half an MHC.
Manifest tell: `mhc 172, b2m 0`.

**3CVH — flag.** Mouse H-2Kb. Its β2m chains (99 res) classify as `None`,
not `b2m`: the `B2M` reference at line 73 is the human sequence and
`kmer_sim(s, B2M) > 0.25` fails for mouse. β2m is dropped. Passes because
`b2m` is not in `want`. Only mouse class I entry in the set, but a latent
bug for any added later.

**MHC/b2m numbering is not comparable across structures.**
`renumber_serial` (line 491) renumbers MHC, b2m and peptide consecutively
from 1, erasing crystallographic gaps. 9PKV has 21 missing residues across
three α3 loops; after prep they are consecutive. The receptor side is
IMGT-numbered and is fine. Any cross-structure comparison of *MHC residue
positions* is unsafe. IMGT defines numbering for G domains — not applied here.

**9d73 / 9d74** have `cdr3_1` = 6 (typical 11-17). Check against SEQRES
whether CDR-H3 is genuinely short or partly unmodelled before using their
contact counts. Same antibody (B1.23.2) as 8TQ6.

## Refusals (run of 48, 2026-10-01)

35 prepared, 13 refused:

- **Class II, peptide fused into β chain** — 6XP6, 8W83, 8W84, 8W85, 8W86
- **Single-chain trimer (peptide + b2m fused)** — 8TNJ
- **scFv** — 8EK5. Line 207 takes `num[0]` only, so VL is dropped
- **VHH / single domain** — 9HKQ. No light chain exists
- **Cross-copy Fv pairing** — 7RE7, 8TQ6. Guard at line 596 (<50 CA contacts).
  `--chain-map` can't be built from the TCR3d CSV (no chain columns).
  Try `--copy 1`
- **Download failed** — 9NFB, 9NFC, 9O55. Retry; if persistent, check
  whether RCSB serves PDB format for them at all

9PKV, 9YTD and 9YTF contain a nanobody alongside a conventional antibody.
All three picked the correct paired Fv — verified.

## Open questions for supervisor

1. ANARCI version — line 198 unpacks 3, mine returns 4. Which version is hers?
2. Should class II be in the TCRm set? If so the chain logic needs changing;
   if not, it should be excluded explicitly rather than relying on failure
   (9L1L slipped through).
3. 8EK5 (scFv) and 9HKQ (VHH) — include or drop?
4. Is cross-structure MHC position comparison intended? If so, G-domain
   IMGT numbering rather than serial.
5. Does she expect all 49 TCR3d entries, or a filtered subset?

## Project phase

Six-phase thesis plan; this is **phase 2** (dataset curation and structure
prep). Phase 3 is geometric interface characterisation (PyMOL, Python),
phase 4 is energetic analysis (HADDOCK3, PyRosetta). Docking is not the
current deliverable.
