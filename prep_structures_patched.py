#!/usr/bin/env python3
"""Fetch, clean, identify, renumber and merge TCR-pMHC crystal structures.

This is the front of the benchmark pipeline: it turns a list of PDB IDs into
the two-chain, band-numbered files every downstream script in this project
expects.

    <out>/<pdbid>_tcr.pdb       chain A: alpha +0, beta +1000   (IMGT)
    <out>/<pdbid>_pmhc.pdb      chain B: MHC +0, b2m +1000, peptide +2000
    <out>/<pdbid>_complex.pdb   both, for sets/refs_capri
    <out>/manifest.tsv          what was found, per target

WHY IMGT NUMBERING IS DONE BY ANARCI AND NOT BY HAND
----------------------------------------------------
CDR3 is IMGT 105-117.  In sequential 1..N numbering those same integers name
C-terminal framework residues, and measuring them produces a confident wrong
answer -- this project has already paid for that once.  IMGT positions cannot
be derived from a structure by counting: they come from an alignment to
germline, with gaps at defined positions and INSERTION CODES (111A, 112A ...)
for long CDR3s.  ANARCI does that alignment properly.  Anything cheaper is a
guess, so this script refuses to run without it rather than approximating.

INSERTION CODES ARE DATA, NOT DECORATION
----------------------------------------
A long CDR3 carries 111A/112A.  Every residue key here is (chain, number,
icode).  Keying on the number alone silently merges 111 and 111A, and their
atoms overwrite one another -- the "residue" you then measure is a blend of
two.  The validation step checks for exactly that.

MULTIPLE COPIES IN THE ASYMMETRIC UNIT
--------------------------------------
8RLT holds two complete TCR-pMHC complexes (chains A-E and F-J).  Keeping
both would double every chain.  Chains are therefore clustered by contact and
one complete unit is taken; --copy chooses which when there is more than one.

Usage
-----
    python3 prep_structures.py --ids targets.txt --out sets/prepared
    python3 prep_structures.py --ids 8rlt,1ao7,7dzm --out /tmp/try --verbose
"""
import argparse
import collections
import gzip
import os
import sys
import urllib.request

CHAIN_MAP = {}
RCSB = "https://files.rcsb.org/download/%s.pdb"

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
       "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
       "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
       "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
# residues that are really standard amino acids wearing a different label
ALIAS = {"MSE": "MET", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
         "CSO": "CYS", "SEC": "CYS", "PTR": "TYR", "SEP": "SER",
         "TPO": "THR", "MLY": "LYS", "KCX": "LYS"}

BANDS = {"rec1": 0, "rec2": 1000, "mhc": 0, "b2m": 1000, "peptide": 2000}

# ANARCI chain type -> which receptor slot it fills.  Covers alpha/beta and
# gamma/delta TCRs and, for the TCR-mimic antibody sets, heavy/light -- a
# TCRm-Ab complex is a pMHC with an antibody Fv where the TCR would be, so
# everything downstream is identical once the chains are in the right slot.
CHAIN_SLOT = {"A": "rec1", "G": "rec1", "H": "rec1",
              "B": "rec2", "D": "rec2", "K": "rec2", "L": "rec2"}
RECEPTOR = {("A", "B"): "TCRab", ("G", "D"): "TCRgd",
            ("H", "K"): "Ab-HK", ("H", "L"): "Ab-HL"}

# beta2-microglobulin is highly conserved; used only to tell it apart from
# an MHC class II chain of similar length
B2M = ("IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLL"
       "YYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM")


# ------------------------------------------------------------------ fetch

def fetch(pdbid, cache):
    """Download from RCSB unless already cached.  Returns the local path."""
    os.makedirs(cache, exist_ok=True)
    dst = os.path.join(cache, "%s.pdb" % pdbid.lower())
    if os.path.exists(dst) and os.path.getsize(dst) > 2000:
        return dst
    url = RCSB % pdbid.upper()
    try:
        with urllib.request.urlopen(url, timeout=60) as r, open(dst, "wb") as fh:
            fh.write(r.read())
    except Exception as e:
        return None
    return dst if os.path.getsize(dst) > 2000 else None


# ------------------------------------------------------------------ clean

class Residue(object):
    __slots__ = ("chain", "num", "icode", "name", "atoms")

    def __init__(self, chain, num, icode, name):
        self.chain, self.num, self.icode, self.name = chain, num, icode, name
        self.atoms = []          # raw ATOM lines, already normalised

    @property
    def aa(self):
        return AA3.get(self.name, "X")


def read_clean(path):
    """First model, altloc A, standard residues only -> [Residue] per chain.

    Cleaning rules, in order:
      - stop at the first ENDMDL (NMR entries hold many models)
      - keep altloc ' ' or 'A' and blank the indicator, so downstream tools
        do not see two copies of the same atom
      - map MSE/HSD/PTR... to their standard parent and rewrite the name,
        because selenomethionine is a methionine for every purpose here
      - drop waters, ligands, ions and hydrogens
    """
    op = gzip.open if path.endswith(".gz") else open
    chains = collections.OrderedDict()
    with op(path, "rt", errors="ignore") as fh:
        for ln in fh:
            if ln.startswith("ENDMDL"):
                break
            if ln[:6] not in ("ATOM  ", "HETATM"):
                continue
            alt = ln[16]
            if alt not in (" ", "A"):
                continue
            name = ln[17:20].strip()
            name = ALIAS.get(name, name)
            if name not in AA3:
                continue                       # water, ligand, ion, nucleotide
            atom = ln[12:16].strip()
            if atom.startswith("H") or (len(atom) > 1 and atom[1] == "H"
                                        and atom[0].isdigit()):
                continue                       # hydrogens
            try:
                num = int(ln[22:26])
            except ValueError:
                continue
            ch, icode = ln[21], ln[26]
            key = (ch, num, icode)
            per = chains.setdefault(ch, collections.OrderedDict())
            if key not in per:
                per[key] = Residue(ch, num, icode, name)
            # normalise: altloc blanked, residue name rewritten, element kept
            per[key].atoms.append(ln[:16] + " " + "%-3s" % name + ln[20:])
    return {c: list(v.values()) for c, v in chains.items()}


# ------------------------------------------------------------------ identify

def seq_of(residues):
    return "".join(r.aa for r in residues)


def kmer_sim(s, ref, k=4):
    """Fraction of the shorter sequence's k-mers shared with `ref`.

    Position-by-position identity fails here: an IMGT-numbered beta2m is
    offset and gapped relative to a plain reference sequence, so a direct
    comparison scored near zero and beta2m came back unclassified.  k-mer
    overlap does not care where the sequence starts.
    """
    if len(s) < k or len(ref) < k:
        return 0.0
    a = set(s[i:i + k] for i in range(len(s) - k + 1))
    b = set(ref[i:i + k] for i in range(len(ref) - k + 1))
    return len(a & b) / float(min(len(a), len(b)))


def identity(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return sum(1 for i in range(n) if a[i] == b[i]) / float(n)


def anarci_domains(chain_seqs):
    """{chain: (chain_type, [((num, icode), aa) ...], start, end)} for TCR V domains.

    chain_type is 'A' or 'B' as ANARCI reports it -- alpha or beta.  Only the
    first domain per chain is used: a TCR chain has one V domain, and the
    constant domain is not numbered by IMGT's V scheme.
    """
    try:
        from anarci import run_anarci
    except ImportError:
        sys.exit("ABORT: ANARCI is not installed.\n"
                 "       conda install -c bioconda anarci   (it needs hmmer)\n"
                 "       IMGT numbering is not derivable without it, and\n"
                 "       guessing it is how CDR3 measurements go wrong.")
    pairs = [(c, s) for c, s in chain_seqs.items() if 70 <= len(s) <= 400]
    if not pairs:
        return {}
    try:
        _seqs, numbering, details, _hits = run_anarci(
            pairs, scheme="imgt", allow=set(["A", "B", "G", "D", "H", "K", "L"]))
    except FileNotFoundError:
        sys.exit("ABORT: ANARCI is installed but hmmscan is not on PATH.\n"
                 "       conda install -c bioconda hmmer")
    out = {}
    for (cid, _s), num, det in zip(pairs, numbering, details):
        if not num or not det:
            continue
        dom, start, end = num[0]
        out[cid] = (det[0]["chain_type"], dom, start, end)
    return out


# Framework positions that DIFFER between the alpha and beta V domains,
# extracted from 1AO7 (chain D is alpha -- J region FGAGTQVVV, TRAJ; chain E is
# beta -- FGPGTRLTV, TRBJ).  62 of the 85 shared framework positions differ, so
# the two are easy to tell apart without ANARCI.  CDR positions are excluded:
# they are hypervariable and would add noise rather than signal.
REF_A = {3:'E', 5:'E', 7:'N', 8:'S', 9:'G', 10:'P', 11:'L', 12:'S', 13:'V',
         14:'P', 15:'E', 17:'A', 18:'I', 19:'A', 20:'S', 22:'N', 24:'T',
         25:'Y', 26:'S', 39:'F', 40:'F', 45:'Y', 46:'S', 48:'K', 49:'S',
         50:'P', 51:'E', 54:'M', 55:'S', 66:'K', 67:'E', 68:'D', 74:'G',
         75:'R', 76:'F', 77:'T', 78:'A', 79:'Q', 80:'L', 81:'N', 83:'A',
         84:'S', 85:'Q', 86:'Y', 87:'V', 88:'S', 90:'L', 91:'I', 92:'R',
         93:'D', 94:'S', 95:'Q', 98:'D', 99:'S', 100:'A', 101:'T', 103:'L',
         120:'A', 123:'Q', 124:'V', 125:'V', 128:'P'}
REF_B = {3:'G', 5:'T', 7:'T', 8:'P', 9:'K', 10:'F', 11:'Q', 12:'V', 13:'L',
         14:'K', 15:'T', 17:'Q', 18:'S', 19:'M', 20:'T', 22:'Q', 24:'A',
         25:'Q', 26:'D', 39:'M', 40:'S', 45:'D', 46:'P', 48:'M', 49:'G',
         50:'L', 51:'R', 54:'H', 55:'Y', 66:'T', 67:'D', 68:'Q', 74:'N',
         75:'G', 76:'Y', 77:'N', 78:'V', 79:'S', 80:'R', 81:'S', 83:'T',
         84:'T', 85:'E', 86:'D', 87:'F', 88:'P', 90:'R', 91:'L', 92:'L',
         93:'S', 94:'A', 95:'A', 98:'Q', 99:'T', 100:'S', 101:'V', 103:'F',
         120:'P', 123:'R', 124:'L', 125:'T', 128:'E'}


def imgt_anchored(residues):
    """True when a chain carries the invariant IMGT V-domain anchors."""
    g = {r.num: r.aa for r in residues if r.icode == " "}
    return g.get(23) == "C" and g.get(41) == "W" and g.get(104) == "C"


def alpha_or_beta(residues):
    """'rec1' or 'rec2' by framework identity to the 1AO7 reference pair.

    DO NOT TRUST THIS ON ITS OWN.  Measured across 222 STCRDab complexes it
    separates 1AO7 perfectly (62-0, the structure it was built from) and
    almost nothing else: median margin 0.04, and 3DXA ties exactly.  TCR V
    genes are far too diverse for one reference pair -- telling alpha from
    beta properly needs germline HMM profiles, which is what ANARCI has.

    It is kept only as a tie-break when a chain map is supplied for some
    chains and not others.  The score margin is returned so the caller can
    see how weak the call is.
    """
    g = {r.num: r.aa for r in residues if r.icode == " "}
    sa = sum(1 for p, a in REF_A.items() if g.get(p) == a)
    sb = sum(1 for p, a in REF_B.items() if g.get(p) == a)
    return ("rec1" if sa >= sb else "rec2"), sa, sb


def load_chain_map(path):
    """{pdbid: {chain: slot}} from a TCR3d / STCRDab style table.

    The database already knows which chain is alpha and which is beta, so the
    authoritative answer is to read it rather than infer it.  Any delimited
    file works: the header is scanned for a PDB column and for columns naming
    the alpha/beta (or heavy/light) chain.
    """
    if not path:
        return {}
    rows = [l.rstrip("\n") for l in open(path) if l.strip()]
    if not rows:
        return {}
    delim = "\t" if "\t" in rows[0] else ","
    head = [h.strip().strip('"').lower() for h in rows[0].split(delim)]

    def find(*keys):
        for i, h in enumerate(head):
            if any(k in h for k in keys):
                return i
        return None

    ipdb = find("pdb", "id")
    ia = find("alpha", "tcra", "heavy", "vh")
    ib = find("beta", "tcrb", "light", "vl")
    if ipdb is None or ia is None or ib is None:
        sys.exit("ABORT: --chain-map %s needs a PDB column plus alpha/beta "
                 "(or heavy/light) chain columns; found headers %s"
                 % (path, head))
    out = {}
    for r in rows[1:]:
        f = [x.strip().strip('"') for x in r.split(delim)]
        if len(f) <= max(ipdb, ia, ib):
            continue
        pid = f[ipdb].lower()
        if len(pid) != 4:
            continue
        m = {}
        if f[ia]:
            m[f[ia][0]] = "rec1"
        if f[ib]:
            m[f[ib][0]] = "rec2"
        if m:
            out[pid] = m
    return out


def classify(chains, verbose=False, numbering="anarci", v_max=128,
             cmap=None):
    """chain id -> one of rec1 / rec2 / mhc / b2m / peptide / None.

    numbering='anarci'  ANARCI identifies and numbers the V domains
    numbering='file'    the file is already IMGT-numbered (STCRDab and the
                        like).  The anchors are verified rather than assumed,
                        and alpha/beta is decided by framework identity,
                        because IMGT positions alone do not say which is which.
    """
    seqs = {c: seq_of(r) for c, r in chains.items()}
    doms = {} if numbering == "file" else anarci_domains(seqs)
    out = {}
    for c, r in chains.items():
        s, n = seqs[c], len(r)

        # TCR3d's ".axis" files carry two poly-glycine pseudo-chains marking
        # the docking axis.  They are plain ATOM records with standard GLY
        # residues, so cleaning keeps them; only their composition gives them
        # away.  They currently fall through every size test by luck -- this
        # rejects them on purpose, so a longer axis in some future file does
        # not get mistaken for a peptide.
        if n and s.count("G") / float(n) > 0.9:
            out[c] = None
            if verbose:
                print("      chain %s: %3d res, %.0f%% glycine -> marker, dropped"
                      % (c, n, 100.0 * s.count("G") / n))
            continue

        if numbering == "file" and imgt_anchored(r):
            if cmap and c in cmap:
                out[c] = cmap[c]            # the database said so
                if verbose:
                    print("      chain %s: %3d res -> %s (from chain map)"
                          % (c, n, out[c]))
            else:
                kind, sa, sb = alpha_or_beta([x for x in r if x.num <= v_max])
                out[c] = kind
                if verbose:
                    print("      chain %s: %3d res -> %s (WEAK: alpha %d / "
                          "beta %d -- supply --chain-map)" % (c, n, kind, sa, sb))
            continue
        if c in doms:
            out[c] = CHAIN_SLOT.get(doms[c][0])
            continue
        if n <= 20:
            out[c] = "peptide"
        elif 85 <= n <= 115 and kmer_sim(s, B2M) > 0.25:
            out[c] = "b2m"
        elif n >= 150:
            out[c] = "mhc"
        else:
            out[c] = None
        if verbose:
            print("      chain %s: %3d res -> %s" % (c, n, out[c]))
    return out, doms


# ------------------------------------------------------------------ copies

def ca_coords(residues):
    import numpy as np
    xyz = []
    for r in residues:
        for a in r.atoms:
            if a[12:16].strip() == "CA":
                xyz.append((float(a[30:38]), float(a[38:46]), float(a[46:54])))
    return np.asarray(xyz, float) if xyz else None


def contacts(p, q, cut=12.0):
    """Number of CA pairs within `cut` -- a cheap proxy for interface size."""
    import numpy as np
    if p is None or q is None or not len(p) or not len(q):
        return 0
    d = np.linalg.norm(p[:, None, :] - q[None, :, :], axis=2)
    return int((d < cut).sum())


def pick_unit(chains, kinds, want, copy_index=0, verbose=False):
    """Assemble one biological complex outward from a peptide.

    WHY NOT CLUSTER BY DISTANCE
      Crystal copies pack closer than any threshold that still holds a single
      complex together -- on 1LP9 a 60 A centroid cut merged both copies into
      one ten-chain blob.  Interfaces are what define a complex, so the unit
      is built by following them:

          peptide -> the MHC it sits in       (most CA contacts)
                  -> the b2m on that MHC
                  -> the alpha and beta V domains docked on that MHC

      Each step takes the best-contacting partner, so a second copy elsewhere
      in the cell is never reached.  Copies are enumerated by peptide, which
      is the one chain guaranteed to appear exactly once per complex.
    """
    coords = {c: ca_coords(r) for c, r in chains.items() if kinds.get(c)}
    mhcs = sorted(c for c in coords if kinds[c] == "mhc")
    if not mhcs:
        return None, "no MHC chain found"

    # GLOBAL ASSIGNMENT, NOT GREEDY PICKING.
    # Taking each MHC's best-contacting partner in turn goes wrong whenever a
    # chain sits between two complexes.  On 2AK4, TCR chain I touches MHC A
    # (171 CA pairs) more than A's true partner D does (104), because I also
    # straddles F; greedy therefore gave A<-I/J and left F with no TCR.
    # Scoring every MHC-partner pair and solving the assignment that maximises
    # the TOTAL recovers A<-D/E and F<-I/J (319 contacts against 238).
    units = [{"mhc": m} for m in mhcs]

    def assign(kind):
        cands = sorted(c for c in coords if kinds[c] == kind)
        if not cands:
            return
        M = [[contacts(coords[m], coords[c]) for c in cands] for m in mhcs]
        rows, cols = None, None
        try:
            from scipy.optimize import linear_sum_assignment
            import numpy as np
            rows, cols = linear_sum_assignment(-np.asarray(M, float))
        except ImportError:
            taken, rows, cols = set(), [], []
            order = sorted(((M[i][j], i, j)
                            for i in range(len(mhcs))
                            for j in range(len(cands))), reverse=True)
            for _s, i, j in order:
                if i in rows or j in taken:
                    continue
                rows.append(i); cols.append(j); taken.add(j)
        for i, j in zip(rows, cols):
            if M[i][j] <= 0:
                continue                      # not in contact: leave unfilled
            units[i][kind] = cands[j]
            if verbose:
                print("      %-7s %s -> mhc %s  (%d CA contacts)"
                      % (kind, cands[j], mhcs[i], M[i][j]))

    for kind in ("peptide", "b2m", "rec1", "rec2"):
        assign(kind)

    complete = [u for u in units if all(k in u for k in want)]
    if not complete:
        have = units[0] if units else {}
        return None, ("no complete copy; best had %s"
                      % ("+".join(sorted(have)) or "nothing"))
    if copy_index >= len(complete):
        return None, "only %d complete copy/copies present" % len(complete)
    if verbose:
        print("      copy %d of %d = %s" % (copy_index, len(complete),
              " ".join("%s:%s" % (k, v)
                       for k, v in sorted(complete[copy_index].items()))))
    return complete[copy_index], ""


# ------------------------------------------------------------------ renumber

def renumber_tcr(residues, dom):
    """Apply ANARCI's IMGT numbering to a V domain.

    ANARCI numbers the SEQUENCE, so its output is walked alongside the residue
    list: alignment gaps ('-') consume a numbering slot but no residue.
    Residues outside the domain (the constant domain, tags) are dropped --
    this pipeline docks V domains.
    """
    _ct, numbering, start, end = dom
    span = residues[start:end + 1]
    out, k = [], 0
    for (pos, icode), aa in numbering:
        if aa == "-":
            continue
        if k >= len(span):
            break
        r = span[k]
        k += 1
        if r.aa != aa:
            return None, ("ANARCI/residue mismatch at IMGT %d%s: %s vs %s"
                          % (pos, icode.strip(), r.aa, aa))
        r.num, r.icode = pos, (icode if icode.strip() else " ")
        out.append(r)
    if not out:
        return None, "empty domain"
    return out, ""


def renumber_serial(residues, start=1):
    for i, r in enumerate(residues):
        r.num, r.icode = start + i, " "
    return residues


# ------------------------------------------------------------------ write

def emit(residues, chain_id, band, serial_start=1):
    """Rewrite residues onto one chain with a band offset applied."""
    lines, n = [], serial_start
    for r in residues:
        for a in r.atoms:
            lines.append("%s%5d%s%s%4d%s%s"
                         % (a[:6], n, a[11:21], chain_id,
                            r.num + band, r.icode, a[27:].rstrip("\n")) + "\n")
            n += 1
    return lines, n


def validate(residues, label):
    """Refuse silently-wrong output: duplicate keys, missing IMGT anchors."""
    keys = [(r.num, r.icode) for r in residues]
    if len(keys) != len(set(keys)):
        dup = [k for k, c in collections.Counter(keys).items() if c > 1]
        return "%s has duplicate residue keys %s" % (label, dup[:4])
    if label in ("rec1", "rec2"):
        by = {(r.num, r.icode): r.aa for r in residues}
        want = {23: "C", 41: "W", 104: "C"}
        bad = {p: by.get((p, " "), "-") for p, a in want.items()
               if by.get((p, " ")) != a}
        if bad:
            return ("%s missing IMGT anchors (got %s) -- numbering is not IMGT"
                    % (label, bad))
        if not any(105 <= r.num <= 117 for r in residues):
            return "%s has no residue in the CDR3 window 105-117" % label
    return ""


# ------------------------------------------------------------------ main

def read_ids(spec):
    """PDB IDs from a plain list, a comma-separated string, or a CSV/TSV.

    TCR3d and STCRDab both hand you a table rather than a list, so rather than
    making you cut a column out first, any delimited file is scanned for the
    column that most looks like PDB IDs: four characters, starting with a
    digit, unique.  A header row is skipped automatically.
    """
    if not os.path.exists(spec):
        return [x.strip() for x in spec.split(",") if x.strip()]
    rows = [l.rstrip("\n") for l in open(spec) if l.strip()
            and not l.startswith("#")]
    if not rows:
        return []
    delim = "\t" if "\t" in rows[0] else ("," if "," in rows[0] else None)
    if delim is None:
        return [r.split()[0] for r in rows]
    table = [r.split(delim) for r in rows]
    ncol = max(len(r) for r in table)
    best, best_score = None, 0
    for j in range(ncol):
        col = [r[j].strip().strip('"') for r in table if len(r) > j]
        hits = [v for v in col
                if len(v) == 4 and v[0].isdigit() and v.isalnum()]
        if len(hits) > best_score:
            best, best_score = j, len(hits)
    if best is None or best_score == 0:
        return []
    out, seen = [], set()
    for r in table:
        if len(r) <= best:
            continue
        v = r[best].strip().strip('"')
        if len(v) == 4 and v[0].isdigit() and v.isalnum() and v.lower() not in seen:
            seen.add(v.lower())
            out.append(v)
    return out


def process(pdbid, args):
    if args.from_dir:
        path = os.path.join(args.from_dir, "%s.pdb" % pdbid.lower())
        if not os.path.exists(path):
            return None, "not found in %s" % args.from_dir
    else:
        path = fetch(pdbid, args.cache)
    if not path:
        return None, "download failed"
    chains = read_clean(path)
    if not chains:
        return None, "no standard residues after cleaning"
    kinds, doms = classify(chains, args.verbose, args.numbering,
                           args.v_max, CHAIN_MAP.get(pdbid.lower()))

    want = ["rec1", "rec2", "mhc", "peptide"]
    parts, why = pick_unit(chains, kinds, want, args.copy, args.verbose)
    if parts is None:
        return None, why

    # rec1 and rec2 are assigned to MHCs independently, so with partial chain
    # labels they can be drawn from different copies -- on 3DXA that produced
    # a "complex" pairing chain D with chain O.  The two receptor domains of a
    # real Fv pack against each other with hundreds of CA pairs, so a weak
    # rec1-rec2 interface means the pairing is wrong, not merely unusual.
    pair_contacts = contacts(ca_coords(chains[parts["rec1"]]),
                             ca_coords(chains[parts["rec2"]]))
    if pair_contacts < 50:
        return None, ("rec1 %s and rec2 %s share only %d CA contacts -- they "
                      "are not a paired Fv (supply --chain-map)"
                      % (parts["rec1"], parts["rec2"], pair_contacts))

    # what kind of receptor this turned out to be, for the manifest
    cts = tuple(doms[parts[k]][0] for k in ("rec1", "rec2")
                if parts.get(k) in doms)
    receptor = RECEPTOR.get(cts, "-".join(cts) if cts else "from-file")

    built = {}
    for k in ("rec1", "rec2"):
        if args.numbering == "file":
            # already IMGT-numbered: keep the V domain, verify, do not renumber
            res = [r for r in chains[parts[k]] if r.num <= args.v_max]
            err = ""
        else:
            res, err = renumber_tcr(chains[parts[k]], doms[parts[k]])
        if res is None:
            return None, "%s: %s" % (k, err)
        e = validate(res, k)
        if e:
            return None, e
        built[k] = res
    for k in ("mhc", "b2m", "peptide"):
        if k not in parts:
            continue
        res = chains[parts[k]]
        if k == "mhc" and args.mhc_groove:
            res = res[:args.mhc_groove]
        built[k] = renumber_serial(res)
        e = validate(built[k], k)
        if e:
            return None, e

    os.makedirs(args.out, exist_ok=True)
    tcr_lines, n = [], 1
    for k in ("rec1", "rec2"):
        ls, n = emit(built[k], "A", BANDS[k], n)
        tcr_lines += ls
    tcr_lines.append("TER\n")
    pm_lines, n = [], 1
    for k in ("mhc", "b2m", "peptide"):
        if k not in built:
            continue
        ls, n = emit(built[k], "B", BANDS[k], n)
        pm_lines += ls
    pm_lines.append("TER\n")

    base = os.path.join(args.out, pdbid.lower())
    open(base + "_tcr.pdb", "w").writelines(tcr_lines + ["END\n"])
    open(base + "_pmhc.pdb", "w").writelines(pm_lines + ["END\n"])
    open(base + "_complex.pdb", "w").writelines(tcr_lines + pm_lines + ["END\n"])
    return {
        "chains": "".join(parts[k] for k in sorted(parts)),
        "receptor": receptor,
        "rec1": len(built["rec1"]), "rec2": len(built["rec2"]),
        "mhc": len(built.get("mhc", [])), "b2m": len(built.get("b2m", [])),
        "peptide": len(built.get("peptide", [])),
        "cdr3a": sum(1 for r in built["rec1"] if 105 <= r.num <= 117),
        "cdr3b": sum(1 for r in built["rec2"] if 105 <= r.num <= 117),
        "icodes": sum(1 for k in ("rec1", "rec2")
                      for r in built[k] if r.icode != " "),
    }, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True,
                    help="file with one PDB ID per line, or a comma-separated list")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default="pdb_cache")
    ap.add_argument("--copy", type=int, default=0,
                    help="which complete complex to take when the asymmetric "
                         "unit holds several (0 = first)")
    ap.add_argument("--mhc-groove", type=int, default=0, metavar="N",
                    help="truncate the MHC heavy chain to its first N residues "
                         "(180 keeps alpha1+alpha2, the groove the TCR sees)")
    ap.add_argument("--numbering", default="anarci", choices=["anarci", "file"],
                    help="anarci: number the V domains with ANARCI (needs "
                         "hmmscan).  file: the input is already IMGT-numbered "
                         "(STCRDab, TCR3d) -- the anchors are verified, not "
                         "assumed, and nothing is renumbered")
    ap.add_argument("--from-dir", default="",
                    help="read <id>.pdb from this directory instead of "
                         "downloading from RCSB")
    ap.add_argument("--v-max", type=int, default=128,
                    help="with --numbering file, keep IMGT positions <= this "
                         "(128 = the V domain; drops the constant domain)")
    ap.add_argument("--chain-map", default="",
                    help="table naming the alpha/beta (or heavy/light) chain "
                         "per PDB entry -- e.g. the CSV TCR3d hands you. "
                         "Strongly recommended with --numbering file: chain "
                         "type cannot be inferred from sequence reliably")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    global CHAIN_MAP
    CHAIN_MAP = load_chain_map(args.chain_map)
    if CHAIN_MAP:
        print("chain map loaded for %d entries" % len(CHAIN_MAP))
    elif args.numbering == "file":
        print("WARNING: --numbering file without --chain-map; alpha/beta is "
              "being guessed from one reference pair and is unreliable")
    ids = read_ids(args.ids)
    if not ids:
        sys.exit("no PDB IDs found in %s" % args.ids)
    print("%d PDB ID(s) to prepare" % len(ids))

    os.makedirs(args.out, exist_ok=True)
    man = open(os.path.join(args.out, "manifest.tsv"), "w")
    man.write("pdbid\treceptor\tchains\trec1\trec2\tmhc\tb2m\tpeptide"
              "\tcdr3_1\tcdr3_2\ticodes\n")
    ok, refused = 0, []
    for i, pid in enumerate(ids, 1):
        if args.verbose:
            print("[%d/%d] %s" % (i, len(ids), pid))
        try:
            info, why = process(pid, args)
        except Exception as e:
            info, why = None, "%s: %s" % (type(e).__name__, e)
        if info is None:
            refused.append((pid, why))
            if args.verbose:
                print("      REFUSED: %s" % why)
            continue
        man.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n"
                  % (pid.lower(), info["receptor"], info["chains"],
                     info["rec1"], info["rec2"],
                     info["mhc"], info["b2m"], info["peptide"],
                     info["cdr3a"], info["cdr3b"], info["icodes"]))
        ok += 1
    man.close()

    print("\nPREPARED %d / %d targets -> %s" % (ok, len(ids), args.out))
    if refused:
        print("refused %d:" % len(refused))
        seen = collections.Counter(w for _p, w in refused)
        for why, k in seen.most_common(8):
            print("   %3d  %s" % (k, why))
        for p, w in refused[:6]:
            print("      e.g. %s: %s" % (p, w))
    if not ok:
        sys.exit("ABORT: nothing prepared")


if __name__ == "__main__":
    main()
