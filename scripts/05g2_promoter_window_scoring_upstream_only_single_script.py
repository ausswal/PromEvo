#!/usr/bin/env python3
"""
STAGE 05g -- Promoter-window scoring for UPSTREAM-ONLY promoter FASTA.

INPUT
-----
outputs/promoter_sequences.fasta from 01e_new9_fetch_family_and_promoters.py,
where every record is the N bp (e.g. 1000 bp) UPSTREAM of the ANNOTATED GENE
START, with no downstream flank. The last base of each record is therefore the
base immediately before the annotated gene start.

IMPORTANT: THIS IS NOT EXPERIMENTAL TSS PREDICTION
--------------------------------------------------
The reference point is the gene start taken from the annotation (GFF). It may
be the start codon or the start of the 5' UTR depending on the annotation, and
it is not an experimentally mapped TSS. This stage learns what the sequence
just upstream of that annotated position looks like and asks where inside each
promoter such sequence is found. Coordinates are relative to the annotated gene
start: -1 = last base upstream. Do not describe the output as "TSS prediction"
unless it has been compared with experimental TSS data (TSS-seq / 5'-end RNA-seq).

WHAT THIS SCRIPT DOES
---------------------
1. Positive example per gene = the WINDOW_SIZE bp immediately upstream of the
   annotated gene start. Negatives = random windows far upstream, "hard"
   windows 50-300 bp upstream of the positive, and a DINUCLEOTIDE-PRESERVING
   shuffled copy of the positive (Altschul-Erickson), so the control keeps the
   dinucleotide structure that dominates chance hits in AT-rich sequence.
2. Two scorers, blended 50/50 after min-max scaling inside each gene: a
   Position Weight Matrix (log-odds) and an MLP on 3/4/5-mer frequencies.
3. Cluster-aware K-fold cross-validation (similar promoters are never split
   between train and test). The whole held-out promoter is scanned and we test
   whether the top-scoring window is the true proximal window, compared with
   chance, for the combined scorer AND each scorer on its own (PWM only, MLP
   only), with an approximate one-sided p-value against chance.
4. Excel export: one global workbook and one workbook per genus. EVERY gene
   gets a row. By default (EXPORT_OUT_OF_FOLD = True) each gene is scored by
   the model from the fold in which it was HELD OUT, so no exported score comes
   from a model that was trained on that gene. Set it to False to score with the
   final all-genes model instead (optimistic; kept only for comparison).

HOW TO READ THE OUTPUT
----------------------
* "Relative score" is min-max scaled WITHIN each gene, so it is a rank, not a
  probability, and every gene will have a best window near 1.0. It cannot be
  used as confidence. The "MLP probability" column is the absolute classifier
  output for the same window; look at both.
* The best window tends to be the one nearest the annotated gene start, because
  every positive example came from there. A real signal needs the CV hit rate
  to be clearly above the printed chance rate AND above what a single scorer
  alone (e.g. composition captured by the MLP) already achieves.
* The CV p-value assumes genes are independent. Related promoters that survive
  the clustering threshold make it optimistic, so treat it as a sanity check.

Requirements:
    pip install biopython scikit-learn numpy pandas openpyxl
"""

import argparse
import math
import os
import re
import sys
import random
import warnings
from itertools import product
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)

# ============================== CONFIGURATION ==============================
FASTA_PATH = "outputs/promoter_sequences.fasta"
LABELS_PATH = "outputs/family_labels.csv"           # needs seq_id, organism (gene optional)
OUTPUT_DIR = "outputs/promoter_window_scoring"       # global workbook goes here
ORG_ROOT = "outputs/promoters_by_organism"           # <Genus>/<Genus>_promoter_window_scoring_summary.xlsx

EXPECTED_PROMOTER_LEN = 1000   # only used for a sanity warning, not for the maths
MAX_AMBIGUOUS_FRACTION = 0.02  # keep a sequence if <= 2% of its bases are N/ambiguous
MIN_GENES_FOR_TRAINING = 10

WINDOW_SIZE = 100              # bp; even number
SCAN_STEP = 5                  # bp between scanned window centres
N_NEG_PER_POS = 4              # random far decoys per positive
MIN_NEG_DISTANCE = 150         # decoy centres >= this far from the positive's centre
N_HARD_NEG_PER_POS = 4         # near-miss decoys per positive
HARD_NEG_MIN_OFFSET = WINDOW_SIZE // 2
HARD_NEG_MAX_OFFSET = WINDOW_SIZE * 3
KMER_KS = (3, 4, 5)
PWM_PSEUDOCOUNT = 0.5
COMBINE_WEIGHT_PWM = 0.5       # 0 = pure MLP, 1 = pure PWM
SIM_KMER_K = 12                # homology clustering (keeps paralogs/orthologs on one side)
SIM_THRESHOLD = 0.70
N_FOLDS = 5
PROMOTER_SCORE_THRESHOLD = 0.80
MIN_GENES_FOR_GENUS_STAT = 5   # per-genus CV stat is flagged unreliable below this
RANDOM_SEED = 42
EXPORT_OUT_OF_FOLD = True      # score each gene with the model that did NOT see it
# ===========================================================================

random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

HALF = WINDOW_SIZE // 2
BASES = "ACGT"
BASE_IDX = {b: i for i, b in enumerate(BASES)}
_LUT = np.full(256, -1, dtype=np.int16)
for _b, _i in BASE_IDX.items():
    _LUT[ord(_b)] = _i

# k-mer feature layout: one block per k, each block normalised separately
KMER_INDEX, KMER_BLOCKS = {}, {}
for _k in KMER_KS:
    _start = len(KMER_INDEX)
    for _p in product(BASES, repeat=_k):
        KMER_INDEX["".join(_p)] = len(KMER_INDEX)
    KMER_BLOCKS[_k] = (_start, len(KMER_INDEX))
N_FEATURES = len(KMER_INDEX)


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_promoters(fasta_path):
    """Returns (records, dropped). records = [(seq_id, SEQ)], all upper-case.
    A sequence is dropped only if it is too short to build windows, or has
    more than MAX_AMBIGUOUS_FRACTION ambiguous bases. (earlier versions dropped a
    sequence for even ONE 'N'.) Every drop is reported with its reason."""
    records, dropped = [], []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper().replace("U", "T")
        if len(seq) < 3 * WINDOW_SIZE:
            dropped.append((rec.id, f"too short ({len(seq)} bp)"))
            continue
        n_bad = sum(1 for c in seq if c not in BASES)
        if n_bad / len(seq) > MAX_AMBIGUOUS_FRACTION:
            dropped.append((rec.id, f"{n_bad / len(seq) * 100:.1f}% ambiguous bases"))
            continue
        records.append((rec.id, seq))
    return records, dropped


def extract_genus(organism):
    organism = (organism or "").strip()
    if not organism or organism.lower() == "nan":
        return "Unknown_genus"
    return organism.split()[0]


def load_label_maps(labels_path):
    """Returns ({seq_id: genus}, {seq_id: gene display name})."""
    genus_map, gene_names = {}, {}
    if not os.path.exists(labels_path):
        print(f"[!] {labels_path} not found -- all genes grouped as 'Unknown_genus'.")
        return genus_map, gene_names
    df = pd.read_csv(labels_path)
    for _, row in df.iterrows():
        sid = str(row["seq_id"])
        genus_map[sid] = extract_genus(str(row.get("organism", "")))
        g = str(row.get("gene", "")).strip()
        gene_names[sid] = g if g and g.lower() != "nan" else sid
    return genus_map, gene_names


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name).strip()).strip("_") or "Unknown"


# --------------------------------------------------------------------------
# Windows, PWM, features
# --------------------------------------------------------------------------
def extract_window(seq, center, size=WINDOW_SIZE):
    start = center - size // 2
    end = start + size
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


def true_window_center(seq):
    """Centre of the window that ENDS at the annotated gene start (= last base of the record)."""
    return len(seq) - HALF


def dinuc_shuffle(seq):
    """Dinucleotide-preserving shuffle (Altschul & Erickson 1985): a random
    Eulerian path through the dinucleotide graph. Same mono- AND dinucleotide
    counts, same first and last base. Uses the module-level `random` state."""
    if len(seq) < 4:
        return seq
    edges = defaultdict(list)
    for x, y in zip(seq, seq[1:]):
        edges[x].append(y)
    last = seq[-1]
    while True:
        shuffled = {v: random.sample(nx, len(nx)) for v, nx in edges.items()}
        last_edge = {v: nx[-1] for v, nx in shuffled.items() if v != last}
        ok = True
        for v in last_edge:                          # last edges must form a tree rooted at `last`
            seen, cur = set(), v
            while cur != last:
                if cur in seen or cur not in last_edge:
                    ok = False
                    break
                seen.add(cur)
                cur = last_edge[cur]
            if not ok:
                break
        if ok:
            break
    ptr = defaultdict(int)
    out, cur = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        nxt = shuffled[cur][ptr[cur]]
        ptr[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def build_training_windows(records):
    positives, negatives = [], []
    for _, seq in records:
        pc = true_window_center(seq)
        pos_win = extract_window(seq, pc)
        if pos_win is None:
            continue
        positives.append(pos_win)

        tries = collected = 0                       # far random decoys
        while collected < N_NEG_PER_POS and tries < 50:
            tries += 1
            c = random.randint(HALF, len(seq) - HALF)
            if abs(c - pc) < MIN_NEG_DISTANCE:
                continue
            w = extract_window(seq, c)
            if w:
                negatives.append(w)
                collected += 1

        tries = collected = 0                       # hard near-miss decoys (upstream only)
        while collected < N_HARD_NEG_PER_POS and tries < 50:
            tries += 1
            c = pc - random.randint(HARD_NEG_MIN_OFFSET, HARD_NEG_MAX_OFFSET)
            if c < HALF:
                continue
            w = extract_window(seq, c)
            if w:
                negatives.append(w)
                collected += 1

        negatives.append(dinuc_shuffle(pos_win))    # dinucleotide-preserving shuffled copy
    return positives, negatives


def background_freqs(seqs):
    counts = {b: 1 for b in BASES}
    for s in seqs:
        for b in BASES:
            counts[b] += s.count(b)
    tot = sum(counts.values())
    return {b: counts[b] / tot for b in BASES}


def build_pwm(positives, bg):
    counts = np.full((4, WINDOW_SIZE), PWM_PSEUDOCOUNT)
    for win in positives:
        idx = _LUT[np.frombuffer(win.encode("ascii"), dtype=np.uint8)]
        for pos, b in enumerate(idx):
            if b >= 0:
                counts[b, pos] += 1
    freqs = counts / counts.sum(axis=0)
    bgv = np.array([bg[b] for b in BASES]).reshape(4, 1)
    return np.log2(freqs / bgv)


def kmer_features(win):
    vec = np.zeros(N_FEATURES)
    for k in KMER_KS:
        s, e = KMER_BLOCKS[k]
        n = 0
        for i in range(len(win) - k + 1):
            j = KMER_INDEX.get(win[i:i + k])
            if j is not None:
                vec[j] += 1
                n += 1
        if n:
            vec[s:e] /= n
    return vec


def train_model(records):
    positives, negatives = build_training_windows(records)
    if len(positives) < MIN_GENES_FOR_TRAINING:
        raise RuntimeError(f"only {len(positives)} positive windows could be built")
    pwm = build_pwm(positives, background_freqs([s for _, s in records]))
    X = np.array([kmer_features(w) for w in positives + negatives])
    y = np.array([1] * len(positives) + [0] * len(negatives))
    scaler = StandardScaler()
    clf = MLPClassifier(hidden_layer_sizes=(64, 32), activation="relu",
                        max_iter=3000, random_state=RANDOM_SEED)
    clf.fit(scaler.fit_transform(X), y)
    return pwm, clf, scaler


def _scan_raw(seq, pwm, clf, scaler):
    """Raw per-window scores across the WHOLE promoter.
    Returns (centres, pwm_log_odds, mlp_probability)."""
    centers = list(range(HALF, len(seq) - HALF + 1, SCAN_STEP))
    if centers[-1] != len(seq) - HALF:              # always include the window ending at the gene start
        centers.append(len(seq) - HALF)
    idx = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    cols = np.arange(WINDOW_SIZE)

    pwm_scores, wins = [], []
    for c in centers:
        s0 = c - HALF
        sl = idx[s0:s0 + WINDOW_SIZE]
        ok = sl >= 0
        pwm_scores.append(float(pwm[sl[ok], cols[ok]].sum()))
        wins.append(seq[s0:s0 + WINDOW_SIZE])

    X = np.array([kmer_features(w) for w in wins])
    mlp = clf.predict_proba(scaler.transform(X))[:, 1]
    return np.array(centers), np.array(pwm_scores), mlp


def _norm(a):
    a = np.asarray(a, dtype=float)
    return np.zeros_like(a) if a.max() == a.min() else (a - a.min()) / (a.max() - a.min())


def _combine(pwm_scores, mlp):
    return COMBINE_WEIGHT_PWM * _norm(pwm_scores) + (1 - COMBINE_WEIGHT_PWM) * _norm(mlp)


def scan_sequence(seq, pwm, clf, scaler):
    """Returns (centres, combined_relative_score, mlp_probability)."""
    centers, pwm_scores, mlp = _scan_raw(seq, pwm, clf, scaler)
    return centers, _combine(pwm_scores, mlp), mlp


def call_promoters(centers, combined, threshold):
    """Merge neighbouring windows scoring >= threshold into discrete calls
    (0-based inclusive start/end)."""
    calls, cur = [], None
    for pos, sc in zip(centers, combined):
        if sc >= threshold:
            ws, we = pos - HALF, pos + HALF - 1
            if cur is None:
                cur = {"start": ws, "end": we, "peak_score": sc, "peak_center": pos}
            else:
                cur["end"] = max(cur["end"], we)
                if sc > cur["peak_score"]:
                    cur["peak_score"], cur["peak_center"] = sc, pos
        elif cur is not None:
            calls.append(cur)
            cur = None
    if cur is not None:
        calls.append(cur)
    return calls


# --------------------------------------------------------------------------
# Homology-aware folds
# --------------------------------------------------------------------------
class _DSU:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def cluster_sequences(records, k=SIM_KMER_K, threshold=SIM_THRESHOLD):
    n = len(records)
    sets = [set(s[i:i + k] for i in range(len(s) - k + 1)) for _, s in records]
    dsu = _DSU(n)
    for i in range(n):
        for j in range(i + 1, n):
            u = len(sets[i] | sets[j])
            if u and len(sets[i] & sets[j]) / u >= threshold:
                dsu.union(i, j)
    groups = defaultdict(list)
    for i in range(n):
        groups[dsu.find(i)].append(i)
    return list(groups.values())


def make_folds(records, n_folds):
    clusters = cluster_sequences(records)
    random.shuffle(clusters)
    clusters.sort(key=len, reverse=True)
    n_folds = max(2, min(n_folds, len(clusters)))
    folds = [[] for _ in range(n_folds)]
    for cl in clusters:                              # biggest cluster -> currently smallest fold
        min(folds, key=len).extend(records[i] for i in cl)
    return folds, len(clusters)


# --------------------------------------------------------------------------
# Cross-validation (honest accuracy) with per-genus breakdown
# --------------------------------------------------------------------------
def _hit_stats(hits, chance):
    """Observed vs expected hit rate. Expected count = sum of per-gene chance
    probabilities (Poisson-binomial); one-sided p from a normal approximation."""
    n = len(hits)
    obs, exp = float(sum(hits)), float(sum(chance))
    var = float(sum(c * (1 - c) for c in chance))
    z = (obs - exp) / math.sqrt(var) if var > 0 else float("nan")
    p = 0.5 * math.erfc(z / math.sqrt(2)) if z == z else float("nan")
    return obs / n, exp / n, z, p


def run_cross_validation(records, genus_map):
    """Cluster-aware CV. Returns {seq_id: (pwm, clf, scaler)} where each model is
    the one from the fold in which that gene was held out (used for the
    out-of-fold export)."""
    folds, n_clusters = make_folds(records, N_FOLDS)
    print(f"\n=== {len(folds)}-fold cluster-aware cross-validation "
          f"({n_clusters} similarity clusters among {len(records)} genes) ===")
    if n_clusters >= 0.95 * len(records):
        print(f"[!] Almost every gene is its own cluster (SIM_KMER_K={SIM_KMER_K}, "
              f"SIM_THRESHOLD={SIM_THRESHOLD}). The clustering is merging almost nothing, so related "
              f"promoters may still fall on both sides of a split. Consider a lower SIM_THRESHOLD.")

    scorers = ("combined", "PWM only", "MLP only")
    errs = {k: [] for k in scorers}
    hits = {k: [] for k in scorers}
    chance, base_errs, per_genus = [], [], defaultdict(list)
    fold_models = {}

    for fi, test in enumerate(folds):
        train = [r for j, f in enumerate(folds) if j != fi for r in f]
        if not test or len(train) < MIN_GENES_FOR_TRAINING:
            continue
        pwm, clf, scaler = train_model(train)
        fold_errs = []
        for sid, seq in test:
            fold_models[sid] = (pwm, clf, scaler)
            centers, pwm_scores, mlp = _scan_raw(seq, pwm, clf, scaler)
            truth = true_window_center(seq)
            dist = np.abs(centers - truth)
            score_sets = {"combined": _combine(pwm_scores, mlp),
                          "PWM only": _norm(pwm_scores),
                          "MLP only": _norm(mlp)}
            for k, sc in score_sets.items():
                e = int(dist[int(np.argmax(sc))])
                errs[k].append(e)
                hits[k].append(e <= WINDOW_SIZE)
            fold_errs.append(errs["combined"][-1])
            base_errs.append(float(dist.mean()))           # expected error of a random pick
            chance.append(float((dist <= WINDOW_SIZE).mean()))
            per_genus[genus_map.get(sid, "Unknown_genus")].append(
                (errs["combined"][-1], hits["combined"][-1]))
        print(f"Fold {fi + 1}/{len(folds)}: {len(test)} genes -> median abs error "
              f"{np.median(fold_errs):.0f} bp")

    if errs["combined"]:
        n = len(errs["combined"])
        print(f"\nTop-scoring window vs the true proximal window (n={n} held-out genes):")
        print(f"  median abs error (combined): {np.median(errs['combined']):.0f} bp   "
              f"(random-pick baseline ~{np.median(base_errs):.0f} bp)")
        print(f"  hit = top window within {WINDOW_SIZE} bp of the true proximal window")
        for k in scorers:
            rate, exp, z, p = _hit_stats(hits[k], chance)
            print(f"  {k:<9}: hit {rate * 100:5.1f}%  chance {exp * 100:5.1f}%  z={z:5.2f}  p~{p:.2g}")
        print("  -> trust the exported scores only if the hit rate is clearly above chance, and check")
        print("     whether PWM only / MLP only already explain it (simple composition near the gene start).")
    print("\n--- per-genus (same global model, NOT trained per genus) ---")
    for g, v in sorted(per_genus.items(), key=lambda kv: -len(kv[1])):
        if len(v) < MIN_GENES_FOR_GENUS_STAT:
            print(f"  {g}: only {len(v)} held-out gene(s) -- too few for a stable estimate.")
        else:
            print(f"  {g}: n={len(v)}, median abs error {np.median([e for e, _ in v]):.0f} bp, "
                  f"within {WINDOW_SIZE} bp: {np.mean([h for _, h in v]) * 100:.0f}%")
    print("=== end cross-validation ===\n")
    return fold_models


# --------------------------------------------------------------------------
# Excel export: one global workbook + one workbook per genus
# --------------------------------------------------------------------------
def unique_sheet_name(base, used):
    name = "".join(ch for ch in base if ch.isalnum() or ch in "_- ")[:31] or "gene"
    cand, i = name, 2
    while cand.lower() in used or cand.lower() == "summary":
        suf = f"_{i}"
        cand = name[:31 - len(suf)] + suf
        i += 1
    used.add(cand.lower())
    return cand


def describe_gene(seq_id, seq, gene, genus, model, scored_by, score_col):
    pwm, clf, scaler = model
    L = len(seq)                                     # annotated gene start = boundary after the last base
    centers, combined, mlp = scan_sequence(seq, pwm, clf, scaler)
    calls = call_promoters(centers, combined, PROMOTER_SCORE_THRESHOLD)

    bi = int(np.argmax(combined))                    # ALWAYS report the best window
    bc, bs = int(centers[bi]), float(combined[bi])
    row = {
        "Gene": gene, "Genus": genus,
        score_col: len(calls),
        "Best window (1-based position in record)": f"{bc - HALF + 1}-{bc + HALF}",
        "Relative score": round(bs, 2),
        "MLP probability": round(float(mlp[bi]), 3),
        "Best window location (relative to annotated gene start)":
            f"{(bc - HALF) - L}\u2192{(bc + HALF - 1) - L}",
        "Relative score class": f">= {PROMOTER_SCORE_THRESHOLD:.2f}" if bs >= PROMOTER_SCORE_THRESHOLD
                                else f"< {PROMOTER_SCORE_THRESHOLD:.2f}",
        "Scored by": scored_by,
    }
    detail = None
    if calls:
        rows = []
        for c in sorted(calls, key=lambda c: c["peak_score"], reverse=True):
            rows.append({"Start-End (local)": f"{c['start'] + 1}-{c['end'] + 1}",
                         "Relative score": round(float(c["peak_score"]), 2),
                         "Location (relative to annotated gene start)":
                             f"{c['start'] - L}\u2192{c['end'] - L}"})
        detail = pd.DataFrame(rows)
    return row, detail


def write_workbook(path, summary_rows, sheets, drop_genus_col):
    df = pd.DataFrame(summary_rows)
    if drop_genus_col and "Genus" in df.columns:
        df = df.drop(columns=["Genus"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Summary", index=False)
        for name, d in sheets.items():
            d.to_excel(writer, sheet_name=name, index=False)


def export_all(records, genus_map, gene_names, final_model, fold_models):
    score_col = f"#Window calls (relative score>={PROMOTER_SCORE_THRESHOLD:.2f})"
    all_rows, all_sheets, used_all = [], {}, set()
    by_genus_rows, by_genus_sheets, used_by_genus = defaultdict(list), defaultdict(dict), defaultdict(set)
    n_fallback = 0

    for sid, seq in records:
        genus = genus_map.get(sid, "Unknown_genus")
        gene = gene_names.get(sid, sid)
        if EXPORT_OUT_OF_FOLD and sid in fold_models:
            model, scored_by = fold_models[sid], "out-of-fold model"
        else:
            model = final_model
            scored_by = "final all-genes model (optimistic)"
            n_fallback += 1
        row, detail = describe_gene(sid, seq, gene, genus, model, scored_by, score_col)
        all_rows.append(row)
        by_genus_rows[genus].append(row)
        if detail is not None:
            all_sheets[unique_sheet_name(gene, used_all)] = detail
            by_genus_sheets[genus][unique_sheet_name(gene, used_by_genus[genus])] = detail

    if n_fallback:
        print(f"[!] {n_fallback} gene(s) were scored with the final all-genes model "
              f"(EXPORT_OUT_OF_FOLD={EXPORT_OUT_OF_FOLD}); their scores are optimistic. "
              f"See the 'Scored by' column.")

    global_path = os.path.join(OUTPUT_DIR, "promoter_window_scoring_summary.xlsx")
    write_workbook(global_path, all_rows, all_sheets, drop_genus_col=False)
    n_hit = sum(1 for r in all_rows if r[score_col] > 0)
    print(f"[\u2713] Global workbook: {len(all_rows)} genes, {n_hit} with >=1 call "
          f"(relative score >= {PROMOTER_SCORE_THRESHOLD}) -> {global_path}")

    for genus, rows in sorted(by_genus_rows.items()):
        p = os.path.join(ORG_ROOT, safe_name(genus), f"{safe_name(genus)}_promoter_window_scoring_summary.xlsx")
        write_workbook(p, rows, by_genus_sheets[genus], drop_genus_col=True)
        print(f"    {genus}: {len(rows)} gene(s), "
              f"{sum(1 for r in rows if r[score_col] > 0)} with >=1 call -> {p}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def parse_args():
    """Command-line overrides for the module-level defaults (no interactive prompts)."""
    global FASTA_PATH, LABELS_PATH, OUTPUT_DIR, ORG_ROOT, RANDOM_SEED, N_FOLDS
    global PROMOTER_SCORE_THRESHOLD, EXPORT_OUT_OF_FOLD
    ap = argparse.ArgumentParser(description="Promoter-window scoring for upstream-only promoter FASTA (stage 05g).")
    ap.add_argument("--fasta", default=FASTA_PATH, help="upstream-only promoter FASTA")
    ap.add_argument("--labels", default=LABELS_PATH, help="CSV with seq_id, organism (gene optional)")
    ap.add_argument("--outdir", default=OUTPUT_DIR, help="folder for the global workbook")
    ap.add_argument("--org-root", default=ORG_ROOT, help="folder for per-genus workbooks")
    ap.add_argument("--seed", type=int, default=RANDOM_SEED)
    ap.add_argument("--folds", type=int, default=N_FOLDS)
    ap.add_argument("--threshold", type=float, default=PROMOTER_SCORE_THRESHOLD,
                    help="relative-score threshold used for window calls")
    ap.add_argument("--no-out-of-fold", action="store_true",
                    help="export scores from the final all-genes model instead of out-of-fold models")
    a = ap.parse_args()
    FASTA_PATH, LABELS_PATH, OUTPUT_DIR, ORG_ROOT = a.fasta, a.labels, a.outdir, a.org_root
    RANDOM_SEED, N_FOLDS, PROMOTER_SCORE_THRESHOLD = a.seed, a.folds, a.threshold
    EXPORT_OUT_OF_FOLD = not a.no_out_of_fold
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)


def main():
    parse_args()
    if not os.path.exists(FASTA_PATH):
        sys.exit(f"{FASTA_PATH} not found -- run 01e_new9_fetch_family_and_promoters.py first.")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    n_in_fasta = sum(1 for _ in SeqIO.parse(FASTA_PATH, "fasta"))
    records, dropped = load_promoters(FASTA_PATH)
    print(f"FASTA records: {n_in_fasta} | usable: {len(records)} | dropped: {len(dropped)}")
    for sid, why in dropped:
        print(f"   dropped {sid}: {why}")
    if len(records) < MIN_GENES_FOR_TRAINING:
        sys.exit(f"Only {len(records)} usable sequences -- need at least {MIN_GENES_FOR_TRAINING}.")

    lengths = Counter(len(s) for _, s in records)
    print(f"Record lengths (most common): {lengths.most_common(5)}")
    near = sum(c for l, c in lengths.items() if abs(l - EXPECTED_PROMOTER_LEN) <= 0.05 * EXPECTED_PROMOTER_LEN)
    if near / len(records) < 0.8:
        print(f"\n[!] WARNING: fewer than 80% of records are ~{EXPECTED_PROMOTER_LEN} bp. This script "
              f"treats the LAST base of every record as the base just before the annotated gene start. "
              f"If your records carry a downstream flank (e.g. 1201 bp = 1000 up + start + 200 down), "
              f"this script would mislabel them.\n")

    genus_map, gene_names = load_label_maps(LABELS_PATH)
    counts = Counter(genus_map.get(sid, "Unknown_genus") for sid, _ in records)
    print(f"Genera ({len(counts)}): {dict(counts)}")

    fold_models = run_cross_validation(records, genus_map)

    print("Training final model on all genes...")
    final_model = train_model(records)
    print("Done.\n")

    export_all(records, genus_map, gene_names, final_model, fold_models)


if __name__ == "__main__":
    main()
