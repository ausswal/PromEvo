#!/usr/bin/env python3
"""
STAGE 05g — Promoter-window prediction for UPSTREAM-ONLY promoter FASTA
(ONE script that replaces 05b + 05d + 05f).

WHEN TO USE THIS SCRIPT
-----------------------
Use it when every record in outputs/promoter_sequences.fasta is the N bp
(e.g. 1000 bp) that lie UPSTREAM of the annotated TSS, with NO downstream
flank. The last base of each record is therefore the base immediately
before the TSS.

05b/05d/05f assume the opposite: a 1,201 bp record (1000 up + TSS + 200
down), where the TSS sits at index 1000. On a 1000 bp record they place the
"TSS" at index 799 (sequence length - 201), mislabel every training example
and silently throw away the last 200 bp. Just editing UPSTREAM_BP /
DOWNSTREAM_BP does not rescue them either: with DOWNSTREAM_BP = 0 the
TSS-centred training window needs 50 bp that do not exist, so zero positive
examples are built and training crashes. If your records DO contain a
downstream flank (1201 bp), keep using 05b/05d/05f.

WHAT THIS SCRIPT DOES (same method as 05b/05d/05f, adapted to the data)
-----------------------------------------------------------------------
1. Positive example per gene = the WINDOW_SIZE bp immediately upstream of
   the annotated TSS (the last WINDOW_SIZE bases of the record). Negatives
   = random windows far upstream, "hard" windows 50-300 bp upstream of the
   positive, and a shuffled copy of the positive.
2. Two scorers, blended 50/50 after min-max scaling inside each gene:
   a Position Weight Matrix (log-odds) and an MLP on 3/4/5-mer frequencies.
3. Cluster-aware 5-fold cross-validation (homologous promoters never split
   between train and test) with a per-genus breakdown. Because the records
   have no flank, the test is: scan the WHOLE held-out promoter, and check
   whether the top-scoring window is the true proximal window. It is
   compared with what random window picking would give (chance).
4. Final model on all genes -> Excel export: one global workbook and one
   workbook per genus. EVERY gene gets a row (best window + Confidence).

HOW TO READ THE OUTPUT -- IMPORTANT LIMITS
------------------------------------------
* The "TSS" is the ANNOTATED TSS from your fetch stage, not an
  experimentally mapped one. The model learns what the sequence just
  upstream of that annotated position looks like.
* Coordinates are relative to that annotated TSS: -1 = last base upstream.
* The 0.80 "Score" is min-max scaled WITHIN each gene, so it is a relative
  rank, not a probability. The "MLP probability" column is the absolute
  classifier output for the same window -- look at both.
* The final model is trained on all genes and then scores those same genes,
  so exported scores are optimistic; the cross-validation numbers printed
  first are the honest accuracy estimate.
* The best window tends to be the one nearest the TSS, because that is where
  every positive example came from. A real promoter signal needs the CV hit
  rate to be clearly above the printed chance rate.

Requirements:
    pip install biopython scikit-learn numpy pandas openpyxl
"""

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
OUTPUT_DIR = "outputs/tss_prediction"                # global workbook goes here
ORG_ROOT = "outputs/promoters_by_organism"           # <Genus>/<Genus>_promoter_prediction_summary.xlsx

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
    more than MAX_AMBIGUOUS_FRACTION ambiguous bases. (05b dropped a
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
    """Centre of the window that ENDS at the annotated TSS (= last base of the record)."""
    return len(seq) - HALF


def shuffle_decoy(win):
    """Mononucleotide-preserving shuffle (NOT dinucleotide-preserving)."""
    chars = list(win)
    random.shuffle(chars)
    return "".join(chars)


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

        negatives.append(shuffle_decoy(pos_win))    # shuffled copy of the positive
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


def scan_sequence(seq, pwm, clf, scaler):
    """Scores every window across the WHOLE promoter. Returns
    (centres, combined_relative_score, mlp_probability)."""
    centers = list(range(HALF, len(seq) - HALF + 1, SCAN_STEP))
    if centers[-1] != len(seq) - HALF:              # always include the window ending at the TSS
        centers.append(len(seq) - HALF)
    idx = _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    cols = np.arange(WINDOW_SIZE)

    pwm_scores, wins = [], []
    for c in centers:
        s = c - HALF
        sl = idx[s:s + WINDOW_SIZE]
        ok = sl >= 0
        pwm_scores.append(float(pwm[sl[ok], cols[ok]].sum()))
        wins.append(seq[s:s + WINDOW_SIZE])

    X = np.array([kmer_features(w) for w in wins])
    mlp = clf.predict_proba(scaler.transform(X))[:, 1]

    def norm(a):
        a = np.asarray(a, dtype=float)
        return np.zeros_like(a) if a.max() == a.min() else (a - a.min()) / (a.max() - a.min())

    combined = COMBINE_WEIGHT_PWM * norm(pwm_scores) + (1 - COMBINE_WEIGHT_PWM) * norm(mlp)
    return np.array(centers), combined, mlp


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
def run_cross_validation(records, genus_map):
    folds, n_clusters = make_folds(records, N_FOLDS)
    print(f"\n=== {len(folds)}-fold cluster-aware cross-validation "
          f"({n_clusters} similarity clusters among {len(records)} genes) ===")
    errs, base_errs, hits, chance, per_genus = [], [], [], [], defaultdict(list)

    for fi, test in enumerate(folds):
        train = [r for j, f in enumerate(folds) if j != fi for r in f]
        if not test or len(train) < MIN_GENES_FOR_TRAINING:
            continue
        pwm, clf, scaler = train_model(train)
        fold_errs = []
        for sid, seq in test:
            centers, combined, _ = scan_sequence(seq, pwm, clf, scaler)
            truth = true_window_center(seq)
            dist = np.abs(centers - truth)
            err = int(dist[int(np.argmax(combined))])
            errs.append(err)
            fold_errs.append(err)
            base_errs.append(float(dist.mean()))           # expected error of a random pick
            hits.append(err <= WINDOW_SIZE)
            chance.append(float((dist <= WINDOW_SIZE).mean()))
            per_genus[genus_map.get(sid, "Unknown_genus")].append((err, err <= WINDOW_SIZE))
        print(f"Fold {fi + 1}/{len(folds)}: {len(test)} genes -> median abs error "
              f"{np.median(fold_errs):.0f} bp")

    if errs:
        print(f"\nTop-scoring window vs the true proximal window (n={len(errs)} held-out genes):")
        print(f"  median abs error     : {np.median(errs):.0f} bp   "
              f"(random-pick baseline ~{np.median(base_errs):.0f} bp)")
        print(f"  within {WINDOW_SIZE} bp of truth  : {np.mean(hits) * 100:.1f}%   "
              f"(chance ~{np.mean(chance) * 100:.1f}%)")
        print("  -> only trust the exported scores if the hit rate is clearly above chance.")
    print("\n--- per-genus (same global model, NOT trained per genus) ---")
    for g, v in sorted(per_genus.items(), key=lambda kv: -len(kv[1])):
        if len(v) < MIN_GENES_FOR_GENUS_STAT:
            print(f"  {g}: only {len(v)} held-out gene(s) -- too few for a stable estimate.")
        else:
            print(f"  {g}: n={len(v)}, median abs error {np.median([e for e, _ in v]):.0f} bp, "
                  f"within {WINDOW_SIZE} bp: {np.mean([h for _, h in v]) * 100:.0f}%")
    print("=== end cross-validation ===\n")


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


def describe_gene(seq_id, seq, gene, genus, pwm, clf, scaler, score_col):
    L = len(seq)                                     # annotated TSS = boundary after the last base
    centers, combined, mlp = scan_sequence(seq, pwm, clf, scaler)
    calls = call_promoters(centers, combined, PROMOTER_SCORE_THRESHOLD)

    bi = int(np.argmax(combined))                    # ALWAYS report the best window
    bc, bs = int(centers[bi]), float(combined[bi])
    row = {
        "Gene": gene, "Seq_ID": seq_id, "Genus": genus,
        score_col: len(calls),
        "Promoter with Highest Prediction Score (Start-End)": f"{bc - HALF + 1}-{bc + HALF}",
        "Score": round(bs, 2),
        "MLP probability": round(float(mlp[bi]), 3),
        "Location in 1kb Sequence (relative to TSS)": f"{(bc - HALF) - L}\u2192{(bc + HALF - 1) - L}",
        "Confidence": f">= {PROMOTER_SCORE_THRESHOLD:.2f}" if bs >= PROMOTER_SCORE_THRESHOLD
                      else f"< {PROMOTER_SCORE_THRESHOLD:.2f}",
    }
    detail = None
    if calls:
        rows = []
        for c in sorted(calls, key=lambda c: c["peak_score"], reverse=True):
            rows.append({"Start-End (local)": f"{c['start'] + 1}-{c['end'] + 1}",
                         "Score": round(float(c["peak_score"]), 2),
                         "Location (relative to TSS)": f"{c['start'] - L}\u2192{c['end'] - L}"})
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


def export_all(records, genus_map, gene_names, pwm, clf, scaler):
    score_col = f"#Promoters (Score>={PROMOTER_SCORE_THRESHOLD:.2f})"
    all_rows, all_sheets, used_all = [], {}, set()
    by_genus_rows, by_genus_sheets, used_by_genus = defaultdict(list), defaultdict(dict), defaultdict(set)

    for sid, seq in records:
        genus = genus_map.get(sid, "Unknown_genus")
        gene = gene_names.get(sid, sid)
        row, detail = describe_gene(sid, seq, gene, genus, pwm, clf, scaler, score_col)
        all_rows.append(row)
        by_genus_rows[genus].append(row)
        if detail is not None:
            all_sheets[unique_sheet_name(gene, used_all)] = detail
            by_genus_sheets[genus][unique_sheet_name(gene, used_by_genus[genus])] = detail

    global_path = os.path.join(OUTPUT_DIR, "promoter_prediction_summary.xlsx")
    write_workbook(global_path, all_rows, all_sheets, drop_genus_col=False)
    n_hit = sum(1 for r in all_rows if r[score_col] > 0)
    print(f"[\u2713] Global workbook: {len(all_rows)} genes, {n_hit} with >=1 call "
          f"(score >= {PROMOTER_SCORE_THRESHOLD}) -> {global_path}")

    for genus, rows in sorted(by_genus_rows.items()):
        p = os.path.join(ORG_ROOT, safe_name(genus), f"{safe_name(genus)}_promoter_prediction_summary.xlsx")
        write_workbook(p, rows, by_genus_sheets[genus], drop_genus_col=True)
        print(f"    {genus}: {len(rows)} gene(s), "
              f"{sum(1 for r in rows if r[score_col] > 0)} with >=1 call -> {p}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    if not os.path.exists(FASTA_PATH):
        sys.exit(f"{FASTA_PATH} not found -- run 05_fetch_promoter_sequences_for_all.py first.")
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
              f"treats the LAST base of every record as the base just before the TSS. If your records "
              f"carry a downstream flank (e.g. 1201 bp = 1000 up + TSS + 200 down), use "
              f"05b/05d/05f instead -- this script would mislabel them.\n")

    genus_map, gene_names = load_label_maps(LABELS_PATH)
    counts = Counter(genus_map.get(sid, "Unknown_genus") for sid, _ in records)
    print(f"Genera ({len(counts)}): {dict(counts)}")

    run_cross_validation(records, genus_map)

    print("Training final model on all genes...")
    pwm, clf, scaler = train_model(records)
    print("Done.\n")

    export_all(records, genus_map, gene_names, pwm, clf, scaler)


if __name__ == "__main__":
    main()
