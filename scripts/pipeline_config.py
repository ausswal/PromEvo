#!/usr/bin/env python3
"""
Shared pipeline-wide constants and consistency checks.

Every stage that assumes a fixed promoter window (06a, 06b, 06c, 07i, 05g,
08, 08c, 010, ...) should import PROMOTER_LEN from here, and the fetch stage
(01e) should use the SAME value for its upstream length instead of letting
the user type an arbitrary one:

    from pipeline_config import PROMOTER_LEN, check_promoter_lengths
    # 01e:        PROMOTER_UPSTREAM_BP = PROMOTER_LEN       (no user prompt)
    # other stages, right after parsing the promoter FASTA:
    check_promoter_lengths([len(r.seq) for r in records], stage="06a")

Promoters SHORTER than PROMOTER_LEN are legitimate (clipped at a contig end)
and never an error. Promoters LONGER than PROMOTER_LEN mean the fetch stage
ran with a different length, so downstream TSS-relative coordinates would
not be comparable -- that stops the stage unless it trims itself (trims=True)
or you set the environment variable PIPELINE_ALLOW_LEN_MISMATCH=1.
"""
import os
import sys

PROMOTER_LEN = 1000   # bp upstream of the TSS; the ONE place this is defined


def check_promoter_lengths(lengths, stage, trims=False):
    lengths = [int(x) for x in lengths]
    if not lengths:
        return True
    longest = max(lengths)
    n_over = sum(1 for x in lengths if x > PROMOTER_LEN)
    n_short = sum(1 for x in lengths if x < PROMOTER_LEN)

    if longest < PROMOTER_LEN:
        print(f"[!] [{stage}] longest promoter is {longest} bp but PROMOTER_LEN="
              f"{PROMOTER_LEN} -- was the fetch stage (01e) run with a different "
              f"upstream length? TSS-relative coordinates assume {PROMOTER_LEN} bp.")
    if n_over and not trims:
        msg = (f"[{stage}] {n_over}/{len(lengths)} promoter(s) are longer than "
               f"PROMOTER_LEN={PROMOTER_LEN} bp (longest = {longest} bp). Re-run the "
               f"fetch stage (01e) with {PROMOTER_LEN} bp, or set "
               f"PIPELINE_ALLOW_LEN_MISMATCH=1 to continue anyway.")
        if os.environ.get("PIPELINE_ALLOW_LEN_MISMATCH") == "1":
            print("[!] " + msg + " (continuing: override set)")
            return False
        sys.exit("[!] " + msg)
    if n_short:
        print(f"[i] [{stage}] {n_short}/{len(lengths)} promoter(s) are shorter than "
              f"{PROMOTER_LEN} bp (clipped at a contig end) -- expected, kept as-is.")
    return True
