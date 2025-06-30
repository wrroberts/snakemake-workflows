"""
The blastx output files must have a customized tabular output:
0-qseqid 1-qlen 2-seqid 3-slen 4-frame 5-pident
6-nident 7-length 8-mismatch 9-gapopen 10-qstart 11-qend
12-start 13-send 14-evalue 15-bitscore

Detects chimeras using:
Step 1: Filter HSPs by length and similarity
Step 2: Check for self chimeras in query-hit blocks
Step 3: Check for multi-gene chimeras in full query blocks

Outputs:
- `.cut` file with regions to cut
- `.info` file with HSP blocks for manual validation
"""

import sys
import os

PIDENT_CUTOFF = 30  # Minimum percent identity to keep an HSP
LENGTH_CUTOFF = 100  # Minimum query coverage length

def qcov(hsp):
    return abs(hsp[11] - hsp[10]) + 1

def separated(hsp1, hsp2):
    length1 = qcov(hsp1)
    length2 = qcov(hsp2)
    start = min(hsp1[10], hsp1[11], hsp2[10], hsp2[11])
    end = max(hsp1[10], hsp1[11], hsp2[10], hsp2[11])
    overlap = length1 + length2 - (end - start) + 1
    return overlap < min(60, 0.2 * min(length1, length2))

def expand_range(hsp1, hsp2):
    if hsp1 == []:
        return hsp2
    if hsp2 == []:
        return hsp1
    start1, end1, start2, end2 = hsp1[10], hsp1[11], hsp2[10], hsp2[11]
    if start1 < end1 and start2 < end2:
        start, end = min(start1, start2), max(end1, end2)
    elif start1 > end1 and start2 > end2:
        start, end = max(start1, start2), min(end1, end2)
    else:
        return hsp1  # Don't merge opposite directions
    hsp1[10], hsp1[11] = start, end
    return hsp1

def check_block(block, multigene):
    if len(block) == 1:
        return True
    pos, neg = [], []
    for hsp in block:
        if str(hsp[4]).startswith("-"):
            neg = expand_range(neg, hsp)
        else:
            pos = expand_range(pos, hsp)
    if (not pos and neg) or (not neg and pos):
        return True
    elif separated(pos, neg):
        if multigene:
            for h in [pos, neg]:
                start, end = min(h[10], h[11]), max(h[10], h[11])
                outfile1.write(f"{h[0]} {int(start)} {int(end)} trans-multi\n")
        else:
            outhsp = pos if qcov(pos) > qcov(neg) else neg
            start, end = min(outhsp[10], outhsp[11]), max(outhsp[10], outhsp[11])
            outfile1.write(f"{outhsp[0]} {int(start)} {int(end)} trans-self\n")
        for i in pos:
            outfile2.write(f"{i}\t")
        outfile2.write("\n")
        for i in neg:
            outfile2.write(f"{i}\t")
        outfile2.write("\n")
        return False
    else:
        return True

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python detect_chimera_from_blastx.py blastx_output output_dir")
        sys.exit()

    blastx_output = sys.argv[1]
    DIR = sys.argv[2]

    if DIR == ".":
        DIR = os.getcwd()
    if not os.path.isabs(DIR):
        DIR = os.path.abspath(DIR)
    if not DIR.endswith("/"):
        DIR += "/"

    blastx_base_name = os.path.splitext(os.path.basename(blastx_output))[0]
    infile = open(blastx_output, "r", encoding="utf-8")
    outfile1 = open(f"{DIR}{blastx_base_name}.cut", "w")
    outfile2 = open(f"{DIR}{blastx_base_name}.info", "w")

    last_query = ""
    for line in infile:
        if len(line) < 3:
            continue
        hsp = line.strip().split("\t")
        for i in [5, 10, 11]:
            hsp[i] = float(hsp[i])
        if hsp[5] < PIDENT_CUTOFF or qcov(hsp) < LENGTH_CUTOFF:
            continue
        query, hit = hsp[0], hsp[2]

        if last_query == "":
            hit_block = [hsp]
            query_block = [hsp]
            good_seq = True
        elif query == last_query:
            query_block.append(hsp)
            if good_seq:
                if hit == last_hit:
                    hit_block.append(hsp)
                else:
                    good_seq = check_block(hit_block, False)
                    hit_block = [hsp]
        else:
            if good_seq:
                good_seq = check_block(hit_block, False)
            if good_seq:
                good_seq = check_block(query_block, True)
            query_block, hit_block = [hsp], [hsp]
            good_seq = True

        last_query, last_hit = query, hit

    if good_seq:
        good_seq = check_block(hit_block, False)
    if good_seq:
        good_seq = check_block(query_block, True)

    infile.close()
    outfile1.close()
    outfile2.close()
