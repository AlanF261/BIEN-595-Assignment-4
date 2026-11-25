import csv
import subprocess
import os
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj

def makeTree(D, ids, outPath):
    dm = DistanceMatrix(D, ids=ids)
    tree = nj(dm)
    tree.write(outPath)

def parse_LMBD(path):
    ELM_identifiers = []
    interactionDomainID = []
    interactionDomainDescription = []
    interactionDomainName = []

    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for line in reader:
            # line = line.strip()
            # terms = line.split("\t")
            # if len(line) == 0:
            #     continue
            if line[0] == ('ELM identifier'):
                continue
            else:
                # terms = [term.strip('"') for term in terms]
                ELM_identifiers.append(line[0])
                interactionDomainID.append(line[1])
                interactionDomainDescription.append(line[2])
                interactionDomainName.append(line[3])
    return ELM_identifiers, interactionDomainID, interactionDomainDescription, interactionDomainName

def buildHHM(interactionDomainID):
    interactionDomainID = set(interactionDomainID)
    interactionDomainID = {ID for ID in interactionDomainID if ID.startswith("PF")}
    # a3mPaths = []
    hhrPaths = []
    IDs = []
    hhmPaths = []
    awkLine = r'$2=="AC" && index($3, id".")==1 {f=1}  f && !/^#/ && NF==2 {print ">"$1; print $2}  /^\/\/$/ {f=0}'
    for ID in interactionDomainID:
        # awk f'$3=="ID" && $3=="{ID}"{f=1} f; /^\/\/$/{f=0}' Pfam-A.seed > f"{ID}.sto"
        # # stoPaths.append(f"{ID}".sto)
        # awk '!/^#/ && NF==2 {print ">"$1; print $2}' f"{ID}".sto > f"{ID}.a3m"
        with open(f"a3m/{ID}.a3m", "w") as out:
            subprocess.run(["awk", "-v", f"id={ID}", awkLine, "data/Pfam-A.seed"], stdout=out, check=True)
        if os.path.getsize(f"a3m/{ID}.a3m") == 0:
            print("not found:", ID)
            continue
        # a3mPaths.append(f"a3m/{ID}.a3m")
        IDs.append(ID)
    for ID in IDs:
        # newPath = a3mPath[:-4] + ".hhm"
        hhmPaths.append(f"hhm/{ID}.hhm")
        subprocess.run(["hhmake", "-i", f"a3m/{ID}.a3m", "-o", f"hhm/{ID}.hhm", "-M", "50", "-name", ID], check=True)
    for ID1 in IDs:
        for ID2 in IDs:
            # if hhmPath != hhmPath2:
            targetPath = f"hhm/{ID2}.hhm"
            hhrPath = f"hhr/{ID1}-{ID2}.hhr"
            subprocess.run(["hhalign", "-i", f"hhm/{ID1}.hhm", "-t", targetPath, "-o", hhrPath, "-v", "0"])
            hhrPaths.append(hhrPath)
    return hhrPaths

def makeTable(hhrPaths):

    scores = dict()
    # evals = dict()
    # probabs = dict()
    for hhrPath in hhrPaths:
        ID1, ID2 = hhrPath[4:-4].split("-")
        with open(hhrPath, "r") as f:
            for line in f:
                if line.startswith("Probab"):
                    line = line.strip()
                    line = line.split()
                    probab = line[0].split("=")
                    evalue = line[1].split("=")
                    score = line[2].split("=")
                    scores[(ID1, ID2)] = (float(score[1]), float(probab[1]), float(evalue[1]))
                    # probabs[(ID1, ID2)] = float(probab[1])
                    # evals[(ID1, ID2)] = float(evalue[1])
    with open("hhalign_scores.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["query", "target", "prob", "evalue", "score"])
        for (ID1, ID2) in scores:
            score, probab, evalue = scores[(ID1, ID2)]
            writer.writerow([ID1, ID2, probab, evalue, score])

def buildMatrix(csvPath):
    rows = []
    with open(csvPath, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            rows.append((row[0], row[1], float(row[4])))
    ids = sorted({r[0] for r in rows} | {r[1] for r in rows})
    index = {ID: i for i, ID in enumerate(ids)}
    S = np.zeros((len(ids), len(ids)))
    for query, target, score in rows:
        S[index[query], index[target]] = score
    return S, ids

def toDistance(S):
    S = (S + S.T) / 2
    selfScores = np.diag(S)
    D = 1 - S / np.sqrt(np.outer(selfScores, selfScores))
    D = np.clip(D, 0, 1)
    np.fill_diagonal(D, 0)
    return D

def main():
    # ELM_identifiers, interactionDomainID, interactionDomainDescription, interactionDomainName = parse_LMBD("elm_interaction_domains.tsv")
    # # interactionDomainID = {"PF00675", "PF00082"}
    # hhrPaths = buildHHM(interactionDomainID)
    # makeTable(hhrPaths)
    S, ids = buildMatrix("hhalign_scores.csv")
    D = toDistance(S)
    makeTree(D, ids, "lmbd_tree.nwk")

if __name__ == "__main__":
    main()