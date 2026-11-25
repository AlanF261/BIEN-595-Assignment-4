import csv
import collections
dom2elm = collections.defaultdict(set)
elm2dom = collections.defaultdict(set)
with open("elm_interaction_domains.tsv", newline="", encoding="utf-8") as f:
    rows = csv.reader(f, delimiter="\t")
    next(rows)                            # header
    for elm, dom, *_ in rows:
        dom2elm[dom].add(elm)
        elm2dom[elm].add(dom)
        
