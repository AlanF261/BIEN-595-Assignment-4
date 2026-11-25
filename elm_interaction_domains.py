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
        
def write_tsv(path, header, mapping):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for key, values in sorted(mapping.items()):
            writer.writerow([key, ",".join(sorted(values))])

write_tsv("dom2elm.tsv", ["domain_id", "elm_ids"], dom2elm)
write_tsv("elm2dom.tsv", ["elm_id", "domain_ids"], elm2dom)