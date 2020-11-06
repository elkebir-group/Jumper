import pandas as pd
import sys
import re

def parse_attributes(attr: str):
    return dict(pair.split(' ') for pair in attr.split("; ") if len(pair) > 0)

pattern = re.compile(r'$.*#(\d+)^')

def name_increment(name: str):
    m = pattern.match(name)
    if m is not None:
        return name.split('#') + str(int(m.group(1)) + 1)
    else:
        return name + "#2"

def get_ORF_region(annotation_file: str):
    annot = pd.read_csv(annotation_file,
                        sep='\t',
                        header=None,
                        names=["id", "source", "feature", "start", "end",
                               "score", "strand", "phase", "attributes"],
                        comment='#')
    regions = []
    names = set()
    for idx, row in annot[annot["feature"] == "CDS"].iterrows():
        infos = parse_attributes(row["attributes"])
        if "gene_id" in infos:
            name = infos["gene_id"]
        elif "transcript_id" in infos:
            name = infos["transcript_id"]
        else:
            print("No label: ", row)
            exit(1)

        if name in names:
            name = name_increment(name)
        regions.append({"id": name, "start": row["start"]-1, "end": row["end"]-1})
        names.update(name)

    return pd.DataFrame(regions)

if __name__ == "__main__":
    print(get_ORF_region(sys.argv[1]))