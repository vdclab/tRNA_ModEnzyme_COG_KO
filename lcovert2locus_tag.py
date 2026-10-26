import csv, sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t", dtype=str, quoting=csv.QUOTE_NONE, keep_default_na=False)

def join_tags(x):
    return ";".join(dict.fromkeys(t for t in x if t))

mat = (df.groupby(["COG Assembly", "COG"])["COG2024 locus tag"]
         .agg(join_tags)
         .unstack(fill_value=""))

mat.to_csv(output_file, sep="\t")
