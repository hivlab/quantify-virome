
import json
import re
import pandas as pd

def collect_stats(fastp, cdhit, rm, rm_good, refgenome_unmapped, refgenome_filtered, outfile): 
    
    # compile regex
    digs = re.compile("\\d+")
    
    # fastp
    with open(fastp, "r") as read_file:
        data = json.load(read_file)
        fastp = {"path": fastp, "summary": data["summary"]}

    # cdhit
    with open(cdhit, "r") as read_file:
        for line in read_file:
            if "total seq" in line:
                cdhit = {"path": cdhit, "sequences": int(re.findall(digs, line)[0])}

    # repatmasker
    repeatmasker = []
    for path in rm:
        with open(path, "r") as read_file:
            for line in read_file:
                if "sequences" in line:
                    entry = {"path": path, "sequences": int(re.findall(digs, line)[0])}
                    repeatmasker.append(entry)

    # repeatmasker_good
    repeatmasker_good = []
    for path in rm_good:
        entry = {"path": path, "sequences": len([1 for line in open(path) if line.startswith(">")])}
        repeatmasker_good.append(entry) 

    # refgenome_unmapped
    refgenome_unmapped = []
    for path in refgenome_unmapped:
        entry = {"path": path, "sequences": len([1 for line in open(path) if line.startswith(">")])}
        refgenome_unmapped.append(entry)

    # megablast
    megablast = []
    for path in refgenome_filtered:
        try:
            hits = pd.read_table(path, sep = "\\s+")
            entry = {"path": path, "sequences": len(set(hits["qseqid"]))}
        except:
            entry = {"path": path, "sequences": 0}
        megablast.append(entry)

    # write stats to json
    stats = {"fastp": fastp, "cdhit": cdhit, "tantan_good": repeatmasker, "repeatmasker_good": repeatmasker_good, "refgenome_unmapped": refgenome_unmapped, "megablast": megablast}
    with open(outfile, "w") as write_file:
        json.dump(stats, write_file)

if __name__ == '__main__':
    # Merge function arguments into dictionary.
    options = dict(snakemake.input)
    options.update(snakemake.output)
    # Unwrap arguments and run parse_blast
    collect_stats(**options)
