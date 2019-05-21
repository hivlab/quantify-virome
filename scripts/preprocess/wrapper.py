
from snakemake.shell import shell

# bbduk adapter trimming parameters.
bbduk = snakemake.params.get("bbduk", "")

# Subsampling parameters.
frac = snakemake.params.get("frac", "1")
seed = snakemake.params.get("seed", "11")
print("Sampling {} of reads using seed {}.".format(frac, seed))

# Preprocessing command to run.
commands = [
            "bbmerge.sh in1={snakemake.input[0]} in2={snakemake.input[1]} outa={snakemake.output.adapters}",
            "bbmerge.sh in1={snakemake.input[0]} in2={snakemake.input[1]} out={snakemake.output.merged} outu={snakemake.output.unmerged} adapters={snakemake.output.adapters}",
            "cat {snakemake.output.merged} {snakemake.output.unmerged} > {snakemake.output.reads}",
            "bbduk.sh in={snakemake.output.reads} out={snakemake.output.trimmed} ref={snakemake.output.adapters} {bbduk}"
            ]

# Run preprocessing commands.
for cmd in commands:
  shell(cmd)

# If sample fraction is given, subsample reads using seed.
# Otherwise symlink trimmed reads to final output.
if frac and frac < 1:
  shell("reformat.sh in={snakemake.output.trimmed} out={snakemake.output.sampled} samplerate={frac} sampleseed={seed}")
else:
  shell("ln -sr {snakemake.output.trimmed} {snakemake.output.sampled}")
