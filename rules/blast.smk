
# Prepare taxonomy annotation tables.
rule prepare_taxonomy_data:
  input: config["names"], config["nodes"], config["division"]
  output:
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
  conda:
    "../envs/tidyverse.yaml"
  script:
    "../scripts/prepare_taxonomy_data.R"

# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
      query = rules.parse_megablast.output.unmapped
    output:
      out = "blast/{sample}_blastn_virus_{n,\d+}.tsv"
    params:
      db = config["virus_nt"],
      task = "blastn",
      evalue = config["blastn_virus"]["evalue"],
      db_soft_mask = config["blastn_virus"]["db_soft_mask"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
      query = rules.parse_megablast.output.unmapped,
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = temp("blast/{sample}_blastn_virus_{n,\d+}_known-viral.tsv"),
      unmapped = temp("blast/{sample}_blastn_virus_{n,\d+}_unmapped.fa")
    params:
      e_cutoff = 1e-5,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped
    output:
      out = temp("blast/{sample}_blastx_virus_{n}.tsv")
    params:
      db = config["virus_nr"],
      word_size = 6,
      evalue = config["blastx_virus"]["evalue"],
      db_soft_mask = config["blastx_virus"]["db_soft_mask"],
      max_hsps = config["blastx_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastx_virus:
    input:
      query = rules.blastx_virus.input.query,
      blast_result = rules.blastx_virus.output.out
    output:
      mapped = temp("blast/{sample}_blastx_virus_{n}_known-viral.tsv"),
      unmapped = temp("blast/{sample}_blastx_virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Filter sequences by division id.
# Saves hits with division id
rule classify_phages:
  input:
    [rules.parse_blastn_virus.output.mapped, rules.parse_blastx_virus.output.mapped] if config["run_blastx"] else rules.parse_blastn_virus.output.mapped,
    nodes = "taxonomy/nodes.csv"
  output:
    division = "results/{sample}_phages_{n}.csv",
    other = temp("blast/{sample}_candidate_viruses_{n}.csv")
  params:
    taxdb = config["vhunter"],
    division_id = 3
  wrapper:
    config["wrappers"]["blast_taxonomy"]

# Filter unmasked candidate virus reads.
rule unmasked_other:
    input:
      rules.classify_phages.output.other,
      rules.refgenome_unmapped.output
    output:
      temp("blast/{sample}_candidate_viruses_{n}_unmasked.fa")
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmasked_viral.py"

# Map reads against bacterial genomes.
rule bwa_mem_refbac:
    input:
      reads = [rules.unmasked_other.output]
    output:
      temp("blast/{sample}_bac_mapped_{n}.bam")
    params:
      index = config["ref_bacteria"],
      extra = "-k 15",
      sort = "none"
    log:
      "logs/{sample}_bwa_map_refbac_{n}.log"
    threads: 2
    wrapper:
      "0.32.0/bio/bwa/mem"

# Extract unmapped reads and convert bam file to fastq file.
rule refbac_unmapped_fastq:
    input:
      rules.bwa_mem_refbac.output
    output:
      temp("blast/{sample}_refbac_unmapped_{n}.fq")
    params:
      "-n -f 4"
    threads: 2
    wrapper:
      "0.32.0/bio/samtools/bam2fq/interleaved"

# Convert fastq file to fasta file.
rule refbac_unmapped:
    input:
      rules.refbac_unmapped_fastq.output
    output:
      temp("blast/{sample}_refbac_unmapped_{n}.fa")
    shell:
      "cat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}"

# Calculate bam file stats
rule refbac_bam_stats:
    input:
      rules.bwa_mem_refbac.output
    output:
      "stats/{sample}_refbac_stats_{n}.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"

# Subset repeatmasker masked reads using unmapped reads.
rule refbac_unmapped_masked:
    input:
      rules.refbac_unmapped.output,
      rules.repeatmasker_good.output.masked_filt
    output:
      temp("blast/{sample}_bac_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmapped_masked_ids.py"

# Megablast against nt database.
rule megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output
    output:
      out = temp("blast/{sample}_megablast_nt_{n,\d+}.tsv")
    params:
      db = config["nt"],
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      max_hsps = config["megablast_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter megablast hits for the cutoff value.
rule parse_megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output,
      blast_result = rules.megablast_nt.output.out
    output:
      mapped = temp("blast/{sample}_megablast_nt_{n,\d+}_mapped.tsv"),
      unmapped = temp("blast/{sample}_megablast_nt_{n,\d+}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Collect stats
rule refbac_unmapped_stats:
  input:
    expand(["blast/{{sample}}_bac_unmapped_{n}_masked.fa", "blast/{{sample}}_megablast_nt_{n}_unmapped.fa"], n = N)
  output:
    "stats/{sample}_refbac.tsv"
  params:
    extra = "-T"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/dfff20d4f55ed7b9e52afa34f57a4556e295680f/bio/seqkit/stats"

# Blastn against nt database.
rule blastn_nt:
    input:
      query = rules.parse_megablast_nt.output.unmapped
    output:
      out = temp("blast/{sample}_blastn_nt_{n,\d+}.tsv")
    params:
      db = config["nt"],
      task = "blastn",
      evalue = config["blastn_nt"]["evalue"],
      max_hsps = config["blastn_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn records for the cutoff value.
rule parse_blastn_nt:
    input:
      query = rules.blastn_nt.input.query,
      blast_result = rules.blastn_nt.output.out
    output:
      mapped = temp("blast/{sample}_blastn_nt_{n,\d+}_mapped.tsv"),
      unmapped = temp("blast/{sample}_blastn_nt_{n,\d+}_unmapped.fa") if config["run_blastx"] else temp("results/{sample}_unassigned_{n,\d+}.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Blastx unmapped sequences against nr database.
rule blastx_nr:
    input:
      query = rules.parse_blastn_nt.output.unmapped
    output:
      out = temp("blast/{sample}_blastx_nr_{n,\d+}.tsv")
    params:
      db = config["nr"],
      word_size = 6,
      evalue = config["blastx_nr"]["evalue"],
      max_hsps = config["blastx_nr"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastx records for the cutoff value.
rule parse_blastx_nr:
    input:
      query = rules.blastx_nr.input.query,
      blast_result = rules.blastx_nr.output.out
    output:
      mapped = temp("blast/{sample}_blastx_nr_{n,\d+}_mapped.tsv"),
      unmapped = "results/{sample}_unassigned_{n,\d+}.fa" if config["run_blastx"] else "{sample}_None_{n}"
    params:
      e_cutoff = 1e-3,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Filter sequences by division id.
# Saves hits with division id
rule classify_phages_viruses:
  input:
    [rules.parse_megablast_nt.output.mapped, rules.parse_blastn_nt.output.mapped, rules.parse_blastx_nr.output.mapped] if config["run_blastx"] else [rules.parse_megablast_nt.output.mapped, rules.parse_blastn_nt.output.mapped],
    nodes = "taxonomy/nodes.csv"
  output:
    division = "results/{sample}_phages_viruses_{n}.csv",
    other = "results/{sample}_non_viral_{n}.csv"
  params:
    taxdb = config["vhunter"],
    division_id = [3, 9] # pool phages and viruses
  wrapper:
    config["wrappers"]["blast_taxonomy"]

# Upload results to Zenodo.
if config["zenodo"]["deposition_id"]:
  rule upload:
    input:
      expand("results/{{sample}}_{{result}}_{n}.{{ext}}", n = N)
    output:
      ZEN.remote(expand("{deposition_id}/files/results/{{sample, [^_]+}}_{{result}}.{{ext}}.tar.gz", deposition_id = config["zenodo"]["deposition_id"]))
    shell:
      "tar zcvf {output} {input}"

# Collect stats
rule blast_stats:
  input:
    expand(["blast/{{sample}}_blastn_virus_{n}_unmapped.fa",
    "blast/{{sample}}_blastx_virus_{n}_unmapped.fa",
    "blast/{{sample}}_candidate_viruses_{n}_unmasked.fa",
    "blast/{{sample}}_megablast_nt_{n}_unmapped.fa",
    "blast/{{sample}}_blastn_nt_{n}_unmapped.fa"], n = N) if config["run_blastx"] else expand("results/{{sample}}_unassigned_{n}.fa", n = N)
  output:
    "stats/{sample}_blast.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]
