
# MegaBlast against reference genome to remove host sequences
rule megablast_refgenome:
    input:
      query = rules.repeatmasker_good.output.masked_filt
    output:
      out = temp("refgenomefilter/{sample}_megablast_{n}.tsv")
    params:
      db = config["ref_genome"],
      task = "megablast",
      perc_identity = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = "'6 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    wrapper:
      config["wrappers"]["blast"]

# Filter megablast records for the cutoff value
rule parse_megablast:
    input:
      blast_result = rules.megablast_refgenome.output.out,
      query = rules.repeatmasker_good.output.masked_filt
    output:
      mapped = temp("refgenomefilter/{sample}_refgenome_megablast_{n}_known-host.tsv"),
      unmapped = temp("refgenomefilter/{sample}_refgenome_megablast_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Collect stats
rule refgenome_unmapped_stats:
  input:
    expand("refgenomefilter/{{sample}}_refgenome_megablast_{n}_unmapped.fa", n = N)
  output:
    "stats/{sample}_refgenomeblast.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]
