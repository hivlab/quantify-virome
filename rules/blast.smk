
## Extract unmapped reads [12a]
rule unmapped_reads:
    input: rules.bwa_mem.output
    output:
      bam = "{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.bam",
      fq = "{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.fq",
      fa = "{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.fa"
    conda:
      "envs/bwa-sam-bed.yml"
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids [12b]
rule unmapped_masked:
    input: rules.unmapped_reads.output.fa, rules.repeatmasker_good.output.masked
    output:
      "{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa"
    conda:
      "envs/biopython.yml"
    script:
      "scripts/unmapped_masked_ids.py"

## MegaBlast against reference genome to remove more host sequences [13]
rule megablast_ref_genome:
    input:
      db = config["ref_genome"],
      query = rules.unmapped_masked.output
    output:
      "{sample}/13_megablast/megablast.{n}.xml"
    params:
      perc_ident = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      num_desc = config["megablast_ref_genome"]["num_descriptions"],
      num_align = config["megablast_ref_genome"]["num_alignments"]
    threads: 8
    conda:
      "envs/biopython.yml"
    script:
      "scripts/megablast_ref_genome.py"

## Filter megablast records for the cutoff value [14]
rule parse_megablast:
    input:
      blastxml = rules.megablast_ref_genome.output,
      query = rules.unmapped_masked.output
    output:
      known = "{sample}/14_megablast_parsed/RefGenome_megablast.{n}.non-viral.out",
      unmapped = "{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa"
    params:
      e_cutoff = 1e-10
    conda:
      "envs/biopython.yml"
    script:
      "scripts/parse_blast_xml.py"

## Blast against virus database [15]
rule blastn_virus_nt:
    input:
      db = config["virus_nt"],
      query = rules.parse_megablast.output.unmapped
    output:
      out = "{sample}/15_blast_virusnt/blast_virusnt.{n}.xml"
    params:
      task = "blastn",
      show_gis = True,
      evalue = 1e-4,
      db_soft_mask = 100,
      num_threads = 8
    conda:
      "envs/biopython.yml"
    script:
      "scripts/blastn_virus_db.py"

## Filter blastn records for the cutoff value [16]
rule parse_virusntblast:
    input:
      blastxml = rules.blastn_virus_nt.output.out,
      query = rules.parse_megablast.output.unmapped
    output:
      known = "{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.known-viral.out",
      unmapped = "{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.unmapped.fa"
    params:
      e_cutoff = 1e-5
    conda:
      "envs/biopython.yml"
    script:
      "scripts/parse_blast_xml.py"

# Download taxonomy names [17a]
rule download_taxonomy:
    output:
      names = "names.csv",
      nodes = "nodes.csv"
    params:
      datadir = config["datadir"]
    conda:
      "envs/tidyverse.yml"
    script:
      "scripts/download_taxonomy_names.R"

def get_knownviral(wildcards):
  path = expand("{sample}/16_blastntvirus_parsed/blastnt_virus.*.known-viral.out", sample = wildcards.sample)
  return glob.glob(*path)

# Add taxonomy to virus nt blast [17b]
# n_ids is a global wildcard determined after split_fasta rule
rule virus_nt_taxonomy:
    input:
      known = get_knownviral,
      vhunter = config["vhunter"],
      names = rules.download_taxonomy.output.names,
      nodes = rules.download_taxonomy.output.nodes
    output:
      "{sample}/17_virus_nt_taxonomy/known_taxa.csv"
    conda:
      "envs/tidyverse.yml"
    script:
      "scripts/munge_taxonomy.R"

# Taxonomy report to virus nt blast [17c]
rule virus_nt_taxonomy_report:
    input:
      rules.virus_nt_taxonomy.output,
      names = rules.download_taxonomy.output.names
    output:
      "{sample}/17_virus_nt_taxonomy/taxonomy_report.html"
    params:
      lambda wildcards: wildcards.sample
    conda:
      "envs/tidyverse.yml"
    script:
      "scripts/taxonomy_report.Rmd"
