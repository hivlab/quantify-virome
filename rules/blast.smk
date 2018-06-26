
## Extract unmapped reads [12a]
rule unmapped_reads:
    input:
      os.path.join(config["outdir"], dynamic("{sample}/11_bwa_mem/mapped.{n}.bam"))
    output:
      bam = os.path.join(config["outdir"], dynamic("{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.bam")),
      fq = os.path.join(config["outdir"], dynamic("{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.fq")),
      fa = os.path.join(config["outdir"], dynamic("{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.fa"))
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids [12b]
rule unmapped_masked:
    input:
      os.path.join(config["outdir"], dynamic("{sample}/12a_unmapped_reads/RefGenome_unmapped.{n}.fa")),
      os.path.join(config["outdir"], dynamic("{sample}/10_repeatmasker_good/masked.{n}.fa"))
    output:
      os.path.join(config["outdir"], dynamic("{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa"))
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmapped_masked_ids.py"

## MegaBlast against reference genome to remove more host sequences [13]
rule megablast_ref_genome:
    input:
      db = config["ref_genome"],
      query = os.path.join(config["outdir"], dynamic("{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa"))
    output:
      os.path.join(config["outdir"], dynamic("{sample}/13_megablast/megablast.{n}.xml"))
    params:
      perc_ident = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      num_desc = config["megablast_ref_genome"]["num_descriptions"],
      num_align = config["megablast_ref_genome"]["num_alignments"]
    threads: 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/megablast_ref_genome.py"

## Filter megablast records for the cutoff value [14]
rule parse_megablast:
    input:
      blastxml = os.path.join(config["outdir"], dynamic("{sample}/13_megablast/megablast.{n}.xml")),
      query = os.path.join(config["outdir"], dynamic("{sample}/12b_unmapped_masked/RefGenome_unmapped.{n}.masked.fa"))
    output:
      known = os.path.join(config["outdir"], dynamic("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.non-viral.out")),
      unmapped = os.path.join(config["outdir"], dynamic("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa"))
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

## Blast against virus database [15]
rule blastn_virus_nt:
    input:
      db = config["virus_nt"],
      query = os.path.join(config["outdir"], dynamic("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa"))
    output:
      out = os.path.join(config["outdir"], dynamic("{sample}/15_blast_virusnt/blast_virusnt.{n}.xml"))
    params:
      task = "blastn",
      show_gis = True,
      evalue = 1e-4,
      db_soft_mask = 100,
      num_threads = 8
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blastn_virus_db.py"

## Filter blastn records for the cutoff value [16]
rule parse_virusntblast:
    input:
      blastxml = os.path.join(config["outdir"], dynamic("{sample}/15_blast_virusnt/blast_virusnt.{n}.xml")),
      query = os.path.join(config["outdir"], dynamic("{sample}/14_megablast_parsed/RefGenome_megablast.{n}.unmapped.fa"))
    output:
      known = os.path.join(config["outdir"], dynamic("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.known-viral.out")),
      unmapped = os.path.join(config["outdir"], dynamic("{sample}/16_blastntvirus_parsed/blastnt_virus.{n}.unmapped.fa"))
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast_xml.py"

