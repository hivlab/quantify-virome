
## Blast input, output, and params keys must match commandline blast option names https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a

## Blast against NT virus database
rule blastn_virus:
    input:
      db = config["virus_nt"],
      query = preprocessing("output/{sample}_refgenome_filtered_{n}_unmapped.fa")
    output:
      out = "output/blast/{sample}_blastn_virus_{n}.xml"
    params:
      task = "blastn",
      evalue = 1e-4,
      db_soft_mask = 100,
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastn_virus:
    input:
      rules.blastn_virus.output.out,
      preprocessing("output/{sample}_refgenome_filtered_{n}_unmapped.fa")
    output:
      known_xml = "output/{sample}_blastn_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastn_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastx unmapped sequences against NR virus database
rule blastx_virus:
    input:
      db = config["virus_nr"],
      query = rules.parse_blastn_virus.output.unmapped
    output:
      out = "output/blast/{sample}_blastx_virus_{n}.xml"
    params:
      evalue = 1e-2,
      db_soft_mask = 100,
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastx_virus:
    input:
      rules.blastx_virus.output.out,
      rules.blastx_virus.input.query
    output:
      known_xml = "output/{sample}_blastx_virus_{n}_known-viral.xml",
      unmapped = "output/{sample}_blastx_virus_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Get unmasked known viral sequences
rule unmasked_viral:
  input:
    rules.parse_blastn_virus.output.known_xml,
    rules.parse_blastx_virus.output.known_xml,
    preprocessing("output/{sample}_refgenome_unmapped_{n}.fa")
  output:
    "output/{sample}_blastn_virus_{n}_known-viral.fa",
    "output/{sample}_blastx_virus_{n}_known-viral.fa"
  conda:
      "../envs/biopython.yml"
  script:
      "../scripts/unmasked_viral.py"

## Merge blast outputs
rule merge_unmasked_viral:
  input:
    rules.unmasked_viral.output
  output:
    "output/{sample}_known-viral_{n}_unmasked.fa"
  shell:
    """
    cat {input} > {output}
    """
## Map against bacterial genomes
rule bwa_mem:
    input:
        config["ref_bacteria"],
        ["output/{sample}_known-viral_{n}_unmasked.fa"]
    output:
        "output/bwa_mem/{sample}_refbacteria_mapped_{n}.bam"
    log:
        "output/logs/{sample}_bactbwa_mem_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
        "(bwa mem -L 100,100 -k 15 -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

## Extract unmapped reads
rule unmapped_reads:
    input: rules.bwa_mem.output
    output:
      bam = temp("output/{sample}_refbacteria_unmapped_{n}.bam"),
      fq = temp("output/{sample}_refbacteria_unmapped_{n}.fq"),
      fa = "output/{sample}_refbacteria_unmapped_{n}.fa"
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        samtools view -b -f 4 {input} > {output.bam}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """

## Subset repeatmasker masked reads using unmapped ids
rule unmapped_masked:
    input: rules.unmapped_reads.output.fa, preprocessing("output/{sample}_repmaskedgood_{n}.fa")
    output:
      temp("output/{sample}_refbacteria_unmapped_{n}_masked.fa")
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/unmapped_masked_ids.py"

## MegaBlast against NT to remove host sequences
rule megablast_nt:
    input:
      db = config["nt"],
      query = rules.unmapped_masked.output
    output:
      out = "output/blast/{sample}_megablast_nt_{n}.xml"
    params:
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      num_descriptions = config["megablast_nt"]["num_descriptions"],
      num_alignments = config["megablast_nt"]["num_alignments"],
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter megablast records for the cutoff value
rule parse_megablast:
    input:
      rules.megablast_nt.output.out,
      rules.unmapped_masked.output
    output:
      "output/blast/{sample}_nt_filtered_{n}_mapped.xml",
      "output/blast/{sample}_nt_filtered_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastn against NT database
rule blastn_nt:
    input:
      db = config["nt"],
      query = "output/blast/{sample}_nt_filtered_{n}_unmapped.fa"
    output:
      out = "output/blast/{sample}_blastn_nt_{n}.xml"
    params:
      task = "blastn",
      evalue = 1e-8,
      num_descriptions = 50,
      num_alignments = 50,
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastn_nt:
    input:
      rules.blastn_nt.output.out,
      rules.blastn_nt.input.query
    output:
      known_xml = "output/{sample}_blastn_nt_{n}_mapped.xml",
      unmapped = "output/{sample}_blastn_nt_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-10
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"

## Blastx unmapped sequences against NR virus database
rule blastx_nr:
    input:
      db = config["nr"],
      query = rules.parse_blastn_nt.output.unmapped
    output:
      out = "output/blast/{sample}_blastx_nr_{n}.xml"
    params:
      evalue = 1e-2,
      num_descriptions = 50,
      num_alignments = 50,
      show_gis = True,
      num_threads = 8,
      outfmt = 5
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/blast.py"

## Filter blastn records for the cutoff value
rule parse_blastx_nr:
    input:
      rules.blastx_nr.output.out,
      rules.blastx_nr.input.query
    output:
      known_xml = "output/{sample}_blastx_nr_{n}_mapped.xml",
      unassigned = "output/{sample}_blastx_nr_{n}_unmapped.fa"
    params:
      e_cutoff = 1e-3
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/parse_blast.py"
