
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])


# Convert reads to interleaved format
rule interleave:
    input:
        lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
    output:
        out = temp("output/{run}/interleaved.fq.gz"),
        bhist = "output/{run}/bhist.txt",
        qhist = "output/{run}/qhist.txt",
        aqhist = "output/{run}/aqhist.txt",
        bqhist = "output/{run}/bqhist.txt",
        lhist = "output/{run}/lhist.txt",
        gchist = "output/{run}/gchist.txt"
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 20,
        mem_mb = 4000
    log:
        "output/{run}/log/interleave.txt"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/reformat"


# Remove PCR and optical duplicates
rule clumpify:
    input:
        input = rules.interleave.output.out
    output:
        out = temp("output/{run}/clumpify.fq.gz")
    params:
        extra = lambda wildcards, resources: f"dedupe optical -Xmx{resources.mem_mb / 1000:.0f}g -da" # -da -- suppress assertions
    resources:
        runtime = lambda wildcards, attempt: 20 + (attempt * 20),
        mem_mb = 4000
    log: 
        "output/{run}/log/clumpify.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/clumpify"


rule filterbytile:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{run}/filterbytile.fq.gz")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = 20,
        mem_mb = 4000
    log: 
        "output/{run}/log/filterbytile.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/filterbytile"


rule trim:
    input:
        input = rules.filterbytile.output.out
    output:
        out = temp("output/{run}/trimmed.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = lambda wildcards, attempt: 20 + (attempt * 20),
        mem_mb = 4000
    log: 
        "output/{run}/log/trim.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


rule artifacts:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{run}/filtered.fq.gz"
    params:
        extra = lambda wildcards, resources: f"k=31 ref=artifacts,phix ordered cardinality -Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = lambda wildcards, attempt: 20 + (attempt * 20),
        mem_mb = 4000
    log: 
        "output/{run}/log/artifacts.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


# Remove host sequences
rule maphost:
    input:
        input = rules.artifacts.output.out,
        ref = HOST_GENOME
    output:
        outu = "output/{run}/unmaphost.fq.gz",
        outm = "output/{run}/maphost.fq.gz",
        statsfile = "output/{run}/maphost.txt"
    params:
        extra = lambda wildcards, resources: f"nodisk -Xmx{resources.mem_mb / 1000:.0f}g"
    log: 
        "output/{run}/log/maphost.log"
    resources:
        runtime = lambda wildcards, attempt: 120 + (attempt * 60),
        mem_mb = 24000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


rule correct1:
    input:
        input = rules.maphost.output.outu
    output:
        out = temp("output/{run}/ecco.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ecco mix vstrict ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct1.log"
    resources:
        runtime = lambda wildcards, attempt: 60 + (attempt * 20),
        mem_mb = 4000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbmerge"


rule correct2:
    input:
        input = rules.correct1.output.out
    output:
        out = temp("output/{run}/eccc.fq.gz")
    params:
        extra = lambda wildcards, resources: f"passes=4 reorder -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct2.log"
    resources:
        runtime = lambda wildcards, attempt: 60 + (attempt * 20),
        mem_mb = 16000
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/clumpify"


rule correct3:
    input:
        input = rules.correct2.output.out
    output:
        out = temp("output/{run}/ecct.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ecc k=62 ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct3.log"
    resources:
        runtime = lambda wildcards, attempt: 60 + (attempt * 20),
        mem_mb = 16000
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/tadpole"


rule merge:
    input:
        input = rules.correct3.output.out
    output:
        out = temp("output/{run}/merged.fq.gz"),
        outu = temp("output/{run}/unmerged.fq.gz"),
        ihist = "output/{run}/ihist.txt"
    params:
        extra = lambda wildcards, resources: f"strict k=93 extend2=80 rem ordered -Xmx{resources.mem_mb / 1000:.0f}g"
    log: 
        "output/{run}/log/merge.log"
    resources:
        runtime = lambda wildcards, attempt: 20 + (attempt * 20),
        mem_mb = 16000
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbmerge"


rule qtrim:
    input:
        input = rules.merge.output.outu
    output:
        out = temp("output/{run}/qtrimmed.fq.gz")
    params:
        extra = lambda wildcards, resources: f"qtrim=r trimq=10 minlen=70 ordered -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = lambda wildcards, attempt: 20 + (attempt * 20),
        mem_mb = 4000
    log: 
        "output/{run}/log/qtrim.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


rule concatenate:
    input:
      rules.merge.output.out, rules.qtrim.output.out
    output:
      temp("output/{run}/concatenated.fq.gz")
    resources:
        runtime = 10
    shell:
      "cat {input} > {output}"


rule fastqtofasta:
    input:
        rules.concatenate.output[0]
    output:
        out = temp("output/{run}/concatenated.fa")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 20
    log:
        "output/{run}/log/fastqtofasta.txt"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/reformat"


# Run cd-hit to find and cluster duplicate sequences
rule cd_hit:
    input:
        rules.fastqtofasta.output.out
    output:
        repres = temp("output/{run}/cdhit.fa"),
        clstr = temp("output/{run}/cdhit.fa.clstr")
    params:
        extra = "-c 0.984 -G 0 -n 10 -d 0 -aS 0.984 -r 1 -M 0"
    log:
        "logs/{run}_cdhit.log"
    threads: 4
    resources:
        runtime = 2880,
        mem_mb = 14000
    wrapper:
        f"{WRAPPER_PREFIX}/master/cdhit"


# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
        rules.cd_hit.output.repres
    output:
        temp("output/{run}/tantan.fa")
    params:
        extra = "-x N" # mask low complexity using N
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        f"{WRAPPER_PREFIX}/master/tantan"


# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
    input:
        masked = rules.tantan.output[0]
    output:
        masked_filt = temp("output/{run}/tantangood.fa")
    params:
        min_length = 50,
        por_n = 40
    resources:
        runtime = 120
    wrapper:
        f"{WRAPPER_PREFIX}/master/filter/masked"


# Split reads to smaller chunks for Repeatmasker
rule split_fasta:
    input:
        rules.tantan_good.output.masked_filt
    output:
        temp(expand("output/{{run}}/repeatmasker_{n}.fa", n = N))
    params:
        config["split_fasta"]["n_files"]
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30) 
    wrapper:
        f"{WRAPPER_PREFIX}/master/split-fasta"


# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
    input:
        fa = "output/{run}/repeatmasker_{n}.fa"
    output:
        masked = temp("output/{run}/repeatmasker_{n}.fa.masked"),
        out = temp("output/{run}/repeatmasker_{n}.fa.out"),
        cat = temp("output/{run}/repeatmasker_{n}.fa.cat"),
        tbl = "output/{run}/repeatmasker_{n}.fa.tbl"
    params:
        extra = "-qq"
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 1440,
        mem_mb = 16000
    singularity:
        "shub://tpall/repeatmasker-singularity"
    script:
        f"{WRAPPER_PREFIX}/master/repeatmasker/wrapper.py"


# Filter repeatmasker output
# 1) Sequences > min_length nt of consecutive sequence without N
# 2) Sequences with >= % of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
    input:
        masked = rules.repeatmasker.output.masked,
        original = rules.repeatmasker.input.fa
    output:
        masked_filt = temp("output/{run}/repmaskedgood_{n}.fa"),
        original_filt = temp("output/{run}/unmaskedgood_{n}.fa")
    params:
        min_length = 100,
        por_n = 30
    resources:
        runtime = 120
    wrapper:
        f"{WRAPPER_PREFIX}/master/filter/masked"


# MegaBlast against reference genome to remove host sequences
rule megablast_host:
    input:
        query = rules.repeatmasker_good.output.masked_filt
    output:
        out = temp("output/{run}/megablast-host_{n}.tsv")
    shadow: "minimal"
    params:
        program = "blastn",
        task = "megablast",
        db = HOST_GENOME,
        evalue = 1e-6,
        max_hsps = 1,
        outfmt = "'6 qseqid sseqid pident length evalue'"
    threads: 8
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30)
    wrapper:
        BLAST_QUERY


# Filter megablast records for the cutoff value
rule parse_megablast_host:
    input:
        blast_result = rules.megablast_host.output.out,
        query = rules.repeatmasker_good.output.masked_filt
    output:
        mapped = "output/{run}/megablast-host_{n}_hits.tsv",
        unmapped = temp("output/{run}/megablast-host_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-9,
        outfmt = rules.megablast_host.params.outfmt
    resources:
        runtime = 120
    wrapper:
        PARSE_BLAST


# QC
fastq_screen_config = {
    "database": {
        "human": HOST_GENOME,
        "SILVA_138_SSU_132_LSU": RRNA_DB,
        "cpn60": CPNDB
    }
}

rule fastq_screen:
    input:
        rules.artifacts.output.out
    output:
        txt = "output/{run}/fastq_screen.txt",
        png = "output/{run}/fastq_screen.png"
    params:
        fastq_screen_config = fastq_screen_config,
        subset = 100000
    resources:
        runtime = 120,
        mem_mb = 8000    
    threads: 4
    wrapper:
        f"{WRAPPER_PREFIX}/master/fastq_screen"


rule fastqc:
    input:
        rules.interleave.output.out
    output:
        html = "output/{run}/fastqc.html",
        zip = "output/{run}/fastqc.zip"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


rule multiqc:
    input:
        "output/{run}/fastq_screen.txt",
        "output/{run}/fastqc.zip",
        "output/{run}/maphost.txt",
        expand("output/{{run}}/mapbact_{n}.txt", n = N)
    output:
        report("output/{run}/multiqc.html", caption = "report/multiqc.rst", category = "Quality control")
    log:
        "output/{run}/log/multiqc.log"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
      f"{WRAPPER_PREFIX}/master/multiqc"
