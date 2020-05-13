
# Helper function to import tables
def safely_read_csv(path, **kwargs):
  try:
    return pd.read_csv(path, **kwargs)
  except pd.errors.EmptyDataError:
    pass

RANKS_OF_INTEREST = ['superkingdom', 'order', 'family', 'genus', 'species']

def concatenate_tables(input, sep = "\s+", cols_to_integer = None):
  frames = [safely_read_csv(f, sep = sep) for f in input]
  frames_concatenated = pd.concat(frames, keys = input, sort=False)
  if cols_to_integer:
    frames_concatenated[cols_to_integer] = frames_concatenated[cols_to_integer].apply(lambda x: pd.Series(x, dtype = "Int64"))
  return(frames_concatenated)

rule get_virus_taxids:
    output: 
        "output/blast/10239.taxids"
    params:
       taxid = 10239
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/blast/query/environment.yaml"
    resources:
        runtime = 20
    shell:
       "get_species_taxids.sh -t {params.taxid} > {output}"

rule get_negative_taxids:
    output: 
        temp("output/blast/{taxid}.taxids")
    params:
       taxid = lambda wildcards: wildcards.taxid
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/blast/query/environment.yaml"
    resources:
        runtime = 20
    shell:
       "get_species_taxids.sh -t {params.taxid} > {output}"


rule merge_taxidlists:
    input:
        expand("output/blast/{taxid}.taxids", taxid = [9606, 2, 12908])
    output:
        "output/blast/negative.taxids"
    resources:
        runtime = 10
    shell:
        "cat {input} > {output}"


rule megablast_virus:
    input:
        query = rules.parse_megablast_host.output.unmapped,
        taxidlist = "output/blast/10239.taxids"
    output:
        out = temp("output/{run}/megablast-virus_{n}.tsv")
    params:
        program = "megablast",
        db = "nt_v5",
        word_size = 16,
        evalue = 10,
        max_hsps = 50,
        outfmt = "'6 qseqid sacc staxid pident length evalue'"
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 480,
        mem_mb = 4000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_megablast_virus:
    input:
        query = rules.parse_megablast_host.output.unmapped,
        blast_result = rules.megablast_virus.output.out
    output:
        mapped = temp("output/{run}/megablast-virus_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/megablast-virus_{n}_unmapped.fa")
    params:
        e_cutoff = 5,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST


# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
        query = rules.parse_megablast_virus.output.unmapped,
        taxidlist = "output/blast/10239.taxids"
    output:
        out = temp("output/{run}/blastn-virus_{n}.tsv")
    params:
        program = "blastn",
        db = "nt_v5",
        max_hsps = 50,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 480,
        mem_mb = 4000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
        query = rules.parse_megablast_virus.output.unmapped,
        blast_result = rules.blastn_virus.output.out
    output:
        mapped = temp("output/{run}/blastn-virus_{n}_hits.tsv"),
        unmapped = temp("output/{run}/blastn-virus_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-5,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST

# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
        query = rules.parse_blastn_virus.output.unmapped,
        taxidlist = "output/blast/10239.taxids"
    output:
        out = temp("output/{run}/blastx-virus_{n}.tsv")
    params:
        program = "blastx",
        task = "Blastx-fast",
        db = "nr_v5",
        evalue = 1e-2,
        db_soft_mask = 100,
        max_hsps = 50,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 480,
        mem_mb = 4000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_blastx_virus:
    input:
        query = rules.blastx_virus.input.query,
        blast_result = rules.blastx_virus.output.out
    output:
        mapped = temp("output/{run}/blastx-virus_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/blastx-virus_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-3,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST

# Filter sequences by division id.
# Saves hits with division id
# pp_sway - maximum difference of hit from top hit pident to consider for consensus taxonomy, percent points.
# ranks_of_interest - taxonomic ranks to be included in results. E.g, subspecies will be demoted to species.
# dbfile - path to etetoolkit taxon.sqlite database. Default location is $HOME/.etetoolkit/taxon.sqlite. Plese see http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html for more information.
rule classify_viruses:
    input:
        [rules.parse_blastn_virus.output.mapped, rules.parse_blastx_virus.output.mapped] if config["run_blastx"] else rules.parse_blastn_virus.output.mapped
    output:
        temp("output/{run}/viruses_{n}.csv")
    params:
        pp_sway = 1,
        ranks_of_interest = RANKS_OF_INTEREST,
        dbfile = TAXON_DB
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        BLAST_TAXONOMY

# Filter unmasked candidate virus reads.
rule unmasked_other:
    input:
        rules.parse_blastx_virus.output.unmapped if config["run_blastx"] else rules.parse_blastn_virus.output.unmapped,
        rules.repeatmasker_good.output.original_filt
    output:
        temp("output/{run}/candidate-viruses_{n}_unmasked.fa")
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        SUBSET_FASTA


# Map reads to bacterial genomes
rule mapbact:
    input:
        input = rules.unmasked_other.output,
        ref = REF_BACTERIA
    output:
        outu = temp("output/{run}/unmapbact_{n}.fa"),
        outm = temp("output/{run}/mapbact_{n}.fa"),
        statsfile = "output/{run}/mapbact_{n}.txt"
    params:
        extra = "nodisk -Xmx60g"
    resources:
        runtime = lambda wildcards, attempt: attempt * 240,
        mem_mb = 60000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbwrap"


# Subset repeatmasker masked reads using unmapped reads.
rule refbac_unmapped_masked:
    input:
        rules.mapbact.output.outu,
        rules.repeatmasker_good.output.masked_filt
    output:
        temp("output/{run}/unmapped_{n}_masked.fa")
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        SUBSET_FASTA

# Megablast against nt database.
rule megablast_nt:
    input:
        query = rules.refbac_unmapped_masked.output,
        negative_taxidlist = "output/blast/negative.taxids"
    output:
         out = temp("output/{run}/megablast-nt_{n}.tsv")
    params:
        program = "blastn",
        db = "nt_v5",
        task = "megablast",
        evalue = 1e-8,
        word_size = 16,
        max_hsps = 50,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 240,
        mem_mb = lambda wildcards, attempt: 20000 + (attempt * 20000)
    wrapper:
        BLAST_QUERY

# Filter megablast hits for the cutoff value.
rule parse_megablast_nt:
    input:
        query = rules.refbac_unmapped_masked.output,
        blast_result = rules.megablast_nt.output.out
    output:
        mapped = temp("output/{run}/megablast-nt_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/megablast-nt_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-10,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120
    wrapper:
        PARSE_BLAST

# Blastn against nt database.
rule blastn_nt:
    input:
        query = rules.parse_megablast_nt.output.unmapped,
        negative_taxidlist = "output/blast/negative.taxids"
    output:
        out = temp("output/{run}/blastn-nt_{n}.tsv")
    params:
        program = "blastn",
        db = "nt_v5",
        task = "blastn",
        evalue = 1e-8,
        max_hsps = 50,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 240,
        mem_mb = 40000
    wrapper:
        BLAST_QUERY

# Filter blastn records for the cutoff value.
rule parse_blastn_nt:
    input:
        query = rules.blastn_nt.input.query,
        blast_result = rules.blastn_nt.output.out
    output:
        mapped = temp("output/{run}/blastn-nt_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/blastn-nt_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-10,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120
    wrapper:
        PARSE_BLAST

# Blastx unmapped sequences against nr database.
rule blastx_nr:
    input:
        query = rules.parse_blastn_nt.output.unmapped,
        negative_taxidlist = "blast/negative.taxids"
    output:
        out = temp("output/{run}/blastx-nr_{n}.tsv")
    params:
        program = "blastx",
        task = "Blastx-fast",
        db = "nr_v5",
        evalue = 1e-2,
        max_hsps = 50,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: attempt * 480,
        mem_mb = 40000
    wrapper:
      BLAST_QUERY

# Filter blastx records for the cutoff value.
rule parse_blastx_nr:
    input:
        query = rules.blastx_nr.input.query,
        blast_result = rules.blastx_nr.output.out
    output:
        mapped = temp("output/{run}/blastx-nr_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/blastx-nr_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-3,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120
    wrapper:
        PARSE_BLAST

# Filter sequences by division id.
# Saves hits with division id
rule classify_all:
    input:
        expand("output/{{run}}/{blastresult}_{{n}}_mapped.tsv", blastresult = BLASTNR)
    output:
        temp("output/{run}/all_{n}.csv")
    params:
        pp_sway = 1,
        ranks_of_interest = RANKS_OF_INTEREST,
        dbfile = TAXON_DB
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
    wrapper:
        BLAST_TAXONOMY

# Split classification rule outputs into viruses and non-viral
rule filter_viruses:
    input:
        expand("output/{{run}}/{classified}_{n}.csv", n = N, classified = ["viruses", "all"])
    output:
        viral = "output/{run}/viruses.csv",
        non_viral = "output/{run}/non-viral.csv"
    params:
        ranks = RANKS_OF_INTEREST
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    run:
        tab = concatenate_tables(input, sep = ",", cols_to_integer = params.ranks)
        vir = tab[tab.superkingdom == 10239]
        non_vir = tab[tab.superkingdom != 10239]
        vir.to_csv(output.viral, index = False)
        non_vir.to_csv(output.non_viral, index = False)


# Merge unassigned sequences
rule merge_unassigned:
    input:
        expand("output/{{run}}/blast{type}_{n}_unmapped.fa", type = "x-nr" if config["run_blastx"] else "n-nt", n = N)
    output:
        "output/{run}/unassigned.fa"
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
    shell:
        "cat {input} > {output}"


# Upload results to Zenodo.
if config["zenodo"]["deposition_id"]:
    rule upload_results:
        input:
            expand("output/{{run}}/{result}", result = RESULTS)
        output:
            ZEN.remote("output/{run}/counts.tgz")
        resources:
            runtime = lambda wildcards, attempt: attempt * 120,
        shell:
            "tar czvf {output} {input}"
