
source("../scripts/helpers.R")

blast_taxonomy <- function(nt_virus, nr_virus, taxdb, nodes, phages, not_phages, div_id = 3) {
  # parse blast results xml
  tab <- parse_blast_xml(nt_virus, nr_virus)
  # map tax_ids to gis
  tab <- gi2taxid(tab, taxdb)
  # filter phages (division_id == 3) and save results to csv files
  filter_division(tab = tab, nodes = nodes, div = phages, not_div = not_phages, div_id = div_id)
}

blast_taxonomy(nt_virus = snakemake@input[[1]],
               nr_virus = snakemake@input[[2]],
               taxdb = snakemake@input[[3]],
               nodes = snakemake@input[[4]],
               phages = snakemake@output[[1]],
               not_phages = snakemake@output[[2]])
