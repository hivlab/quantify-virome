
source("../scripts/common/helpers.R")

blast_taxonomy <- function(..., taxdb, nodes, phages, viruses, div_id = 3) {

  message("Parse blast results xml\n")
  dots <- c(...)
  xmls <- grep("xml$", dots, value = TRUE)
  tab <- parse_blast_xml_safe(xmls)

  # if it's not empty result then stop error
  if (!is.null(tab$error) && !str_detect(tab$error$message, "No BLAST hits to parse")) stop(tab$error)

  message("Map tax_ids to gis\n")
  tab <- gi2taxid(tab$result, taxdb)

  message("Filter phages (division_id == 3) and save results to csv files\n")
  filter_division(tab = tab, nodes = nodes, div = phages, not_div = viruses, div_id = div_id)
}

do.call(blast_taxonomy, c(snakemake@input, snakemake@output))

