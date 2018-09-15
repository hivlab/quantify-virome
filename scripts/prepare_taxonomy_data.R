library(readr)
library(dplyr)
library(stringr)

taxonomy_db <- function (names_dmp, nodes_dmp, division_dmp, names_out, nodes_out, division_out) {
  # Import names for tax_id
  # Taxonomy names file has these fields:
  #
  # tax_id -- the id of node associated with this name
  # name_txt -- name itself
  # unique name -- the unique variant of this name if name not unique
  # name class -- (synonym, common name, ...)

  names <- read_delim(names_dmp, delim = "|",
                      escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  names <- mutate_all(names, str_replace_all, "\\t", "")
  names <- select(names, which(vapply(names, function(x) !all(is.na(x)), logical(1))))
  stopifnot(ncol(names) == 4)
  colnames(names) <- c("tax_id", "name_txt", "unique_name", "name_class")
  names <- mutate_at(names, "tax_id", parse_integer)
  write_csv(names, names_out)

  # Import taxonomy nodes
  # This file represents taxonomy nodes. The description for each node includes
  # the following fields:
  #
  # tax_id					-- node id in GenBank taxonomy database
  # parent tax_id				-- parent node id in GenBank taxonomy database
  # rank					-- rank of this node (superkingdom, kingdom, ...)
  # embl code				-- locus-name prefix; not unique
  # division id				-- see division.dmp file
  # inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
  # genetic code id				-- see gencode.dmp file
  # inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
  # mitochondrial genetic code id		-- see gencode.dmp file
  # inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
  # GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
  # hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
  # comments				-- free-text comments and citations

  nodes <- read_delim(nodes_dmp, delim = "|",
                      escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  nodes <- mutate_all(nodes, str_replace_all, "\\t", "")
  nodes <- select(nodes, which(vapply(nodes, function(x) !all(is.na(x)), logical(1))))
  stopifnot(ncol(nodes) == 13)
  colnames(nodes) <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                       "inherited_div_flag", "genetic_code_id", "inherited_GC_flag",
                       "mitochondrial_genetic_code_id", "inherited_MGC_flag", "GenBank_hidden_flag",
                       "hidden_subtree_root_flag", "comments")
  nodes <- mutate_at(nodes, vars(ends_with("id"), ends_with("flag")), parse_integer)
  write_csv(nodes, nodes_out)

  # division.dmp
  # ------------
  # Divisions file has these fields:
  # division id				-- taxonomy database division id
  # division cde				-- GenBank division code (three characters)
  # division name				-- e.g. BCT, PLN, VRT, MAM, PRI...
  # comments

  division <- read_delim(division_dmp, delim = "|",
                         escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  division <- mutate_all(division, str_replace_all, "\\t", "")
  division <- select(division, which(vapply(division, function(x) !all(is.na(x)), logical(1))))
  stopifnot(ncol(division) == 4)
  colnames(division) <- c("division_id", "division_cde", "division_name", "comments")
  write_csv(division, division_out)
}

taxonomy_db(snakemake@input[[1]],
            snakemake@input[[2]],
            snakemake@input[[3]],
            snakemake@output[[1]],
            snakemake@output[[2]],
            snakemake@output[[3]])
