
library(readr)
library(dplyr)
library(stringr)

download.file("ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz",
              destfile = "data/new_taxdump.tar.gz")

untar("data/new_taxdump.tar.gz", files = "names.dmp", exdir = "data")
untar("data/new_taxdump.tar.gz", files = "nodes.dmp", exdir = "data")

# Import names for tax_id
# Taxonomy names file has these fields:
#
# tax_id -- the id of node associated with this name
# name_txt -- name itself
# unique name -- the unique variant of this name if name not unique
# name class -- (synonym, common name, ...)

names <- read_delim("data/names.dmp", delim = "|",
                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names <- mutate_all(names, str_replace_all, "\\t", "")
names <- select(names, -5)
colnames(names) <- c("tax_id", "name_txt", "unique_name", "name_class")
names <- mutate_at(names, "tax_id", parse_integer)
write_csv(names, "data/names.csv")

# Import taxonomy nodes
# This file represents taxonomy nodes. The description for each node includes
# the following fields:
#
# tax_id -- node id in GenBank taxonomy database
# parent tax_id -- parent node id in GenBank taxonomy database
# rank -- rank of this node (superkingdom, kingdom, ...)
# embl code -- locus-name prefix; not unique
# division id -- see division.dmp file
# inherited div flag  (1 or 0) -- 1 if node inherits division from parent
# genetic code id -- see gencode.dmp file
# inherited GC  flag  (1 or 0) -- 1 if node inherits genetic code from parent
# mitochondrial genetic code id -- see gencode.dmp file
# inherited MGC flag  (1 or 0) -- 1 if node inherits mitochondrial gencode from parent
# GenBank hidden flag (1 or 0) -- 1 if name is suppressed in GenBank entry lineage
# hidden subtree root flag (1 or 0) -- 1 if this subtree has no sequence data yet
# comments -- free-text comments and citations
# plastid genetic code id -- see gencode.dmp file
# inherited PGC flag  (1 or 0) -- 1 if node inherits plastid gencode from parent
# specified_species -- 1 if species in the node's lineage has formal name
# hydrogenosome genetic code id -- see gencode.dmp file
# inherited HGC flag  (1 or 0) -- 1 if node inherits hydrogenosome gencode from parent

nodes <- read_delim("data/nodes.dmp", delim = "|",
                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
nodes <- mutate_all(nodes, str_replace_all, "\\t", "")
nodes <- select(nodes, -19)
colnames(nodes) <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                     "inherited_div_flag", "genetic_code_id", "inherited_GC_flag",
                     "mitochondrial_genetic_code_id", "inherited_MGC_flag", "GenBank_hidden_flag",
                     "hidden_subtree_root_flag", "comments", "plastid_genetic_code_id",
                     "inherited_PGC_flag", "specified_species", "hydrogenosome_genetic_code_id",
                     "inherited_HGC_flag")
nodes <- mutate_at(nodes, vars(ends_with("id"), ends_with("flag")), parse_integer)
write_csv(nodes, "data/nodes.csv")

# Remove dmp files and archive
file.remove("data/names.dmp", "data/nodes.dmp", "data/new_taxdump.tar.gz")
