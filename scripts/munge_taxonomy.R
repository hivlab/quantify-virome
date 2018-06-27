
library(tidyverse)

# Import blast hits and munge to data frame
files <- list.files("output/I1164_12629_Harvard_SIV_196_06_2_24_12/16_blastntvirus_parsed/", full.names = TRUE)
known <- data_frame(files) %>%
  filter(str_detect(files, "known"))
known <- mutate(known, lines = map(files, read_lines))

get_out <- function(out) {
  query <- str_trim(out[str_detect(out, "Query:")])
  query <- str_trim(str_remove(query, "Query:"))
  hit <- str_trim(out[str_detect(out, "Hit:")])
  hit <- str_trim(str_remove(hit, "Hit:"))
  stats <- str_trim(out[str_detect(out, "Quick stats:")])
  stats <- str_trim(str_remove(stats, "Quick stats:"))
  out_data <- data_frame(query, hit, stats)
  out_data <- separate(out_data, hit, c("out1","gi","out2","gb", "description"), sep = "\\|") %>%
    select(-starts_with("out")) %>%
    mutate_at("description", str_trim)
  separate(out_data, stats, c("evalue", "bitscore"), sep = ";") %>%
    mutate_at(vars(evalue, bitscore), str_remove, "^[[:alpha:] ]+") %>%
    mutate_at(vars(evalue, bitscore), parse_double)
}

known <- mutate(known, data = map(lines, get_out))
known_unnested <- select(known, files, data) %>%
  unnest() %>%
  mutate(gi = parse_integer(gi))

#
db_path <- "/gpfs/software/VirusSeeker/databases/taxdump_300817/vhunter.db"
db <- src_sqlite(db_path, create = TRUE)
nucl_db <- tbl(db, "gi_taxid_nucl")
gi_tab <- nucl_db %>%
  filter(gi %in% known_unnested$gi) %>%
  collect()

known_unnested <- left_join(known_unnested, gi_tab)

# Import names for tax_id
# Taxonomy names file has these fields:
#
# tax_id					-- the id of node associated with this name
# name_txt				-- name itself
# unique name				-- the unique variant of this name if name not unique
# name class				-- (synonym, common name, ...)
names <- read_delim("~/databases/taxdump/names.dmp", delim = "|",
                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names <- names %>%
  mutate_all(str_remove_all, "\\t")
names <- select(names, -5)
colnames(names) <- c("tax_id", "name_txt", "unique_name", "name_class")
names <- mutate_at(names, "tax_id", parse_integer)

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
# plastid genetic code id                 -- see gencode.dmp file
# inherited PGC flag  (1 or 0)            -- 1 if node inherits plastid gencode from parent
# specified_species			-- 1 if species in the node's lineage has formal name
# hydrogenosome genetic code id           -- see gencode.dmp file
# inherited HGC flag  (1 or 0)            -- 1 if node inherits hydrogenosome gencode from parent

nodes <- read_delim("~/databases/taxdump/nodes.dmp", delim = "|",
                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
nodes <- nodes %>%
  mutate_all(str_remove_all, "\\t")
nodes <- select(nodes, -19)
colnames(nodes) <- c("tax_id", "parent_tax_id", "rank", "embl_code", "division_id",
                     "inherited_div_flag", "genetic_code_id", "inherited_GC_flag",
                     "mitochondrial_genetic_code_id", "inherited_MGC_flag", "GenBank_hidden_flag",
                     "hidden_subtree_root_flag", "comments", "plastid_genetic_code_id",
                     "inherited_PGC_flag", "specified_species", "hydrogenosome_genetic_code_id",
                     "inherited_HGC_flag")
nodes <- mutate_at(nodes, vars(ends_with("id"), ends_with("flag")), parse_integer)

# Merge names by tax_id to known sequences
known_names <- left_join(known_unnested, names)
known_taxa <- left_join(known_names, nodes)

known_taxa_ids <- select(known_taxa, query, gi, tax_id, parent_tax_id, rank,
                          division_id, inherited_div_flag,
                          name_txt, unique_name, name_class)

# Keep only scientific names
unique_taxa <- filter(known_taxa_ids, str_detect(name_class, "scientific"))

# Number of unique queries per tax_id
unique_taxa %>%
  group_by(tax_id, name_txt) %>%
  summarise(N = n()) %>%
  ungroup() %>%
  arrange(desc(N))

# Drop query and get distinct gi-s
(unique_taxa <- unique_taxa %>%
  select(-query) %>%
  mutate_if(is.character, str_to_lower) %>%
  distinct() %>%
  select(gi, tax_id, name_txt, rank, everything())
  )

# Add parent names
parent_names <- filter(names, tax_id %in% unique_taxa$parent_tax_id,
                       str_detect(name_class, "scientific")) %>%
  select(tax_id, name_txt, unique_name)
colnames(parent_names) <- c("parent_tax_id", "parent_name_txt", "parent_unique_name")
unique_taxa <- left_join(unique_taxa, parent_names)

# Let's try to build a taxonomy tree
unique_taxa <- select(unique_taxa, -contains("div"), -name_class) %>%
  mutate_if(is.character, str_to_lower)

# Summarise the number of observations on each unique taxon and parent level
# First, number of rows
(unique_tax_id <- unique_taxa %>%
  group_by(parent_tax_id, parent_name_txt, tax_id, name_txt) %>%
  summarise(N = n()) %>%
  arrange(desc(N))
  )

# Approximate number of taxons belonging to virus/bacteriophage classes
unique_tax_id %>%
  mutate(what = case_when(
    (str_detect(parent_name_txt, "vir") & str_detect(parent_name_txt, "bac|esch|coccu|chla|mona|erw|lis|vibr|shi|yers|salm|phi|clos|chla")) ~ "phage",
    (str_detect(parent_name_txt, "vir") & !str_detect(parent_name_txt, "bac")) ~ "virus",
    str_detect(parent_name_txt, "phag") ~ "phage",
    TRUE ~ "other"
  )) %>%
  group_by(what) %>%
  summarise(N = n(),
            `%` = round((N / nrow(.)) * 100, 1))

# Approximate distribution of parent taxons to virus/bacteriophage classes
(unique_parent_id <- unique_tax_id %>%
    summarise(N = n()) %>%
    summarise(N = n()) %>%
    arrange(desc(N)))

unique_parent_id %>%
  mutate(what = case_when(
    (str_detect(parent_name_txt, "vir") & str_detect(parent_name_txt, "bac|esch|coccu|chla|mona|erw|lis|vibr|shi|yers|salm|phi|clos|chla")) ~ "phage",
    (str_detect(parent_name_txt, "vir") & !str_detect(parent_name_txt, "bac")) ~ "virus",
    str_detect(parent_name_txt, "phag") ~ "phage",
    TRUE ~ "other"
  )) %>%
  group_by(what) %>%
  summarise(N = n(),
            `%` = round((N / nrow(.)) * 100, 1))

