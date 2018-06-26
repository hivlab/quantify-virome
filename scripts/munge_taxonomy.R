library(tidyverse)

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
  out_data <- separate(out_data, stats, c("evalue", "bitscore"), sep = ";") %>%
    mutate_at(vars(evalue, bitscore), str_remove, "^[[:alpha:] ]+") %>%
    mutate_at(vars(evalue, bitscore), parse_double)
  return(out_data)
}

known <- mutate(known, data = map(lines, get_out))
known_unnested <- select(known, files, data) %>%
  unnest() %>%
  mutate(gi = parse_integer(gi))
db_path <- "/gpfs/software/VirusSeeker/databases/taxdump_300817/vhunter.db"
db <- src_sqlite(db_path, create = TRUE)
nucl_db <- tbl(db, "gi_taxid_nucl")
gi_tab <- nucl_db %>%
  filter(gi %in% known_unnested$gi) %>%
  collect()

known_unnested <- left_join(known_unnested, gi_tab)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomeInfoDbData")
library(GenomeInfoDbData)
data(specData)
known_unnested <- left_join(known_unnested, specData)
known_unnested %>%
  mutate_at(vars(genus), str_replace_all, c("[^[:alnum:]]$" = "",  "s$" = "", "(\\(\\d*)" = "\\1\\)" )) %>%
  mutate_at(vars(genus, species), str_to_lower) %>%
  select(tax_id, genus, species) %>%
  distinct()
names <- read_delim("~/databases/taxdump/names.dmp", delim = "|",
                    escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
names <- names %>%
  mutate_all(str_remove_all, "\\t")
names <- select(names, -5)
colnames(names) <- c("tax_id", "name_txt", "unique name", "name class")
names
names <- mutate_at(names, "tax_id", parse_integer)
known_names <- left_join(known_unnested, names)
known_names_ids <- select(known_names, query, gi, tax_id, name_txt, `unique name`, `name class`)
unique_taxa <- filter(known_names_ids, str_detect(`name class`, "scientific")) %>%
  select(-query) %>%
  distinct()
unique_taxa
