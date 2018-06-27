
setClass("snakemake", representation(input = "list", output = "list", params = "list"))
snakemake <- new("snakemake",
                 input = list(known = "output/I1164_12629_Harvard_SIV_196_06_2_24_12/16_blastntvirus_parsed/",
                              vhunter = "/gpfs/software/VirusSeeker/databases/taxdump_300817/vhunter.db",
                              names = "~/databases/taxdump/names.dmp",
                              nodes = "~/databases/taxdump/nodes.dmp"),
                 output = list(known_taxa = "output/I1164_12629_Harvard_SIV_196_06_2_24_12/17_munge_taxonomy/known_taxa.csv"))

library(tidyverse)

# Import blast hits and munge to data frame
files <- list.files(snakemake@input[["known"]], full.names = TRUE)
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

# Merge blast output with taxonomy info
# Import nucleotide gi to tax_id data
db_path <- snakemake@input[["vhunter"]]
db <- src_sqlite(db_path, create = TRUE)
nucl_db <- tbl(db, "gi_taxid_nucl")
gi_tab <- nucl_db %>%
  filter(gi %in% known_unnested$gi) %>%
  collect()

# Join taxonomy by gi
known_unnested <- left_join(known_unnested, gi_tab)



# Merge names by tax_id to known sequences
known_names <- left_join(known_unnested, names)
known_taxa <- left_join(known_names, nodes)
write_csv(known_taxa, snakemake@output[["known_taxa"]])
