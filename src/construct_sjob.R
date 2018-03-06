library(yaml)
library(glue)
library(tidyverse)

input <- yaml_load_file_utf8("sbatch.yaml")
input

parse_key_value_pairs <- function(x) {
  key <- names(x)
  value <- x
  data_frame(key, value) %>%
    mutate(value = map(value, as.character)) %>% 
    unnest(value)
}

parse_sbatch_yaml <- function(input) {
  input <- yaml_load_file_utf8(input)
  nodes <- names(input)
  sbat <- input[["SBATCH"]]
  mod <- input[["modules"]]
  srun <- input["srun"]
  cm <- input[setdiff(nodes, c("SBATCH", "modules", "srun"))]
  sbat <- parse_key_value_pairs(sbat)
  mod <- parse_key_value_pairs(mod)
  srun <- parse_key_value_pairs(srun)
  cm <- parse_key_value_pairs(cm)

  sbatch <- glue_data(sbat, "#SBATCH --{key}={value}") %>% glue::collapse(sep = "\n")
  modules <- glue_data(mod, "modules {key} {value}") %>% glue::collapse(sep = "\n")
  srun <- glue_data(srun, "{key} {value}") %>% glue::collapse(sep = "\n")
  cm <- glue_data(cm, "{key} {value}") %>% glue::collapse(sep = "\n")
  script <- glue("#!/bin/bash\n\n# Arguments to SBATCH\n{sbatch}\n\n# Load modules\n{modules}\n\n# Run commands\n{srun}\n{cm}")
  con <- file("scr/test.sh")
  write_file(script, con)
  close(con)
}

parse_sbatch_yaml("sbatch.yaml")


