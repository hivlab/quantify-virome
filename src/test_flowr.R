library(flowr)

ex <- file.path(system.file(package = "flowr"), "pipelines")
flowmat <- as.flowmat(file.path(ex, "sleep_pipe.tsv"))
flowdef <- as.flowdef(file.path(ex, "sleep_pipe.def"))

fobj <- to_flow(x = flowmat, 
                def = flowdef,
                flowname = "example1", ## give it a name
                platform = "slurm")    ## override platform mentioned in fl

# plot_flow(fobj)
# plot_flow(flowdef)
submit_flow(fobj)
submit_flow(fobj, execute = TRUE)
# source("~/Projects/vs/src/test_flowr.R")
