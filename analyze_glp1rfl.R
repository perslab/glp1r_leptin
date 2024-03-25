## Load your packages, e.g. library(targets).
source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)
options(tidyverse.quiet = TRUE)

tar_option_set(
  format = "qs",
  error = "null",
  retrieval = "worker",
  storage = "worker"
)


knockout_glp1 <- 
  list(
    tar_target(
      agrpko,
      integrated_neuron_agrpko %>%
        map_ref(., dims=30, column_name = "internal_labs")
    ),
    tar_target(
      agrpko_scdist, 
      run_scdist(obj = agrpko, fixed.effects = "geno",  assay="SCT", column = "internal_labs",
                 random.effects = c("hash.mcl.ID", "hash_pool", "sex_predicted", "orig.ident"), 
                 d = 20, subset = celltypes)
    ),
    tar_target(
      celltypes_to_analyze,
      names(table(agrpko$newlabs))[table(glp1rko$newlabs)>200]
    ),
    tar_target(
      agrpko_pb,
      voom_wrapper(agrpko, mincells = 5, column = "internal_labs", 
                   pb_group = c("hash.mcl.ID"), filter_var = "Agrp")
    )
  )