
library(tidyverse)
load(file=paste0(ddatadir, "cleaned_gene_dat.Rdata"))


for(irmd in c(1,3,7:12)) #
{
  rmarkdown::render("data/by_gene_survival.Rmd",
                    output_file = paste0("results/gn")
                    params = list(final = TRUE))
}
