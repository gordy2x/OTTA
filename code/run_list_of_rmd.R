
library(tidyverse)
for(irmd in c(1,3,7:12)) #
{
  rmarkdown::render(paste0(folder,"/",all_rmds$filenm[irmd]), 
                    params = list(final = TRUE))
}
