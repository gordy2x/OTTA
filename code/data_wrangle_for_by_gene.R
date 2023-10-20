

library(survival)
library(splines)
library(ggplot2)
library(gridExtra)
library(survival)
library(survminer)
library(ggpubr) 
library(dplyr)
library(tibble)



#directories 
datadir="data"
codedir<-"code"
ddatadir="data/derived"
if(!dir.exists(ddatadir)) dir.create(ddatadir)

#data
dat_fn<-"NanoII_n2539_final_post_QC_chips corrected_20200831"
subtypefn<-"Subtype GE n1462"

# truncate at 10 years
trunc.time = 10 * 365

#------------ Read in data and process-----------#

#data with patient characteristics from all patients
gene_data<-read.csv(file=paste0(datadir,"/",dat_fn,".csv")) %>% 
  select(OTTA.ID = OTTA_ID,
         site,
         finalstatus,
         timelastfu,
         stagenew,
         timeint_revised,
         refage_revised)

#for this analysis we want only the patients in the following file
subtype_dat<-read.csv(file=paste0(datadir,"/",subtypefn,".csv")) %>%  #only these paitents
  select(OTTA.ID=ID,mstype=Subtype,everything())

# join data, keeping only patients in subtype
cdat1 <- left_join(subtype_dat,gene_data) %>% 

# All observations are HGSC according to Path_Review_Diagnosis
# event, 0 alive, 1 dead
  mutate(event = finalstatus - 1) %>% 
  mutate_at(c("timelastfu", "timeint_revised", "refage_revised"), as.numeric) %>% # not sure why
  
# Right censor all data at 10 years
  mutate(event = ifelse(timelastfu > trunc.time, 0, event), #if alive till trunk_time then change to alive 
         timelastfu = ifelse(timelastfu > trunc.time, trunc.time, timelastfu), # truncate lastfu to trunc time
         
# A couple of weird things
# resign ".", to 0, not needed
#        timeint_revised = ifelse(timeint_revised == ".", 0, timeint_revised), 
# if first time is later than last then make first time 0?????       
         timeint_revised = ifelse(timeint_revised  > timelastfu, 0, timeint_revised),
# change 0 timelastfu to 1 ### seems arbitrary 
         timelastfu = ifelse(timelastfu <1, 1, timelastfu)) %>% 
  
# this is weird but needed to make the code match
# the 8's are actually missing
  mutate(stage = ifelse( !stagenew %in% c(1, 2), 8, stagenew )) %>% 

# remove observations with missing timelastfu
  filter(!is.na(timelastfu)) 
  

# try categorical and b-spline age variables
qtls = quantile( cdat1[, "refage_revised"], c( .3333, .6666 ), na.rm=TRUE )
age.bs = bs( cdat1[, "refage_revised"], knots = qtls ) 
mybasis = bs( cdat1[, "refage_revised"], knots = median(cdat1[, "refage_revised"], na.rm=TRUE) )

age.bs1 = as.data.frame( mybasis )
names(age.bs1) = paste( "bs", 1:ncol(age.bs1), sep="" )
cdat1 = cbind( cdat1, age.bs1 )
cdat1[, "age.cat.f" ] = cut( cdat1[, "refage_revised"], 5 )


#save data
save(cdat1, file=paste0(ddatadir, "/cleaned_gene_dat.Rdata"))
