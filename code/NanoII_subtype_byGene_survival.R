#########################################################
# Program Name: byGene_survival.R
# Purpose: Survival by gene and subtype
# controlling for age, stage, subtype, statified by site
# Based on overall_surv_fulldataset_MStypes_v3.R by Joshua Millstein, removng progression and race
# Need to run NanoII_overall_survival.R first to obtain predictions.

# Statistician: Gordana Popovic g.popovic@unsw.edu.au
# Date: 01/10/2020
# R version 4.0.2

# For each gene seperately using coxph: 
# Hazard ratio of gene controlling for subtype 
# Hazard ratio of interaction of gene and subtype
# Hazard ratio of gene subsetting by subtype 
# Create KM and violin plot by molecular subtype for each gene

library(survival)
library(splines)
library(ggplot2)
library(gridExtra)
library(survival)
library(survminer)
library(ggpubr) 
library(dplyr)
library(tibble)


## functions

## result function to clean up results of survival model

rslts.f = function(mod_main, mod_int = NULL) {
  res <- summary(mod_main)$conf.int %>% # start with summary
    as.data.frame() %>%  
    transmute(hr = `exp(coef)`,
              hr.L95 = `lower .95`,
              hr.U95 = `upper .95`) %>%  #rename and select variables
    rownames_to_column(var = "pred") %>%
    filter(pred == "gene") %>%  # filter just the gene row
    mutate(gene = gn,  # add other things
           p.gene = anova(mod_main)[ "gene", "Pr(>|Chi|)"],
           p.ms.gene.int = NA,
           n = mod_main$n) %>% 
    relocate(gene) %>% 
    select(-pred)
  
  if(is.null(mod_int)){
    res 
  } else{
    res %>% 
      mutate(p.GxMS.os = anova(mod_int)[ "mstype:gene", "Pr(>|Chi|)" ])
  }
}

## survival by gene and subtype
surv_single_gene = function(gn, dat, ids.exclude = NULL,
                            surv_var = c("timeint_revised", "timelastfu", "event"),
                            covar = c("bs1", "bs2", "bs3", "bs4", "stage.f", "site","mstype")){
  print(gn)
  
  # creates a dataset for analysis, 
  tmpdat <- dat %>%
    filter(!OTTA.ID %in% ids.exclude) %>% #remove id.exclude
    rename(gene = gn) %>% #select gene
    mutate_at(c("stage.f", "mstype", "site"), as.factor) %>%
    mutate(stage.f = relevel(stage.f, ref = "2"),# make factors and relevel
           gene = scale(gene)) %>% # scale gene
    select(gene, surv_var, covar)  # select relevant variables
  
  
  # Survival object with interval censoring
  my.surv.tmp = with(tmpdat,
                     Surv(
                       time = tmpdat[, "timeint_revised"],
                       time2 = tmpdat[, "timelastfu"],
                       event = tmpdat[, "event"]
                     ))
  
  ## Survival by gene controlling for all others
  mod_main = coxph( my.surv.tmp ~ bs1 + bs2 + bs3 + bs4 +
                      stage.f  + mstype + gene + strata(site),
                    data = tmpdat,
                    na.action = na.exclude
  )
  
  # add interaction between mstype and gene
  mod_int = update(mod_main, ~ . + mstype:gene)
  
  # summarise results of both models
  # Note: main effect of gene is from model without interaction
  rslts.f(mod_main, mod_int)
  
} 
  





kmplot = function( tmp.dat, bm.name  ){
  tmp.dat1 = as.data.frame(tmp.dat)
  names(tmp.dat1) = c("time", "event", "bm")
  bm = tmp.dat1[, "bm" ]
  ttls = quantile( bm, c( .2, .4, .6, .8 ), na.rm=TRUE )
  brks = c( min(bm,na.rm=TRUE) - 1, ttls, max(bm,na.rm=TRUE) + 1 )
  tgene = cut( bm, breaks=brks )
  tgn = ifelse( is.element( tgene, levels(tgene)[1] ), "Q1", tgene )
  tgn = ifelse( is.element( tgene, levels(tgene)[2] ), "Q2", tgn )
  tgn = ifelse( is.element( tgene, levels(tgene)[3] ), "Q3", tgn )
  tgn = ifelse( is.element( tgene, levels(tgene)[4] ), "Q4", tgn )
  tgn = ifelse( is.element( tgene, levels(tgene)[5] ), "Q5", tgn )
  strata = factor(tgn, levels=c("Q1", "Q2", "Q3", "Q4", "Q5"))	
  tmp.dat2 = cbind(tmp.dat1, strata)
  names(tmp.dat2) = c("time1", "event1", "bm1", "strata1")
  tmp.dat2[,"time1"] = tmp.dat2[,"time1"] / 365
  myfit = survfit( Surv(time1, event1 ) ~ strata1, data=tmp.dat2 ) 
  myplot = ggsurvplot(myfit, data=tmp.dat2, submain=bm.name, 
                      pval = TRUE, conf.int = FALSE, # logrank p-value
                      risk.table = FALSE, 
                      risk.table.y.text.col = TRUE,
                      surv.median.line = "hv",  # add the median survival pointer.
                      xlab="Time (years)",
                      linetype = "strata",
                      ggtheme = theme_bw(), legend="bottom",
                      legend.labs = c("Q1", "Q2", "Q3", "Q4", "Q5"),
                      palette = c( "#99CCCC", "#0099CC", "#6666FF", "#0033CC", "#663366" ) ) 
  return( myplot )
} # End kmplot

vplot = function( pdat, gn.nm ){
  names(pdat) = c("mstype", "gene")
  pdat[,"gene"] = scale(pdat[,"gene"])
  myplot = ggplot( pdat, aes( mstype, gene ) ) + 
    xlab("molecular subtype") + ylab("normalized expression") + theme_bw() +
    geom_violin(aes(fill=mstype)) + 
    geom_point(size=.05) +
    geom_jitter(width = 0.05, height=0.0, alpha=.2) + 
    labs(title=gn.nm) + guides(fill=FALSE) + 
    geom_boxplot(width=0.05, notch=TRUE, outlier.shape=NA) +
    stat_compare_means()
  return(myplot)
}


## compute bh on observed p-values
bh = function(pvals){             
  m = length(pvals)
  po = pvals[order(pvals)]
  i = 1:m
  new.order = i[order(pvals)]
  p.adj = po * m /  i
  for(i in 1:m){
    p.adj[ 1:i ] = ifelse( p.adj[ 1:i ] > p.adj[ i ], p.adj[ i ], p.adj[ 1:i ] )
    if( p.adj[ i ] > .9 ) break
  }
  ind = order( new.order )  # reorder back to original
  p.adj = p.adj[ ind ]
  return(pmin(p.adj,1))
}



#directories 
resdir="results/Subtype"
if(!dir.exists(resdir)) dir.create(resdir)
datadir="data"
codedir<-"code"
ddatadir="data/derived"


#data
nanoIIfn<-"with_predict_NanoII_n2539_final_post_QC_chips corrected_20200831"
subtypefn<-"Subtype GE n1462"


#------------ Read in data and process-----------#

#data with patient characteristics from all patients
pred_dat<-read.csv(file=paste0(ddatadir,"/",nanoIIfn,".csv")) %>% 
  select(OTTA.ID=OTTA_ID,
         site,
         # predict,
         finalstatus,
         timelastfu,
         stagenew,
         timeint_revised,
         refage_revised,
         event)

#for this analysis we want only the patients in the following file
subtype_dat<-read.csv(file=paste0(datadir,"/",subtypefn,".csv")) %>%  #only these paitents
  select(OTTA.ID=ID,mstype=Subtype,everything())

# join data, keeping only patients in subtype
cdat<-left_join(subtype_dat,pred_dat)

# All observations are HGSC according to Path_Review_Diagnosis
# event, 0 alive, 1 dead
cdat[, "event"] = cdat[, "finalstatus"] - 1

# Right censor all data at 10 years
cdat[, "timelastfu" ] = as.numeric(cdat[, "timelastfu" ])

trunc.time = 10 * 365
cdat[, "event" ] = ifelse( cdat[, "timelastfu" ] > trunc.time, 0, cdat[, "event" ] )
cdat[, "timelastfu" ] = as.numeric(cdat[, "timelastfu" ])
cdat[, "timelastfu" ] = ifelse( cdat[, "timelastfu" ] > trunc.time, trunc.time, cdat[, "timelastfu" ] )

# use "stagenew" for stage
# use "refage_revised" for age
# eliminate small groups
# cdat[, "race" ] = ifelse( is.element( cdat[, "race" ] , c(".")), "8", cdat[, "race" ] )
cdat[, "stage" ] = ifelse( !is.element( cdat[, "stagenew" ] , c("1", "2")), "8", cdat[, "stagenew" ] )

# make stage and grade factor variables
cdat[, "stage.f" ] = as.factor( cdat[, "stage" ] ) # missing for 95 individuals, BRO site
# cdat[, "race.f" ] = as.factor( cdat[, "race" ] )
cdat[, "stage.f" ] = relevel( cdat[, "stage.f" ], "2" )

cdat[, "timeint_revised" ] = ifelse( is.element( cdat[, "timeint_revised" ] , "."), "0", cdat[, "timeint_revised" ] ) # resign ".", to 0
cdat[, "timeint_revised" ] = as.numeric( cdat[, "timeint_revised" ] )
cdat[, "timeint_revised" ] = ifelse( cdat[, "timeint_revised" ] >= cdat[, "timelastfu" ] , 0, cdat[, "timeint_revised" ] )

cdat[, "refage_revised" ] = as.numeric( cdat[, "refage_revised" ] )

# remove observations with missing timelastfu
cdat1 = cdat[ !is.na( cdat[, "timelastfu" ] ), ] 

# change 0 timelastfu to 1
cdat1[, "timelastfu" ] = ifelse( cdat1[, "timelastfu" ] < 1, 1, cdat1[, "timelastfu" ] )

# try categorical and b-spline age variables
qtls = quantile( cdat1[, "refage_revised"], c( .3333, .6666 ), na.rm=TRUE )
age.bs = bs( cdat1[, "refage_revised"], knots = qtls ) 
mybasis = bs( cdat1[, "refage_revised"], knots = median(cdat1[, "refage_revised"], na.rm=TRUE) )

age.bs1 = as.data.frame( mybasis )
names(age.bs1) = paste( "bs", 1:ncol(age.bs1), sep="" )
cdat1 = cbind( cdat1, age.bs1 )
cdat1[, "age.cat.f" ] = cut( cdat1[, "refage_revised"], 5 )

#save data
save(cdat1, file=paste0(ddatadir, "cdat1.Rdata"))


#---------------------Survival by subtype and gene----------------------#

genes_all=colnames(cdat1)[3:337]
rslts_all = mysurv_sub( genes_all )

#adjust p-values
p.nms = c("p.os", "p.GxMS.os")
q.nms = c("q.os", "q.GxMS.os")
for(j in 1:length(p.nms)) rslts_all[ , q.nms[j]] = bh( rslts_all[ , p.nms[j]] )


#save 
write.csv(rslts_all, file=paste0(resdir, "/subtype_by_gene.csv"))


#------------ Create KM and violin molecular subtype plots for each gene-----------#

gene.names=genes_all
mssubtypes = as.character(unique(cdat1$mstype))


outc="os"   #no progression 
  for(mst in mssubtypes){
    mstdir<-paste0(resdir,"/",mst) #directory
    if(!dir.exists(mstdir)) dir.create(mstdir)
    msdat = cdat1[is.element(cdat1$mstype, mst), ]
      msdat[,"time.plot"] = msdat[, "timelastfu" ] - msdat[, "timeint_revised" ]
      msdat[,"event.plot"] = msdat[, "event"]

    for(gn in gene.names){
      msdat[,"gene"] = msdat[, gn]
      p1 = kmplot( msdat[,c("time.plot", "event.plot", "gene")], gn )
      outf = paste0(resdir,"/", mst, "/KM_", mst, "_", gn,".pdf")
      pdf(outf, width = 6, height=6)
        grid.arrange( p1$plot, ncol=1 )
      dev.off()
      # print(paste(c(outc,mst,gn)))
    }
  }

vioresdir<-paste0(resdir,"/violin")
if(!dir.exists(vioresdir)) dir.create(vioresdir)

for(gn in gene.names){
  msdat = cdat1[, c("mstype", gn)]
  p1 = vplot( msdat, gn ) 
  outf =paste0(resdir,"/violin/","violin_",gn,".pdf", sep="")
  pdf(outf, width = 6, height=6)
    grid.arrange( p1, ncol=1 )
  dev.off()
  # print(gn)
}

