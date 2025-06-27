#packages installed
install.packages("ggplot2")
install.packages("foreign")
install.packages("graphics") 
install.packages("tidyr")
install.packages("emmeans")
install.packages("lmerTest")
install.packages("mixtools")
install.packages("reshape2")

# used libraries
library(foreign)
require(graphics)
library(tidyr)
require(emmeans)
library(lmerTest)
require(mixtools)
library(reshape2)

# clear environment
rm(list = ls())

# set working directory
setwd("/Volumes/education/FHML_MHeNs/Intern_Naia_Barney_Machado/Dataset/results/")

# set data directory
data_dir="/Volumes/education/FHML_MHeNs/Intern_Naia_Barney_Machado/Dataset/data/"  #VB van mijn pad: '/Users/bettytijms/Documents/AD_11-12/projecten/MemoProteomicsAD/data/EMIF/'

# Read and load the protein data in 20210407_EMIF_tryp_Proteins_patient_id_test
my_data2=read.delim(paste(data_dir,'20210407_EMIF_tryp_Proteins_patient_id_test.txt',sep=''))

# split the variable description into multiple columns 
bla=my_data2$Description
bla=gsub(' OS=',';OS=',bla)
bla=gsub('GN=',';GN=',bla)
bla=gsub(' PE=',';PE=',bla)
bla=gsub('SV=',';SV=',bla)
bla2=strsplit(bla,';')
my_data2=cbind.data.frame(full.name=sapply(bla2,'[[',1),GENE=sapply(bla2,'[[',3),my_data2)

my_data2$GENE=as.vector(my_data2$GENE)
my_data2$Accession=as.vector(my_data2$Accession)
my_data2$GENE=gsub('GN=','',my_data2$GENE)

# if no gene is known enter full name
my_data2$GENE[my_data2$GENE=='PE=1 ']=my_data2$Accession[my_data2$GENE=='PE=1 ']
my_data2$GENE[my_data2$GENE=='PE=5 ']=my_data2$Accession[my_data2$GENE=='PE=5 ']
my_data2$GENE[my_data2$GENE=='PE=2 ']=my_data2$Accession[my_data2$GENE=='PE=2 ']

# when there is a double gene name but own uniprot -> do the combination
my_data2$GENE[which(duplicated(my_data2$GENE))]=paste(my_data2$GENE[which(duplicated(my_data2$GENE))],my_data2$Accession[which(duplicated(my_data2$GENE))],sep='.')
which(duplicated(my_data2$GENE))

# Column MAA.492.1 => should be MAA.492
# Column MAA.492 ==> should be MAA.502
names(my_data2)[which(names(my_data2)=='MAA.492')]='MAA.502'
names(my_data2)[which(names(my_data2)=='MAA.492.1')]='MAA.492' 

# check if any individuals don't have any measurements
na_c=colSums(!is.na(my_data2[,16:697]))
any(na_c==0)

load(paste(data_dir,'prot_mat.RData',sep=''))

#### define contrasts ####

# SNAP.contrast
prot_mat$SNAP.contrast=NA 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='0')]='CN A-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='0')]='CN A+ Hyp-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='1')]='CN A+ Hyp+' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='MCI' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='0')]='MCI A+ Hyp-'
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='MCI' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='1')]='MCI A+ Hyp+' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='AD' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='0')]='AD A+ Hyp-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='AD' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$Hypertension=='1')]='AD A+ Hyp+' 

prot_mat$SNAP.contrast=factor(prot_mat$SNAP.contrast,levels=c("CN A-", "CN A+ Hyp-", "CN A+ Hyp+", "MCI A+ Hyp-", "MCI A+ Hyp+", "AD A+ Hyp-", "AD A+ Hyp+"))
table(prot_mat$SNAP.contrast)

#What I get:
#CN A-    CN A+ Hyp-    CN A+ Hyp+  MCI A+ Hyp-   MCI A+ Hyp+   AD A+ Hyp-    AD A+ Hyp+ 
#  171            24            13          48            31             8             9

prot_mat$hoofd.contrast=NA


####### Protein comparison ######
prot_ind=2:3104
# check if all values ​​are positive (=TRUE) because it could be a problem for log transformation
all(prot_mat[,prot_ind]>0,na.rm=T)

# take the log
prot_mat[,prot_ind]=apply(prot_mat[,prot_ind],2,log)
all_prot_mat=prot_mat

# Standardization (z-score) in reference to the control group
for(i in seq_along(prot_ind)){
  tmean=mean(prot_mat[prot_mat$SNAP.contrast=='CN A-', ] [[prot_ind[i]]],na.rm=T)
  tsd=sd(prot_mat[prot_mat$SNAP.contrast=='CN A-', ] [[prot_ind[i]]],na.rm=T)
  prot_mat[[prot_ind[i]]]=(prot_mat[[prot_ind[i]]]-tmean)/tsd  
}

#check if it worked in the first 10 colums
apply(prot_mat[ ,prot_ind[1:10]],2,function(x) mean(x,na.rm=T))
#A1BG          A2M        A2ML1       A4GALT       A6NIZ1       A6NNZ2 
#0.008939461 -0.008863143 -2.168525392  0.054839510 -0.642465490 -0.087749670 
#A8MTW9        AARS1         ABAT        ABCA2 
#0.194628177  0.197449254  0.159615613  0.015688602 

by(prot_mat[[prot_ind[1]]],prot_mat$SNAP.contrast,function(x) mean(x,na.rm=T))
   
prot_mat=cbind.data.frame(prot_mat,Central.NFL=rep(NA, nrow(prot_mat)),Central.NRGN=rep(NA, nrow(prot_mat)),Central.YKL40=rep(NA, nrow(prot_mat)),Central.AB38=rep(NA, nrow(prot_mat)),Central.AB40=rep(NA, nrow(prot_mat)),Central.AB42=rep(NA, nrow(prot_mat)))

to_z=c("Central_CSF_NFL","Central_CSF_Neurogranin", "Central_CSF_YKL40", "Central_CSF_AB38", "Central_CSF_AB40", "Central_CSF_AB42")

z_ed=3160:3165
names(prot_mat)[z_ed]

for(i in seq_along(z_ed)){
  
  tmean= mean(prot_mat[prot_mat$SNAP.contrast=='CN A-',to_z[i]],na.rm=T)
  tsd=sd(prot_mat[prot_mat$SNAP.contrast=='CN A-',to_z[i]],na.rm=T)
  prot_mat[,z_ed[i]]= (prot_mat[,to_z[i]]-tmean)/tsd
}
   
# Results matrix (column names)
res_mat = data.frame(matrix(NA, nrow = length(prot_ind), ncol =66))
names(res_mat)=c('Protein',
                 'Uniprot',
                 'number of observations',
                 'CN A- count',
                 'CN A+ Hyp- count',
                 'CN A+ Hyp+ count',
                 'MCI A+ Hyp- count',
                 'MCI A+ Hyp+ count',
                 'AD A+ Hyp- count',
                 'AD A+ Hyp+ count',
                 'CN A- (controls) mean (sd)',
                 'CN A- (controls) emm (se)',
                 'CN A+ Hyp- mean (sd)',
                 'CN A+ Hyp- emm (se)', 
                 'CN A+ Hyp+ mean (sd)',
                 'CN A+ Hyp+ emm(se)',
                 'MCI A+ Hyp- mean (sd)',
                 'MCI A+ Hyp- emm (se)',
                 'MCI A+ Hyp+ mean (sd)',
                 'MCI A+ Hyp+ emm (se)',
                 'AD A+ Hyp- mean (sd)',
                 'AD A+ Hyp- emm (se)',
                 'AD A+ Hyp+ mean (sd)',
                 'AD A+ Hyp+ emm (se)',
                 
                 'diff CN A- (controls) vs CN A+ Hyp-',
                 'p CN A- (controls) vs CN A+ Hyp-',
                 'diff CN A- (controls) vs CN A+ Hyp+',
                 'p CN A- (controls) vs NC A+ Hyp+',
                 'diff CN A- (controls) vs MCI A+ Hyp-',
                 'p CN A- (controls) vs MCI A+ Hyp-',
                 'diff CN A- (controls) vs MCI A+ Hyp+',
                 'p CN A- (controls) vs MCI A+ Hyp+',
                 'diff CN A- (controls) vs AD A+ Hyp-',
                 'p CN A- (controls) vs AD A+ Hyp-',
                 'diff CN A- (controls) vs AD A+ Hyp+',
                 'p CN A- (controls) vs AD A+ Hyp+'
                 
                 'diff CN A+ Hyp- vs CN A+ Hyp+',
                 'p CN A+ Hyp- vs CN A+ Hyp+',
                 'diff CN A+ Hyp- vs MCI A+ Hyp-',
                 'p CN A+ Hyp- vs MCI A+ Hyp-',
                 'diff CN A+ Hyp- vs MCI A+ Hyp+',
                 'p CN A+ Hyp- vs MCI A+ Hyp+',
                 'diff CN A+ Hyp- vs AD A+ Hyp-',
                 'p CN A+ Hyp- vs AD A+ Hyp-',
                 'diff CN A+ Hyp- vs AD A+ Hyp+',
                 'p CN A+ Hyp- vs AD A+ Hyp+',
                                  
                 'diff CN A+ Hyp+ vs MCI A+ Hyp-',
                 'p CN A+ Hyp+ vs MCI A+ Hyp-',
                 'diff CN A+ Hyp+ vs MCI A+ Hyp+',
                 'p CN A+ Hyp+ vs MCI A+ Hyp+',
                 'diff CN A+ Hyp+ vs AD A+ Hyp-',
                 'p CN A+ Hyp+ vs AD A+ Hyp-',
                 'diff CN A+ Hyp+ vs AD A+ Hyp+',
                 'p CN A+ Hyp+ vs AD A+ Hyp+',
                 
                 'diff MCI A+ Hyp- vs MCI A+ Hyp+',
                 'p MCI A+ Hyp- vs MCI A+ Hyp+',
                 'diff MCI A+ Hyp- vs AD A+ Hyp-',
                 'p MCI A+ Hyp- vs MCI A+ Hyp-',
                 'diff MCI A+ Hyp- vs AD A+ Hyp+',
                 'p MCI A+ Hyp- vs MCI A+ Hyp+',
                 
                 'diff MCI A+ Hyp+ vs AD A+ Hyp-',
                 'p MCI A+ Hyp- vs AD A+ Hyp-',
                 'diff MCI A+ Hyp- vs AD A+ Hyp+',
                 'p MCI A+ Hyp- vs MCI A+ Hyp+',
                 
                 'diff AD A+ Hyp- vs AD A+ Hyp+',
                 'p MCI A+ Hyp- vs MCI A+ Hyp+')

rownames(res_mat)=names(prot_mat)[prot_ind]

# center for AGE                 
prot_mat$AGE.centered=scale(prot_mat$Age,center=T,scale=F)
# Check if gender and apoe are a factor
prot_mat$Gender=factor(prot_mat$Gender)
prot_mat$APOEdich=factor(prot_mat$APOEdich)

# take a subset of prot_mat with only the groups of interest
sub_prot_mat=prot_mat[which(!is.na(prot_mat$SNAP.contrast)), ]

# Determine the differences between the groups ONLY if there are >9 observations per group.
for(i in 1:length(prot_ind)){
  
  #wat is de uniprot?
  if(i >6){
    res_mat[i,2]=as.vector(my_data2$Accession[which(my_data2$GENE==names(sub_prot_mat)[prot_ind[i]])])
  }
  
  res_mat[i,3]=sum(!is.na(sub_prot_mat[,prot_ind[i]]))
  
  # only continue if there are at least 10 observations!
  if(res_mat[i,3]>9){
    
    # Count number of observations
    group_counts = tapply(!is.na(sub_prot_mat[,prot_ind[i]]), sub_prot_mat$SNAP.contrast, sum)
    
    # Store counts in columns 4 to 10, matching the SNAP.contrast levels
    res_mat[i, 4:10] = group_counts[levels(sub_prot_mat$SNAP.contrast)]
    
    tprot_mat=sub_prot_mat
    
    tprot=tprot_mat[,prot_ind[i]] 
    tres=by(tprot,tprot_mat$SNAP.contrast,function(x) mean(x, na.rm=T))
    tres2=by(tprot,tprot_mat$SNAP.contrast,function(x) sd(x, na.rm=T))
    
    # Note columns adjust depending on contrast you run
    res_mat[i,c(11,13,15,17,19, 21, 23)]=paste(round(tres,2), ' (', round(tres2,2),')',sep='')
    
    # ANCOVA with age,sex, APOE correction
    if(all(table(tprot_mat$SNAP.contrast[!is.na(tprot)])>2)){
      
      # Compare groups with linear model; covariates age + gender+ apoe
      tres=lm(tprot~AGE.centered+Gender+APOEdich+SNAP.contrast,data=tprot_mat)
      # Extract the marginal means
      temmeans=test(emmeans(tres,pairwise~SNAP.contrast),adjust='none')
      
      # add the emm (se) to the table
      res_mat[i, c(12,14,16,18,20, 22, 24)]=paste(round(temmeans$emmeans$emmean,2),' (',
                                               round(temmeans$emmeans$SE,2),
                                               ')',sep='')
      
      # add the differences & the pvalues between the groups
      res_mat[i, c(25,27,29,31,33,35,37,39,41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65)]=round(temmeans$contrasts$estimate,2)
      res_mat[i, c(26,28,30,32,34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66)]=round(temmeans$contrasts$p.value,5)
      
    }
  }
  
}


# save the results
write.table(res_mat,'resultsproteomics2_hypertension.txt',sep='\t')

