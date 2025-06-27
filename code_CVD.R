

#################################################################################
#
#     Doel : proteomic verschillen voor bepaald contrast bepalen
#
#
#################################################################################

###################################################
##### DATA PREPARATION
###################################################


##install packages: hoeft maar 1 keer (na updates)
install.packages("ggplot2")
install.packages("foreign")
install.packages("graphics") #not available for 3.5.1 version?
install.packages("tidyr")
install.packages("emmeans")
install.packages("lmerTest")
install.packages("mixtools")
install.packages("reshape2")

## load packages: Hoeft alleen als je script opstart
# require(ggplot2) # als je een package nog niet hebt, kun je eenmalig de volgende regel runnen - zie boven -
# install.packages("ggplot2")
library(foreign)
require(graphics)
library(tidyr)
require(emmeans)
library(lmerTest)
require(mixtools)
library(reshape2)

## clear environment
rm(list = ls())

# Pas het pad aan naar de directory (=map) waar je de resultaten wil opslaan
setwd("/Volumes/education/FHML_MHeNs/Intern_Naia_Barney_Machado/Dataset/results/")

# Geef de directory waar je de data vandaan haalt
data_dir="/Volumes/education/FHML_MHeNs/Intern_Naia_Barney_Machado/Dataset/data/"  #VB van mijn pad: '/Users/bettytijms/Documents/AD_11-12/projecten/MemoProteomicsAD/data/EMIF/'

# lees de eiwitten data in 
#my_data=read.delim('/Users/stephanievos/Desktop/MCISNAPproteomics/data/20180329_EMIF_tryp_Proteins_patient_id_test.txt')
##### laad de data
my_data2=read.delim(paste(data_dir,'20210407_EMIF_tryp_Proteins_patient_id_test.txt',sep=''))

# split de variabele description in meerdere kolommen: Zo kunnen we gen naam gebruiken ipv uniprot (is makkelijker te lezen)
bla=my_data2$Description
bla=gsub(' OS=',';OS=',bla)
bla=gsub('GN=',';GN=',bla)
bla=gsub(' PE=',';PE=',bla)
bla=gsub('SV=',';SV=',bla)
bla2=strsplit(bla,';')

my_data2=cbind.data.frame(full.name=sapply(bla2,'[[',1),GENE=sapply(bla2,'[[',3),my_data2)
# nu GN= weg
my_data2$GENE=as.vector(my_data2$GENE)
my_data2$Accession=as.vector(my_data2$Accession)
my_data2$GENE=gsub('GN=','',my_data2$GENE)
# NB soms is er geen gen bekend -> doe dan de vollegide naam erin
my_data2$GENE[my_data2$GENE=='PE=1 ']=my_data2$Accession[my_data2$GENE=='PE=1 ']
my_data2$GENE[my_data2$GENE=='PE=5 ']=my_data2$Accession[my_data2$GENE=='PE=5 ']
# Er is ook een PE=2 --> welke is dit?
my_data2$GENE[my_data2$GENE=='PE=2 ']=my_data2$Accession[my_data2$GENE=='PE=2 ']

# NB soms dubbele gen naam maar eigen uniprot -> doe de combi
my_data2$GENE[which(duplicated(my_data2$GENE))]=paste(my_data2$GENE[which(duplicated(my_data2$GENE))],my_data2$Accession[which(duplicated(my_data2$GENE))],sep='.')
which(duplicated(my_data2$GENE))

# Column MAA.492.1 => moet zijn MAA.492
# Column MAA.492 ==> moet zijn MAA.502
names(my_data2)[which(names(my_data2)=='MAA.492')]='MAA.502'
names(my_data2)[which(names(my_data2)=='MAA.492.1')]='MAA.492' 

# zijn er mensen bij die geen enkele meting hebben?
na_c=colSums(!is.na(my_data2[,16:697]))
any(na_c==0) # nee, er is blijkbaar altijd een getal

load(paste(data_dir,'prot_mat.RData',sep=''))

#### definieer de contrasten ####

# SNAP.contrast
prot_mat$SNAP.contrast=NA # intitialiseer eerst
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='0')]='CN A-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='0')]='CN A+ CVD-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='CN' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='1')]='CN A+ CVD+' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='MCI' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='0')]='MCI A+ CVD-'
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='MCI' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='1')]='MCI A+ CVD+' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='AD' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='0')]='AD A+ CVD-' 
prot_mat$SNAP.contrast[which(prot_mat$Diagnosis.Aurore=='AD' & prot_mat$Local_AB42_Abnormal.betty=='1' & prot_mat$CardiovascularDis=='1')]='AD A+ CVD+' 

prot_mat$SNAP.contrast=factor(prot_mat$SNAP.contrast,levels=c("CN A-", "CN A+ CVD-", "CN A+ CVD+", "MCI A+ CVD-", "MCI A+ CVD+", "AD A+ CVD-", "AD A+ CVD+"))
table(prot_mat$SNAP.contrast)

#What I get with CVD:
#CN A-    CN A+ CVD-    CN A+ CVD+  MCI A+ CVD-   MCI A+ CVD+   AD A+ CVD-    AD A+ CVD+ 
#  171            11            3           45            10            20            5 

#What I get with hypertension:
#CN A-    CN A+ Hyp-    CN A+ Hyp+  MCI A+ Hyp-   MCI A+ Hyp+   AD A+ Hyp-    AD A+ Hyp+ 
#  171            24            13          48            31             8             9

#What I get with WMH:
#CN A-    CN A+ WMH-    CN A+ WMH+  MCI A+ WMH-   MCI A+ WMH+   AD A+ WMH-    AD A+ WMH+ 
#  171            39             6          48            23             0             0 

#obesity
#CN A-  CN A+ Obe-  CN A+ Obe+ MCI A+ Obe- MCI A+ Obe+  AD A+ Obe-  AD A+ Obe+ 
#  171          17           1          75           4          15           3 

prot_mat$hoofd.contrast=NA


####### Hier gaan we nu de eiwitten vergelijken ######
prot_ind=2:3104 # dit zijn de kolommen met eiwitten
# check of all waarden positief zijn (=TRUE) - anders probleem zijn voor log transformatie
all(prot_mat[,prot_ind]>0,na.rm=T)

# neem de log (de ln = natural logarithm, net als Spellman)
prot_mat[,prot_ind]=apply(prot_mat[,prot_ind],2,log)
all_prot_mat=prot_mat

# schaal alles op reference groep. Standardization (z-score)
for(i in seq_along(prot_ind)){
  tmean=mean(prot_mat[prot_mat$SNAP.contrast=='CN A-', ] [[prot_ind[i]]],na.rm=T)
  tsd=sd(prot_mat[prot_mat$SNAP.contrast=='CN A-', ] [[prot_ind[i]]],na.rm=T)
  
  prot_mat[[prot_ind[i]]]=(prot_mat[[prot_ind[i]]]-tmean)/tsd  
}


#why only the first 10 columns?
apply(prot_mat[ ,prot_ind[1:10]],2,function(x) mean(x,na.rm=T))


by(prot_mat[[prot_ind[1]]],prot_mat$SNAP.contrast,function(x) mean(x,na.rm=T))
#I get negative means from line 135, is this because the log is applied? and why is it applied?


#### Vergelijk eiwit concentraties voor contrasten  #####
# Ik voeg hier ook Z getransformeerde extra gemeten eiwitten aan prot_mat toe zodat je die ook meteen kunt vergelijken

prot_mat=cbind.data.frame(prot_mat,Central.NFL=rep(NA, nrow(prot_mat)),Central.NRGN=rep(NA, nrow(prot_mat)),Central.YKL40=rep(NA, nrow(prot_mat)),Central.AB38=rep(NA, nrow(prot_mat)),Central.AB40=rep(NA, nrow(prot_mat)),Central.AB42=rep(NA, nrow(prot_mat)))

to_z=c("Central_CSF_NFL","Central_CSF_Neurogranin", "Central_CSF_YKL40", "Central_CSF_AB38", "Central_CSF_AB40", "Central_CSF_AB42")

z_ed=3160:3165
# check of is z_ed naar juist kolommen verwijst
names(prot_mat)[z_ed]

# loop door de kolommen & normaliseer tov controle groep
for(i in seq_along(z_ed)){
  
  tmean= mean(prot_mat[prot_mat$SNAP.contrast=='CN A-',to_z[i]],na.rm=T)
  tsd=sd(prot_mat[prot_mat$SNAP.contrast=='CN A-',to_z[i]],na.rm=T)
  prot_mat[,z_ed[i]]= (prot_mat[,to_z[i]]-tmean)/tsd
  
}

# paste these proteins to the prot_ind
prot_ind=c(z_ed, prot_ind)

# initialiseer een resultaten matrix --> ncol aanpassen naar contrast wat je wil
# idem met de names; dit is alles voor pj
# Ik sla ruwe (=mean (sd) & estimated marginal means (=emm (se)) waarden op om invloed van covariaten te zien
res_mat = data.frame(matrix(NA, nrow = length(prot_ind), ncol =66))

names(res_mat)=c('Protein',
                 'Uniprot',
                 'number of observations',
                 'CN A- count',
                 'CN A+ CVD- count',
                 'CN A+ CVD+ count',
                 'MCI A+ CVD- count',
                 'MCI A+ CVD+ count',
                 'AD A+ CVD- count',
                 'AD A+ Hyp+ count',
                 'CN A- (controls) mean (sd)',
                 'CN A- (controls) emm (se)',
                 'CN A+ CVD- mean (sd)',
                 'CN A+ CVD- emm (se)', 
                 'CN A+ CVD+ mean (sd)',
                 'CN A+ CVD+ emm(se)',
                 'MCI A+ CVD- mean (sd)',
                 'MCI A+ CVD- emm (se)',
                 'MCI A+ CVD+ mean (sd)',
                 'MCI A+ CVD+ emm (se)',
                 'AD A+ CVD- mean (sd)',
                 'AD A+ CVD- emm (se)',
                 'AD A+ CVD+ mean (sd)',
                 'AD A+ CVD+ emm (se)',
                 
                 'diff CN A- (controls) vs CN A+ CVD-',
                 'p CN A- (controls) vs CN A+ CVD-',
                 'diff CN A- (controls) vs CN A+ CVD+',
                 'p CN A- (controls) vs CN A+ CVD+',
                 'diff CN A- (controls) vs MCI A+ CVD-',
                 'p CN A- (controls) vs MCI A+ CVD-',
                 'diff CN A- (controls) vs MCI A+ CVD+',
                 'p CN A- (controls) vs MCI A+ CVD+',
                 'diff CN A- (controls) vs AD A+ CVD-',
                 'p CN A- (controls) vs AD A+ CVD-',
                 'diff CN A- (controls) vs AD A+ CVD+',
                 'p CN A- (controls) vs AD A+ CVD+',
                 
                 
                 'diff CN A+ CVD- vs CN A+ CVD+',
                 'p CN A+ CVD- vs CN A+ CVD+',
                 'diff CN A+ CVD- vs MCI A+ CVD-',
                 'p CN A+ CVD- vs MCI A+ CVD-',
                 'diff CN A+ CVD- vs MCI A+ CVD+',
                 'p CN A+ CVD- vs MCI A+ CVD+',
                 'diff CN A+ CVD- vs AD A+ CVD-',
                 'p CN A+ CVD- vs AD A+ CVD-',
                 'diff CN A+ CVD- vs AD A+ CVD+',
                 'p CN A+ CVD- vs AD A+ CVD+',
                 
                 
                 'diff CN A+ CVD+ vs MCI A+ CVD-',
                 'p CN A+ CVD+ vs MCI A+ CVD-',
                 'diff CN A+ CVD+ vs MCI A+ CVD+',
                 'p CN A+ CVD+ vs MCI A+ CVD+',
                 'diff CN A+ CVD+ vs AD A+ CVD-',
                 'p CN A+ CVD+ vs AD A+ CVD-',
                 'diff CN A+ CVD+ vs AD A+ CVD+',
                 'p CN A+ CVD+ vs AD A+ CVD+',
                 
                 'diff MCI A+ CVD- vs MCI A+ CVD+',
                 'p MCI A+ CVD- vs MCI A+ CVD+',
                 'diff MCI A+ CVD- vs AD A+ CVD-',
                 'p MCI A+ CVD- vs MCI A+ CVD-',
                 'diff MCI A+ CVD- vs AD A+ CVD+',
                 'p MCI A+ CVD- vs MCI A+ CVD+',
                 
                 'diff MCI A+ CVD+ vs AD A+ CVD-',
                 'p MCI A+ CVD- vs AD A+ CVD-',
                 'diff MCI A+ CVD- vs AD A+ CVD+',
                 'p MCI A+ CVD- vs AD A+ CVD+',
                 
                 'diff AD A+ CVD- vs AD A+ CVD+',
                 'p AD A+ CVD- vs AD A+ CVD+')


rownames(res_mat)=names(prot_mat)[prot_ind]

# center AGE                 
prot_mat$AGE.centered=scale(prot_mat$Age,center=T,scale=F)
# Zorg dat gender een factor is
prot_mat$Gender=factor(prot_mat$Gender)
prot_mat$APOEdich=factor(prot_mat$APOEdich)

# take a subset of prot_mat with only the groups of interest
sub_prot_mat=prot_mat[which(!is.na(prot_mat$SNAP.contrast)), ]

# walk through the proteins and determine the differences between the groups: NB only if there are >9 observations per group.
for(i in 1:length(prot_ind)){
  
  #wat is de uniprot?
  if (i>6){
    
    res_mat[i,2]=as.vector(my_data2$Accession[which(my_data2$GENE==names(sub_prot_mat)[prot_ind[i]])])
  }
  
  
  #how many people have an observation in the entire group?
  res_mat[i,3]=sum(!is.na(sub_prot_mat[,prot_ind[i]]))
  
  
  
  
  # only continue if there are at least 10 observations!
  if(res_mat[i,3]>9){
    # compare this protein!
    
    # Count number of observations per SNAP.contrast group for this protein
    group_counts = tapply(!is.na(sub_prot_mat[,prot_ind[i]]), sub_prot_mat$SNAP.contrast, sum)
    
    # Store counts in columns 4 to 10, matching the SNAP.contrast levels
    res_mat[i, 4:10] = group_counts[levels(sub_prot_mat$SNAP.contrast)]
    
    tprot_mat=sub_prot_mat
    
    tprot=tprot_mat[,prot_ind[i]] # take the protein
    # first the raw values for everything and everyone:
    # 1 so first per ab_tau_prof category
    tres=by(tprot,tprot_mat$SNAP.contrast,function(x) mean(x, na.rm=T))
    tres2=by(tprot,tprot_mat$SNAP.contrast,function(x) sd(x, na.rm=T))
    
    # Note columns adjust depending on contrast you run
    res_mat[i,c(11,13,15,17,19, 21, 23)]=paste(round(tres,2), ' (', round(tres2,2),')',sep='')
    
    # Now the statistics: so age/sex correction anyway
    # only possible if everyone has observations
    # So adjust to variable with the contrast you test
    if(all(table(tprot_mat$SNAP.contrast[!is.na(tprot)])>2)){
      
      # Compare groups with linear model; covariates age + gender+ apoe
      tres=lm(tprot~AGE.centered+Gender+APOEdich+SNAP.contrast,data=tprot_mat)
      # With emmeans we can extract the marginal means
      temmeans=test(emmeans(tres,pairwise~SNAP.contrast),adjust='none')
      
      # add the emm (se) to the table: Adjust columns
      res_mat[i, c(12,14,16,18,20, 22, 24)]=paste(round(temmeans$emmeans$emmean,2),' (',
                                                  round(temmeans$emmeans$SE,2),
                                                  ')',sep='')
      
      # add the differences & the pvalues between the groups
      res_mat[i, c(25,27,29,31,33,35,37,39,41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65)]=round(temmeans$contrasts$estimate,2)
      res_mat[i, c(26,28,30,32,34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66)]=round(temmeans$contrasts$p.value,5)
      
    }
  }
  
}


# save the results as a text file
write.table(res_mat,'resultsproteomics_CVD.txt',sep='\t')