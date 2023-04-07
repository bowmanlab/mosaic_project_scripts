## This code represents a secondary QC for all Bowman Lab MOSAiC DNA samples post phylogenetic placement and metabolic inference using paprica v 0.7.2 and ROPE-RDP classification of paprica edges (https://github.com/avishekdutta14/ROPE) In this script, we build a comprehensive taxonomy map and remove extraneous mitochondrial sequences and chloroplasts. 

## Code written by Emelia J Chamberlain

## Code updated April 2023

########## Packages  Required ###########
library(vegan)
library(dplyr)
library(tidyr)
library(readr)

########## Load Data ##########

## load R workspace from script 1
load("data_QC/sequences/generated_data/mosaic_dna_QC_output.RData")

########## Create Full Taxon Maps from Paprica Output ########## 

## prokaryotes
#isolate edges associated with ASVs
edges.pro <- map.16s$global_edge_num

#pull out domains of said ASVs, fill in blanks
prodomain <- taxa.16s[edges.pro,1]
dobac <- prodomain[1:26373]
doarc <- prodomain[26374:27231]
dobac[which(dobac == "")] <- "Bacteria"
dobac[is.na(dobac)] <- "Bacteria" 
doarc[which(doarc == "")] <- "Archaea"
doarc[is.na(doarc)] <- "Archaea "
map.16s$domain <- c(dobac,doarc)

#clean environment
rm(prodomain,doarc,dobac)

#pull out phylum, fill in blanks with last assignment
prophylum <- taxa.16s[edges.pro,2]
phybac <- prophylum[1:26373]
phyarc <- prophylum[26374:27231]
phybac[which(phybac == "")] <- map.16s[which(phybac == ""),4]
phybac[is.na(phybac)] <- map.16s[is.na(phybac),4]
phyarc[which(phyarc == "")] <- map.16s[which(phyarc == ""),4]
phyarc[is.na(phyarc)] <- map.16s[is.na(phyarc),4]
map.16s$phylum <- c(phybac,phyarc)

#clean environment
rm(prophylum, phybac, phyarc)

#pull out class, fill in blanks with last assignment 
proclass <- taxa.16s[edges.pro,4]
classbac <- proclass[1:26373]
classarc <- proclass[26374:27231]
classbac[which(classbac == "")] <- map.16s[which(classbac == ""),5]
classbac[is.na(classbac)] <- map.16s[is.na(classbac),5]
classarc[which(classarc == "")] <- map.16s[which(classarc == ""),5]
classarc[is.na(classarc)] <- map.16s[is.na(classarc),5]
map.16s$class <- c(classbac,classarc)

#clean environment
rm(proclass, classbac, classarc)

#pull out family, fill in blanks with last assignment 
profam <- taxa.16s[edges.pro,6]
fambac <- profam[1:26373]
famarc <- profam[26374:27231]
fambac[which(fambac == "")] <- map.16s[which(fambac == ""),6]
fambac[is.na(fambac)] <- map.16s[is.na(fambac),6]
famarc[which(famarc == "")] <- map.16s[which(famarc == ""),6]
famarc[is.na(famarc)] <- map.16s[is.na(famarc),6]
map.16s$family <- c(fambac,famarc)

#clean environment
rm(profam, fambac, famarc)

#pull out genus, fill in blanks with last assignment 
progenus <- taxa.16s[edges.pro,7]
genbac <- progenus[1:26373]
genarc <- progenus[26374:27231]
genbac[which(genbac == "")] <- map.16s[which(genbac == ""),7]
genbac[is.na(genbac)] <- map.16s[is.na(genbac),7]
genarc[which(genarc == "")] <- map.16s[which(genarc == ""),7]
genarc[is.na(genarc)] <- map.16s[is.na(genarc),7]
map.16s$genus <- c(genbac,genarc)

#clean environment
rm(progenus, genbac, genarc)

#pull out taxon (or closest association), fill in blanks with last assignment
protaxon <- taxa.16s[edges.pro,10]
taxonbac <- protaxon[1:26373]
taxonarc <- protaxon[26374:27231]
taxonbac[which(taxonbac == "")] <- map.16s[which(taxonbac == ""),8]
taxonbac[is.na(taxonbac)] <- map.16s[is.na(taxonbac),8] 
taxonarc[taxonarc == ""] <- map.16s[which(taxonarc == ""),8]
taxonarc[is.na(taxonarc)] <- map.16s[is.na(taxonarc),8] 
map.16s$taxon <- c(taxonbac,taxonarc)

#clean environment
rm(protaxon, taxonbac, taxonarc)
rm(edges.pro)

## eukaryotes
#isolate edges associated with ASVs
edges.euk <- map.18s$global_edge_num

#provide domain of said ASVs
map.18s$domain <- "Eukaryota"

#pull out phylum (division from paprica), fill in blanks with last assignment 
phyeuk <- taxa.18s[edges.euk,3]
phyeuk[which(phyeuk == "")] <- map.18s[which(phyeuk == ""),4]
phyeuk[is.na(phyeuk)] <- map.18s[is.na(phyeuk),4]
map.18s$phylum <- phyeuk

#clean environment
rm(phyeuk)

#pull out class, fill in blanks with last assignment
classeuk <- taxa.18s[edges.euk,4]
classeuk[which(classeuk == "")] <- map.18s[which(classeuk == ""),5]
classeuk[is.na(classeuk)] <- map.18s[is.na(classeuk),5]
map.18s$class <- classeuk

#clean environment
rm(classeuk)

#pull out family, fill in blanks with last assignment 
fameuk <- taxa.18s[edges.euk,6]
fameuk[which(fameuk == "")] <- map.18s[which(fameuk == ""),6]
fameuk[is.na(fameuk)] <- map.18s[is.na(fameuk),6]
map.18s$family <- fameuk

#clean environment
rm(fameuk)

#pull out genus, fill in blanks with last assignment 
geneuk <- taxa.18s[edges.euk,7]
geneuk[which(geneuk == "")] <- map.18s[which(geneuk == ""),7]
geneuk[is.na(geneuk)] <- map.18s[is.na(geneuk),7]
map.18s$genus <- geneuk

#clean environment
rm(geneuk)

#pull out taxon, fill in blanks with last assignment  
taxeuk <- taxa.18s[edges.euk,9]
taxeuk[which(taxeuk == "")] <- map.18s[which(taxeuk == ""),8]
classeuk[is.na(classeuk)] <- map.18s[is.na(taxeuk),8]
map.18s$taxon <- taxeuk

#clean environment
rm(taxeuk)
rm(edges.euk)

#backup these new maps as a csv file
write.csv(map.16s, "data_QC/sequences/generated_data/mosaic_taxaQC_16s_seq_edge_map.csv")
write.csv(map.16s, "data_QC/sequences/generated_data/mosaic_taxaQC_18s_seq_edge_map.csv")

########## Compare with ROPE Classification ##########
#### bacteria

#reading taxon map from ROPE which contains unique ID mapped to taxonomy from ROPE
ROPE_tax = read.csv("~/mosaic/new_2023/Data_QC/DNA/ROPEbac/taxa_map_ROPE_unique.csv", header=TRUE)
names(ROPE_tax)[names(ROPE_tax) == "Unique."] <- "UniqueID"

# reading ROPE unique tally which contains the unique ID mapped to unique sequence
ROPE_seq = read.csv ("~/mosaic/new_2023/Data_QC/DNA/ROPEbac/unique_ID_tally.csv", header= TRUE)
ROPE_seq = ROPE_seq[,1:2]
merge2 = merge(ROPE_seq,ROPE_tax,by= "UniqueID",all=TRUE)

## rope comparison
compare <- read.csv("~/mosaic/new_2023/Data_QC/DNA/ROPEbac/comparison_phylum_RP.csv")
rownames(compare) <- compare$sequences
merge2$comparison <- compare[merge2$sequences,7]

## final bacteria dataframe
bacteria_rope <- merge2

## archaea
#reading taxon map ROPE which contains unqiue ID mapped to taxonomy from ROPE
ROPE_tax = read.csv("~/mosaic/new_2023/Data_QC/DNA/ROPE_arc/taxa_map_ROPE_unique.csv", header=TRUE)
names(ROPE_tax)[names(ROPE_tax) == "Unique."] <- "UniqueID"

# reading ROPE unique tally which contains the unique ID mapped to unique sequence
ROPE_seq = read.csv ("~/mosaic/new_2023/Data_QC/DNA/ROPE_arc/unique_ID_tally.csv", header= TRUE)
ROPE_seq = ROPE_seq[,1:2]
merge2 = merge(ROPE_seq,ROPE_tax,by= "UniqueID",all=TRUE)

## rope comparison
compare <- read.csv("~/mosaic/new_2023/Data_QC/DNA/ROPE_arc/comparison_phylum_RP.csv")
rownames(compare) <- compare$sequences
merge2$comparison <- compare[merge2$sequences,7]

## final archaea dataframe
archaea_rope <- merge2

rm(compare,merge2,ROPE_seq,ROPE_tax) #clean environment

## finalize ROPE files
ROPE.16s <- rbind(archaea_rope,bacteria_rope) #merge bacteria and archaea
rm(archaea_rope,bacteria_rope) #clean environment
colnames(ROPE.16s) <- paste0("rope_",colnames(ROPE.16s)) #distinguish colnames from paprica results

## combine ROPE and paprica files by sequence reads

map.16s$sequence <- rownames(map.16s) #add sequence column to mapfile

fullmap.16s <- left_join(map.16s, ROPE.16s, by = c("sequence" = "rope_sequences")) #join files

##19142 ASV placement matches between paprica and ROPE
##8089 ASV placement mis-matches between paprica and ROPE

#clean environment
rm(map.16s, ROPE.16s)

#### eukaryotes
#reading taxon map ROPE which contains unique ID mapped to taxonomy from ROPE
ROPE_tax = read.csv("~/mosaic/new_2023/Data_QC/DNA/ROPE18S/taxa_map_ROPE_unique_18S.csv", header=TRUE)
names(ROPE_tax)[names(ROPE_tax) == "Unique."] <- "UniqueID"

# reading ROPE unique tally which contains the unique ID mapped to unique sequence
ROPE_seq = read.csv ("~/mosaic/new_2023/Data_QC/DNA/ROPE18S/unique_ID_tally.csv", header= TRUE)
ROPE_seq = ROPE_seq[,1:2]
ROPE.18s = merge(ROPE_seq,ROPE_tax,by= "UniqueID",all=TRUE)

colnames(ROPE.18s) <- paste0("rope_",colnames(ROPE.18s)) #distinguish colnames from paprica results

rm(ROPE_seq,ROPE_tax) #clean environment

## combine ROPE and paprica files by sequence reads
map.18s$sequence <- rownames(map.18s)#add sequence column to mapfile
fullmap.18s <- left_join(map.18s, ROPE.18s, by = c("sequence" = "rope_sequences")) #join files

## no comparison script for 18s yet

## write final files
write.csv(fullmap.16s, "data_QC/sequences/generated_data/mosaic_taxaQC_full16s_seq_edge_map.csv")
write.csv(fullmap.16s, "data_QC/sequences/generated_data/mosaic_taxaQC_full18s_seq_edge_map.csv")

#clean environemnt
rm(map.18s, ROPE.18s)

########## Remove Suspected Mitochondrial Amplicons from prokaryotes ##########

#remove all suspect taxa candidatus nasuia deltocephalinicola or carsonella ruddi   

#isolate just suspect ASVs which are likely be mitochondrial DNA from paprica assignments
Mito1 <- grep('Carsonella', fullmap.16s[, 9]) # 2061 ASVs assigned
Mito2 <- grep('Nasuia', fullmap.16s[, 9]) # no contam.

##combine
allMitos1 <- c(Mito1,Mito2) #generate list of all suspected mitochondrial ASVs
rm(Mito1,Mito2) #clean environment

##remove from map
fullmap.16s.select <- fullmap.16s[-allMitos1,]

## check ROPE classifications 
Mito1 <- grep('Carsonella', fullmap.16s.select[, 18]) # 533 ASVs assigned
Mito2 <- grep('Nasuia', fullmap.16s.select[, 18]) # no contam
##combine
allMitos2 <- c(Mito1,Mito2) #generate list of all suspected mitochondrial ASVs
rm(Mito1,Mito2) #clean environment

##remove from map
fullmap.16s.select2 <- fullmap.16s.select[-allMitos2,]

## Suspect taxa strains that were identified as candidatus carsonella by ROPE-RDP (isolated from command line run of fullmap.16s.select[allMitos2,9] and researched)
#Candidatus Nardonella dryophthoridicola
#Mycoplasmopsis pulmonis/Mycoplasmataceae
#Candidatus Gromoviella agglomerans
#Neorickettsia
#Candidatus Hydrogenosomobacter endosymbioticus
#Candidatus Fokinia solitaria


## check if there are more of these suspect assignments from paprica 
Mito1 <- grep('Candidatus Nardonella dryophthoridicola', fullmap.16s.select2[, 9]) # 15 more ASVs - classified by ROPE to eukaryotes so definitely some sort of organelle sequence
Mito2 <- grep('Mycoplas', fullmap.16s.select2[, 9]) # 194 more ASVs - classified by ROPE across a wide range of taxa with fairly low confidence levels but including lots o feukaryotes, extremeophiles, and Nitrosopumilus
Mito3 <- grep('Neorickettsia', fullmap.16s.select2[, 9]) #86 more ASVs - classified by ROPE across a wide range of taxa, including marine bacterium, but all with relatively low confidence values
Mito4 <- grep('Candidatus Hydrogenosomobacter endosymbioticus', fullmap.16s.select2[, 9]) #55 more ASVs with similar ROPE assignments to Mito3
Mito5 <- grep('Candidatus Fokinia solitaria', fullmap.16s.select2[, 9]) # 47 mor ASVs mostly assigned by ROPE to the endosymbiont genus Lyticum, indicating high likelyhood of organelle sequences. 

##combine 
allMitos3 <- c(Mito1,Mito2,Mito3,Mito4,Mito5) #generate list of all suspected mitochondrial ASVs
rm(Mito1,Mito2,Mito3,Mito4,Mito5) #clean environment

##remove from map
fullmap.16s.select3 <- fullmap.16s.select2[-allMitos3,]

#### Here we will go with a highly conservative QC removing all 3 levels of suspected mitochondrial/organelle sequences from our final map file (total of 2991 ASVs removed)
map.16s.select <- fullmap.16s.select3 #final 16s map file


## Save information for all suspected sequences for potential use in later analyses
mitofile <- rbind(fullmap.16s[allMitos1,],fullmap.16s.select[allMitos2,],fullmap.16s.select2[allMitos3,])

#taxonomic placement
write.csv(mitofile, file = "data_QC/sequences/generated_data/mosaic_suspected_mitochondria_seq_edge_map.csv") 
#unique tallies
write.csv(as.data.frame(unique.16s.select[,mitofile$sequence]), file = "data_QC/sequences/generated_data/mosaic_suspectmitochondria_ASVtally.csv") 

#remove suspect sequences from ASV tally information
unique.16s.select <- unique.16s.select[,map.16s.select$sequence]

## clean environment
rm(fullmap.16s, fullmap.16s.select,fullmap.16s.select2,fullmap.16s.select3, mitofile)
rm(allMitos1,allMitos2,allMitos3)

########## Isolate Suspected Chloroplast Amplicons from prokaryotes ##########

## save mapfile containing chloroplast and cyanobacterial sequences 
map.16s.select.with.chloroplasts <- map.16s.select 
write.csv(map.16s.select.with.chloroplasts, file = "data_QC/sequences/generated_data/mosaic_16S_select_withchloroplasts_seq_edge_map.csv") 
unique.16s.select.with.chloroplasts <- unique.16s.select

#### Isolate Cyanobacteria and Chloroplasts
chloro1 <- grep('Cyanobacteria', map.16s.select[, 5]) #1510 cyanobacteria sequences identified by paprica

map.16s.select2 <- map.16s.select[-chloro1,] #remove sequences

chloro2 <- grep('Cyanobacteria', map.16s.select2[, 18]) #593 more cyanobacteria/chloroplast sequences identified by ROPE
#These were assigned mostly paprica taxonomies with small genomes (i.e. Glaciecola amylolytica) or other endosymbiont strains like Candidatus Gromoviella agglomerans, Wolbachia endosymbiont of Corcyra cephalonica, Elusimicrobium minutum Pei191, etc. as well as some small genome archaeal groups. 

map.16s.select3 <- map.16s.select2[-chloro2,] #remove sequences

#### Here we will go with a highly conservative QC isolating both levels of cyanobacteria and chloroplast sequences (total of 2103 sequences

## Save information for all cyanos/chloroplasts for potential use in later analyses
chlorofile <- rbind(map.16s.select[chloro1,],map.16s.select2[chloro2,])

#taxonomic placement
write.csv(chlorofile, file = "data_QC/sequences/generated_data/mosaic_cyanobacteria_chloroplasts_seq_edge_map.csv") 
#unique tallies
write.csv(as.data.frame(unique.16s.select[,chlorofile$sequence]), file = "data_QC/sequences/generated_data/mosaic_cyanobacteria_chloroplasts_ASVtally.csv") 

## remove chloroplasts and cyanos from final ASV tally information
map.16s.select <- map.16s.select3
unique.16s.select <- unique.16s.select[,map.16s.select$sequence]

## clean environment
rm(map.16s.select2, map.16s.select3, chlorofile)
rm(chloro1,chloro2)

########## eukaryotes ##########

#at this time there is no formal procedure or specific taxa to look out for when QC'ing 18s data within our lab

#there is also no mismatch script for comparison between ROPE-RDP and paprica output. 

map.18s.select <- fullmap.18s



########## Explore remaining mis-matches between the two pipelines
mismatch.16s <- map.16s.select[which(map.16s.select$rope_comparison == "mis-match"),] ## 4859 remaining mis-matches wihtin the 16s data, some of which appear to be suspect sequences? Should these also be removed? 
write.csv(mismatch.16s, "data_QC/sequences/generated_data/mosaic_16s_paprica_ROPE_mismatch.csv")

########## Clean Environment and Save Workspace ##########
##clean environement 
rm(mismatch.16s, fullmap.18s) 

##save environment
save.image(file = "data_QC/sequences/generated_data/mosaic_taxa_QC_output.RData")

## move to script 3

