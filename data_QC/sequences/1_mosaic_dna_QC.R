## This code represents an initial QC for all Bowman Lab MOSAiC DNA samples. Raw sequence reads (available at https://www.ncbi.nlm.nih.gov/bioproject/895866) have already been denoised and run through the paprica pipeline (https://github.com/bowmanjeffs/paprica) for phylogenetic placement and metabolic inference (v 0.7.2)

## Code written by Emelia J Chamberlain, adapted from tutotials by Jeff Bowman (polarmicrobes.org)

## Code updated April 2023

########## Packages  Required ##########
library(vegan)
#library(oce)

########## Read in data ##########

## load R workspace
load("~/mosaic/mosaic_16S_18S/data/round2/20220919_mosaic.Rdata")

########## General QC ##########

## Join the bacteria and archaea datasets - since we are using a cross-domain primer, it therefore makes sense to analyze these data together
unique.16S <- merge(unique.bac, unique.arc, by = 0, all = T)
tally.16S <- merge(tally.bac, tally.arc, by = 0, all = T)

row.names(unique.16S) <- unique.16S$Row.names
row.names(tally.16S) <- tally.16S$Row.names

unique.16S$Row.names <- NULL
tally.16S$Row.names <- NULL

## convert NAs to 0
unique.16S[is.na(unique.16S)] <- 0
tally.16S[is.na(tally.16S)] <- 0

## Generally we want to avoid libraries with fewer than 5000 reads to normalize sequencing depth among samples. However, instead of bulk removing at this stage, we will create flags for libraries below our chosen size threshold to inform our later analyses. These flags will be based off of the unique ASV tallies only. 

#first we will remove the average 16s gene copy number correction 
unique16.uc <- unique.16S * data.bac$n16S.mean 

#extract samples with fewer than 5000
bad.16slibs <- unique16.uc[rowSums(unique16.uc) < 5000,]
bad.euklibs <- unique.euk[rowSums(unique.euk) < 5000,]

# Flag A represents libraries between 4000 and 5000
# Flag B represents libraries between 1000 and 4000
# Flag C represents libraries < 1000

FA.16s <- rownames(bad.16slibs[rowSums(bad.16slibs) > 4000,])
FC.16s <- rownames(bad.16slibs[rowSums(bad.16slibs) < 1000,])
FB.16s <- bad.16slibs[rowSums(bad.16slibs) < 4000,]
FB.16s <- FB.16s[rowSums(FB.16s) > 1000,]
FB.16s <- rownames(FB.16s)

FA.18s <- rownames(bad.euklibs[rowSums(bad.euklibs) > 4000,])
FC.18s <- rownames(bad.euklibs[rowSums(bad.euklibs) < 1000,])
FB.18s <- bad.euklibs[rowSums(bad.euklibs) < 4000,]
FB.18s <- FB.18s[rowSums(FB.18s) > 1000,]
FB.18s <- rownames(FB.18s)

## Remove singleton ASVs
unique.16s.select <- unique.16S[,colSums(unique.16S) > 1]
unique.18s.select <- unique.euk[,colSums(unique.euk) > 1]

## clean environment
rm(unique16.uc)

########## Metadata Mapping ############
## read in metadata and synchronize with sample names
meta <- read.csv("~/mosaic/mosaic_16S_18S/data/FILTRATION2_DNA_ONLY_updated1122.csv", header = TRUE)

meta$libraryname_pro <- NA
for (i in 1:1297) {
  meta[i,32] <- paste(meta[i,4], "_16S.exp.", sep = "")
}

meta$libraryname_euk <- NA
for (i in 1:1297) {
  meta[i,33] <- paste(meta[i,4], "_18S.exp.", sep = "")
}

## Add QC flags
meta$proQC <- NA
meta$eukQC <- NA

#prokaryotes
rownames(meta) <- meta$libraryname_pro
meta[FA.16s, 34] <- "A"
meta[FB.16s, 34] <- "B"
meta[FC.16s, 34] <- "C"

#eukaryotes
rownames(meta) <- meta$libraryname_euk
meta[FA.18s, 35] <- "A"
meta[FB.18s, 35] <- "B"
meta[FC.18s, 35] <- "C"

## sanity check, reconcile numbers of samples 
prolibnamestally <- rownames(unique.16s.select)
prolibnamesmeta <- meta$libraryname_pro

euklibnamestally <- rownames(unique.18s.select)
euklibnamesmeta <- meta$libraryname_euk

#isolate old, incorrect library builds 
tabspro <- prolibnamestally %in% prolibnamesmeta
tabseuk <- euklibnamestally %in% euklibnamesmeta
baddies_pro <- prolibnamestally[which(tabspro == FALSE)]
baddies_euk <- euklibnamestally[which(tabseuk == FALSE)]


prolibnamestally <- prolibnamestally[-which(tabspro == FALSE)] #remove bad libraries from list
euklibnamestally <- euklibnamestally[-which(tabseuk == FALSE)] #remove bad libraries from list

#determine which samples are missing if any 
tabspro2 <- prolibnamesmeta %in% prolibnamestally
table(tabspro2) #1 sample officially missing

tabseuk2 <- euklibnamesmeta %in% euklibnamestally
table(tabseuk2) #1 sample officially missing

#extract these missing libraries
promisssing <- prolibnamesmeta[which(tabspro2 == FALSE)]
eukmisssing <- euklibnamesmeta[which(tabseuk2 == FALSE)]

## flag on data sheet
rownames(meta) <- meta$libraryname_pro
meta[promisssing, 34] <- "librarymissing"
rownames(meta) <- meta$libraryname_euk
meta[eukmisssing, 35] <- "librarymissing"

## save backup of final metadata
write.csv(meta, "data_QC/sequences/generated_data/mosaic_dnaQC_metadata_output.csv")

## Remove missing samples from working version
promisssing #read
eukmisssing #read
meta <- meta[grep('20200405_OC_N8_379_16S.exp.', row.names(meta), invert = T),]
meta <- meta[grep('20200907_FYI_DMSCore2_sec6_1220_18S.exp.', row.names(meta), invert = T),]

## match file orders between meta and tally sheets, removing bad/old libraries

#prokaryotes
rownames(meta) <- meta$libraryname_pro
unique.16s.select <- unique.16s.select[rownames(meta),]
tally.16s <- tally.16S[rownames(meta),]

#eukaryotes
rownames(meta) <- meta$libraryname_euk
unique.18s.select <- unique.18s.select[rownames(meta),]
tally.18s <- tally.euk[rownames(meta),]

## Match sample data files as well
data.bac.select <- data.bac[row.names(unique.16s.select),]
data.arc.select <- data.arc[row.names(unique.16s.select),]

## Add bacteria paprica output to metadata
namesbac <- colnames(data.bac.select)
bac <- paste0("bac_",namesbac)
namesarc <- colnames(data.arc.select)
arc <- paste0("arc_",namesarc)
colnames(data.bac.select) <- bac
colnames(data.arc.select) <- arc

meta <- cbind(meta, data.bac.select)
meta <- cbind(meta, data.arc.select)


## read in seq edge map and create a unique name for each ASV (good for downstream analysis)

#prokaryotes (bacteria & archaea)
map.16s <- rbind(map.bac, map.arc)
map.16s <- map.16s[colnames(unique.16s.select),]

uniquenamespro <- vector(mode = "character", length = nrow(map.16s))
for (i in 1:length(uniquenamespro)) {
  uniquenamespro[i] <- paste("Unique", i, sep = "")
}
map.16s$UniqueID <- uniquenamespro


write.csv(map.16s, "data_QC/sequences/generated_data/mosaic_dnaQC_16s_seq_edge_map.csv") #save a backup of this file

#eukaryotes
map.18s <- map.euk[colnames(unique.18s.select),]

uniquenameseuk <- vector(mode = "character", length = nrow(map.18s))
for (i in 1:length(uniquenameseuk)) {
  uniquenameseuk[i] <- paste("Unique", i, sep = "")
}
map.18s$UniqueID <- uniquenameseuk
write.csv(map.18s, "data_QC/sequences/generated_data/mosaic_dnaQC_18s_seq_edge_map.csv") #save a backup of this file

############## Merge Taxa Files ##################
taxa.16s <- rbind(taxa.bac, taxa.arc)
taxa.18s <- taxa.euk


############### Clean up environment ######################
#Here we will remove unnecessary items from the environment for future use
rm(bad.16slibs, bad.euklibs)
rm(data.arc, data.bac)
rm(data.arc.select,data.bac.select)
rm(map.arc, map.bac, map.euk)
rm(tally.arc, tally.bac, tally.euk)
rm(taxa.arc, taxa.bac, taxa.euk)
rm(unique.16S)
rm(unique.arc,unique.bac,unique.euk)
rm(arc,bac,baddies_euk,baddies_pro,euklibnamesmeta,euklibnamestally,eukmisssing)
rm(FA.16s,FA.18s,FB.16s,FB.18s,FC.16s,FC.18s)
rm(i,namesarc,namesbac,prolibnamesmeta,prolibnamestally,promisssing,tabseuk,tabseuk2,tabspro,tabspro2)
rm(uniquenameseuk,uniquenamespro)
rm(tally.16S)

########## Save final workspace for cont. analyses #########

save.image(file = "data_QC/sequences/generated_data/mosaic_dna_QC_output.RData")

## move to script 2