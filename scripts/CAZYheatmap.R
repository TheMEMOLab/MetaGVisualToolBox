#!/usr/bin/RScript
###############################################################################
#   Visualization of CAZy genes information of Metagenome
# Assembled Genomes (MAGS) annotated by DRAM.
# 
# This script parse and sort the information of DRAM annotation and 
# destillation to create a sorted heatmap of CAZy genes present in 
# the MAGs. 
#
# Dependencies:
#           -R (tidyverse, readxl, gtools, vectrs, pheatmap, paletteer)
# Input tables:
#           -CheckM tabular 
#           -GTDBTk tabular bacterial classification
#           -GTDBTk tabular archeal classification
#           -metabolism_summary.xlsx summary from DRAM
#
# Author: Arturo Vera-Ponce de Leon
# contact: arturo.vera.ponce.de.leon@nmbu.no
# v 1.0
#toto
##############################################################################

## ----Libraries---------------------------------------------------------------
if (!require("tidyverse")){
  install.packages("tidyverse",
                   repos = "http://cran.us.r-project.org")
  library(tidyverse)
}

if (!require("readxl")){
  install.packages("readxl",
                   repos = "http://cran.us.r-project.org")
  library(readxl)
}

if (!require("vctrs")){
  install.packages("vctrs",
                   repos = "http://cran.us.r-project.org")
  library(vctrs)
}

if (!require("gtools")){
  install.packages("gtools",
                   repos = "http://cran.us.r-project.org")
  library(gtools)
}

if (!require("pheatmap")){
  install.packages("pheatmap",
                   repos = "http://cran.us.r-project.org")
  library(pheatmap)
}

if (!require("paletteer")){
  install.packages("paletteer",
                   repos = "http://cran.us.r-project.org")
  library(paletteer)
}

## ----Reading arguments-----------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

Dram <- args[1] #Metabolic Summary from DRAM
CHECK <- args[2] #CheckM results table
BAC <- args[3]  #GTDBTk Bacterial classification
ARCH <- args[4] #GTDBTk Archeal classification
OUT <- args[5] #Output name


## ----Reading DRAM------------------------------------------------------------
print("Reading DRAM....")

DRAM <- read_excel(Dram,
                   sheet = "carbon utilization") %>%
  filter(grepl("CAZY",header))

## ----Reading checkM-----------------------------------------------------
print("Reading checkM....")

CheckM <- read_tsv(CHECK) %>%
  dplyr::select(matches("Bin Id|Com|Con|Strain")) %>%
  rename_all(.,~str_replace_all(.," ","_")) %>%
  mutate(Genome=str_remove_all(Bin_Id,".fasta$")) %>%
  select(Genome, everything())

## ----load_Taxo_info----------------------------------------------------------

## ----reading_GTDBTK-----------------------------------------------------

print("Reading GTDBTk....")

gtdbtkBac <- read_tsv(BAC) %>%
  dplyr::select(user_genome,classification) %>%
  rename("Genome"=user_genome)

gtdbtkArch <- read_tsv(ARCH)%>%
  dplyr::select(user_genome,classification)%>%
  rename("Genome"=user_genome)

GTDBTKTotal <- gtdbtkBac %>%
  bind_rows(gtdbtkArch)


## -----Merge GTDBTK with CheckM------------------------------------------

GenoTaxoInfo <- GTDBTKTotal%>%
  left_join(CheckM,by="Genome")

GenoTaxoInfo <- GenoTaxoInfo %>%
  select(Bin_Id,Genome,
         classification,
         Completeness,
         Contamination,
         Strain_heterogeneity)%>%
  separate(classification,sep=";",
           into=c("D",
                  "P",
                  "C",
                  "O",
                  "Fa",
                  "G",
                  "S"),
           remove = FALSE) %>%
  mutate_at(.vars=vars(S),
            .funs = ~
              str_remove(.,
                         pattern = ".*[[:space:]]")) %>%
  mutate_at(.vars=vars(S),
            .funs = ~
              str_remove(.,
                         pattern = "(?<=sp)\\d.*")) %>%
  mutate_at(.vars = vars(D),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*d__")) %>%
  mutate_at(.vars = vars(P),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*p__")) %>%
  mutate_at(.vars = vars(C),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*c__")) %>%
  mutate_at(.vars = vars(O),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*o__")) %>%
  mutate_at(.vars = vars(Fa),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*f__")) %>%
  mutate_at(.vars = vars(G),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*g__")) %>%
  mutate_at(.vars = vars(S),
            .funs = ~
              str_remove_all(.,
                             pattern = ".*s__")) %>%
  mutate_at(.vars = vars(D,P,C,O,Fa,G,S),
            .funs = ~
              str_remove_all(.,
                             pattern = "_..*")) %>%
  mutate_at(.vars = vars(D,P,C,O,Fa,G,S),
            .funs = ~
              str_replace(.,
                          pattern = "^$",
                          replacement = ";")) %>%
  unite("classification",D:P:C:O:Fa:G:S,sep = "_") %>%
  mutate_at(.vars = vars(classification),
            .funs = ~
              str_replace(.,
                          pattern = "_;.*",
                          replacement = ""))%>%
  separate(classification,sep = "_",
           remove = FALSE,
           into = c("Domain",
                    "Phylum",
                    "Class",
                    "Order",
                    "Family",
                    "Genus",
                    "Species")) %>%
  mutate(Genus=coalesce(Genus,Order)) %>%
  dplyr::arrange(Phylum)



#Select the genomes 

print("Selecting and sorting CAZY genes from DRAM....")

CAZY <- DRAM %>%
  select(matches(str_c(vctrs::vec_c(colnames(DRAM)[1:5],
                                    GenoTaxoInfo$Genome),
                       collapse = "|")))%>%
  select(!module) %>%
  select(!header)



## ----sorting_CAzy------------------------------------------------------------

CAZYMod <- CAZY[gtools::mixedorder(CAZY$gene_id),]


## ----cazy_target-------------------------------------------------------------

subheader <- CAZYMod %>%
  select(gene_id,gene_description,subheader) %>%
  mutate_at(.vars = vars(subheader), 
            .funs = ~ gsub("\\,..*","",.)) %>%
  mutate(Glycan=if_else(grepl("Cellulose",subheader),"Cellulose",
                       if_else(grepl("Pectin|Rhamnose|pectic",subheader),"Pectin",
                               if_else(grepl("Chitin Backbone|Chitin Oligo",subheader),"Chitin",
                                       if_else(grepl("Starch",subheader),"Starch",
                                               if_else(grepl("Crystalline Cellulose",subheader),"LPMO",
                                                       if_else(grepl("Hemicellulose",subheader),"Hemicellulose",
                                                               if_else(grepl("Mucin|Fucose",subheader),"Mucin",
                                                                       if_else(grepl("sialidase",gene_description),"Mucin",
                                                                               "Other"))))))))) %>%
  select(gene_id,subheader,Glycan)

myorder <- c("Starch",
             "Cellulose",
             "Hemicellulose",
             "Pectin",
             "Chitin",
             "Mucin",
             "Other")




## ----log2 to adjust values-------------------------------------------

CAZYforM <- CAZYMod %>%
  select(-gene_description,-subheader) %>%
  filter(grepl("AA|GH|PL",gene_id)) %>%
  filter(if_any(where(is.numeric), ~ .x > 0))

CAzyForMMoudles <- CAZYforM %>%
  inner_join(subheader,by="gene_id") %>%
  arrange(match(Glycan,myorder)) %>%
  mutate_if(is.numeric,~log2(.+1)) 



## ----------------------------------------------------------------------------
CAzyMatrix <- CAzyForMMoudles %>%
  select(-c(subheader,Glycan)) %>%
  column_to_rownames("gene_id")


## -----Sort GenoTaxo and Matrix------------------------------------------------------------

GenoTaxoInfo <- GenoTaxoInfo %>% 
  filter(Genome %in% colnames(CAzyMatrix))

CAzyMatrix <- CAzyMatrix %>%
  select_(.dots = vars(GenoTaxoInfo$Genome))



## ----Phylum_colors by rows---------------------------------------------
heatTaxa <- GenoTaxoInfo %>%
  select(Genome,Phylum) %>%
  arrange(Phylum) %>%
  column_to_rownames("Genome")

colorheat <- paletteer_c("grDevices::Dark 3",
                         length(as.vector(na.omit(unique(GenoTaxoInfo$Phylum)))))

heatmap.colors <- as.character(colorheat)
names(heatmap.colors) <- sort(unique(GenoTaxoInfo$Phylum))


## ----Glycan colors---------------
phAnnot <- CAzyForMMoudles %>%
  select(gene_id,Glycan) %>%
  column_to_rownames("gene_id")


phAnnot$Glycan <- factor(phAnnot$Glycan, levels = myorder)
myGlycanColors <-  rev(RColorBrewer::brewer.pal(7,"Set3"))




## ----------------------------------------------------------------
annot_colors <- list(Glycan=c(Starch=myGlycanColors[1],
                             Cellulose=myGlycanColors[2],
                             Hemicellulose=myGlycanColors[3],
                             Pectin=myGlycanColors[4],
                             Chitin=myGlycanColors[5],
                             Mucin=myGlycanColors[6],
                             Other=myGlycanColors[7]),
                     Phylum=heatmap.colors)


## ----Use pheatmap to plot-----------------------------------
print("Ploting heatmap...")

Colors <- RColorBrewer::brewer.pal(n=8,"Greys")

TotalCAZyPH <- pheatmap(CAzyMatrix,
               color =Colors,
               cluster_rows = F,
               cluster_cols = F,
               annotation_row =  phAnnot,
               annotation_col =  heatTaxa,
               annotation_colors = annot_colors,
               cellwidth = 8,
               border_color = F,
               show_colnames = F)

## ----Save into a pdf-----------------------------------
print(paste0("Saving figure to: ",OUT,".pdf"))

ggsave(TotalCAZyPH,
       file=paste0(OUT,".pdf"),
       width = 30,
       height = 30)



