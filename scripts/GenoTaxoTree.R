#!/usr/bin/RScript
###############################################################################
#   Visualization of Phylogenomic and quality information of Metagenome
# Assembled Genomes (MAGS).
# 
# This script gathers the quality information (completeness and contamination), 
# from CheckM and merge it with the phylogenetic tree obtained by Phyloplan.
# 
# Dependencies:
#           -R (tidyverse,GGtree,ape,tidytree,GGtreeExtra,paletteer)
# Input tables:
#           -CheckM tabular 
#           -GTDBTk tabular bacterial classification
#           -GTDBTk tabular archeal classification
#           -RAxML tree from Phyloplhan
#
# Author: Arturo Vera-Ponce de Leon
# contac: arturo.vera.ponce.de.leon@nmbu.no
# v 1.0
#
##############################################################################


## -----------------------------------------------------------------------

#Libraries

if (!require("tidyverse")){
  install.packages("tidyverse",
                   repos = "http://cran.us.r-project.org")
  library(tidyverse)
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager",
                   repos = "http://cran.us.r-project.org")
}
  
if(!require("ggtree")){
  BiocManager::install("ggtree")
  library(ggtree)
}

if (!require("tidytree")){
  install.packages("tidytree",
                   repos = "http://cran.us.r-project.org")
  library(tidytree)
}

if(!require("ggtreeExtra")){
  BiocManager::install("ggtreeExtra")
  library(ggtree)
}

if (!require("ape")){
  install.packages("ape",
                   repos = "http://cran.us.r-project.org")
  library(ape)
}

if (!require("paletteer")){
  install.packages("paletteer",
                   repos = "http://cran.us.r-project.org")
  library(paletteer)
}

## ----Reading arguments-----------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
CHECK <- args[1] #CheckM results table
BAC <- args[2]  #GTDBTk Bacterial classification
ARCH <- args[3] #GTDBTk Archeal classification
TREE <- args[4] #Phyloplhlan RAXMl tree
OUT <- args[5] #Output name

## ----Reading checkM-----------------------------------------------------
print("Reading checkM....")

CheckM <- read_tsv(CHECK) %>%
  dplyr::select(matches("Bin Id|Com|Con|Strain")) %>%
  rename_all(.,~str_replace_all(.," ","_")) %>%
  mutate(Genome=str_remove_all(Bin_Id,".fasta$")) %>%
  select(Genome, everything())



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




## -----------------------------------------------------------------------
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



## -----------------------------------------------------------------------
print("Reading tree...")

Tree <- read.tree(TREE)

GGT <- ggtree(Tree,
              layout ="circular")  %<+% GenoTaxoInfo

GGT


## -----------------------------------------------------------------------
Domains <- GenoTaxoInfo %>%
  group_by(Domain) %>%
  summarise(across(,list)) %>%
  select(Domain,Bin_Id)

Dom <- map(Domains$Bin_Id, function(x){
  unlist(strsplit(x,split=","))
  })

names(Dom) <- Domains$Domain


## -----------------------------------------------------------------------
print("Plotting tree...")

GG <- groupOTU(Tree, Dom, "Domains")

GGT <- ggtree(GG,aes(color=Domains),
              layout = "fan",
              open.angle=12,
              size=0.2) %<+% GenoTaxoInfo +
  scale_color_brewer(palette = "Set2") +
    geom_treescale(fontsize = 3)

## -----------------------------------------------------------------------
#Adding quality (complenteness and contamination) heatmaps

GGT2 <- GGT + 
  ggtreeExtra::geom_fruit(
  geom=geom_tile,
  mapping=aes(fill=Completeness),
  width=0.1,
  color="white",
  pwidth=0.1,
  offset=0.10)+ 
  scale_fill_gradientn(name="GenomeCompletness %",
                       colors = RColorBrewer::brewer.pal(9, "YlOrRd"))+
  ggnewscale::new_scale_fill()+
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Contamination),
    width=0.1,
    color="white",
    offset=0.12
  )

GGT2


## -----------------------------------------------------------------------
##Preparing heatmap palette colors for the Phyla from GTDBTK 

colorheat <- paletteer_c("grDevices::Dark 3",
                         length(as.vector(na.omit(unique(GGT$data$Phylum)))))

heatmap.colors <- as.character(colorheat)
names(heatmap.colors) <- as.factor(na.omit(unique(GGT$data$Phylum)))


## -----------------------------------------------------------------------

GGT3 <- GGT2 +
  ggnewscale::new_scale_fill() +
geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Phylum),
    width=0.1,
    color="white",
    offset=0.14)+
  scale_fill_manual(breaks =as.factor(na.omit(unique(GGT$data$Phylum))),
                    values = heatmap.colors,
                    name="Taxonomy_Phylum")


## -----------------------------------------------------------------------
print(paste0("Saving figure to: ",OUT,".pdf"))

ggsave(GGT3,
       file=paste0(OUT,".pdf"),
       width = 15,
       height = 12)


