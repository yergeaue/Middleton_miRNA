###Load raw data, manipulate, normalize and save as intermediate files for other scripts

#mapping files
CAN_map <- read.table(file = here("data", "raw", "CAN_miRNA_map.tsv"), row.names = 1, header = T, sep = "\t", comment.char = "") #104 obs in 6 var
FR_ara_map <- read.table(file = here("data", "raw", "FR_ara_mapping_file.tsv"), row.names = 1, header = T, sep = "\t", comment.char = "") #35 obs in 5 var
FR_bra_map <- read.table(file = here("data", "raw", "FR_bra_mapping_file.tsv"), row.names = 1, header = T, sep = "\t", comment.char = "") #22 obs in 4 var

#miRNA abundance tables 
CAN_bra <- read.table(file = here("data", "raw", "CAN_Bdistachyon_Ref.tsv"), header = T, row.names = 1, sep = "\t")
CAN_ara <- read.table(file = here("data", "raw", "CAN_Athaliana_Ref.tsv"), header = T, row.names = 1, sep = "\t")
FR_bra <- read.table(file = here("data", "raw", "FR_bra_miRNA_abundance.tsv"), header = T, row.names = 1, sep = "\t")
FR_ara <- read.table(file = here("data", "raw", "FR_ara_miRNA_abundance.tsv"), header = T, row.names = 1, sep = "\t")

#Normalize miRNA abundance tables
CAN_bra_norm <-  data.frame(apply(CAN_bra, 1, "/", colSums(CAN_bra))) #104 obs. of 524 variables
CAN_ara_norm <-  data.frame(apply(CAN_ara, 1, "/", colSums(CAN_ara))) #104 obs. of 428 variables
FR_bra_norm <- data.frame(apply(FR_bra, 1, "/", colSums(FR_bra))) #5 obs. of 524 variables
FR_ara_norm <- data.frame(apply(FR_ara, 1, "/", colSums(FR_ara))) #35 obs. of 428 variables

#mature miRNA sequences (from mirBASE)
miRNA.fa <- read.table(file = here("data","raw", "mature.fa"), sep = "\t") #97 770 obs of 1 var.
miRNA.identity <- miRNA.fa[c(TRUE,FALSE),]
miRNA.seq <- miRNA.fa[c(FALSE,TRUE),]
miRNA.df <- data.frame(Identity=miRNA.identity, Sequence=miRNA.seq)
miRNA.df$Identity <- gsub(">","", miRNA.df$Identity)
miRNA.df$Identity <- gsub(" MIMAT","_MIMAT",miRNA.df$Identity)
miRNA.df <- separate(miRNA.df, Identity, into = c("Identity", "Genus", "species", "miRNA"), sep = " ")

#Microscopy images
vario <- readTIFF(source =  here("data", "raw", "Variovorax_zoom.tiff")) #Variovorax
vario.g <- rasterGrob(vario, interpolate = TRUE)
bacillus <- readTIFF(source =  here("data", "raw", "Bacillus_zoom.tiff")) #Bacillus
bacillus.g <- rasterGrob(bacillus, interpolate = TRUE)

#Transcriptomics
vario.transcripto <- read.table(file = here("data", "raw", "vario_merged_gene_abundance.tsv"), row.names = 1, header = TRUE)#variovorax 6047 obs of 12 variables
vario.info <- read.table(file = here("data", "raw", "info-table-vario.csv"), sep= ";", row.names = 1, header = TRUE)#variovorax 12 obs of 1 variables
vario.annotations <- read.table(here("data", "raw", "vario_annotations.tsv"), quote = "", comment.char = "", sep="\t", header = T)#variovorax 5936 obs of 33 variables
bacillus.transcripto <- read.table(file = here("data", "raw", "bacillus_merged_gene_abundance.csv"), sep = ";", row.names = 1, header = TRUE)#bacillus 4899 obs of 12 variables
bacillus.info <- read.table(file = here("data", "raw", "info-table-bacillus.csv"), sep= ";", row.names = 1, header = TRUE)#bacillus 12 obs of 1 variables

#miPEP and variovorax RT-qPCR
vario.qpcr.miPEP <- read.table(here("data", "raw", "qPCR-miPEP-variovorax.txt"), sep = "\t", header = T)

#mutant-community
mut.biom_file<-import_biom(here("data","raw", "16S.mut.otu_table_filtered.biom"), here("data","raw","16S.mut.tree.fasttree"), here("data", "raw", "16S.mut.otus.fasta"), parseFunction=parse_taxonomy_greengenes)
mut.map <- read.table(here("data","raw", "16S.mut.mapping_file.tsv"), sep="\t", quote = "", row.names = 1, comment.char = "", header = TRUE)
mut.map <- mut.map %>% mutate(Type = replace(Type, Type== "WT6.5", "WT"))
mut.map <- mut.map %>% mutate(Type = replace(Type, Type== "ago1", "ago1-27"))
mut.map <- mut.map %>% mutate(Type = replace(Type, Type== "dcl1", "dcl1-2"))
mut.map <- mut.map %>% mutate(Type = replace(Type, Type== "hen1", "hen1-4"))
mut.map <- mut.map %>% mutate(Type = replace(Type, Type== "Control", "Unplanted_soil"))
mut.map <- sample_data(mut.map)
mut.phylo_filtered <- merge_phyloseq(mut.biom_file,mut.map)

#miPEP-community
miPEP.biom_file <- import_biom(here("data","raw", "16S.miPEP.otu_table_filtered.biom"), here("data","raw","16S.miPEP.tree.fasttree"), here("data","raw","16S.miPEP.otus.fasta"), parseFunction=parse_taxonomy_greengenes)
miPEP.map <- data.frame(fread(here("data","raw", "16S.miPEP.mapping_file.tsv"), sep="\t"), check.names=F)
miPEP.map <- sample_data(miPEP.map)
rownames(miPEP.map) <- miPEP.map$`#SampleID`
miPEP.phylo_filtered <- merge_phyloseq(miPEP.biom_file,miPEP.map)

#miPEP-16S qPCR
miPEP.qPCR.16S <- read.table(here("data", "raw", "qPCR-miPEP-16S.txt"), header = T, sep="\t")

#Growth curves miRNA in mix AA
Dosedf3<-read.table(file.path("data","raw", "L-AA_screen3.tsv"), header=T, sep="\t", comment.char = "") # 4452 obs 8 variables
Dosedf3$Time<-as.factor(Dosedf3$Time)
Dosedf3$Pair<-as.factor(Dosedf3$Pair)
Dosedf3$Dose.uM.<-as.factor(Dosedf3$Dose.uM.)
Dosedf3$Treatment<-recode_factor(Dosedf3$Treatment, scrambled_miRNA= "Scramble miRNAs",plant_miRNA="Plant miRNAs") #Renaming variable

#Community composition miRNA in mix AA
Jess_ASV_df<-read.table(file=here("data", "raw", "Jess_feature_table_filtered.tsv"), row.names=1, header=T, sep="\t", comment.char = "", check.names = F) #80 obs 171 variables
Jess_ASV_df2<-Jess_ASV_df[!grepl('Chloroplast', Jess_ASV_df$taxonomy),]#80 obs 171 variables #remove Chloroplasts
Jess_ASV_df3<-Jess_ASV_df2[!grepl('Mitochondria', Jess_ASV_df2$taxonomy),]##80 obs 171 variables variables #remove Mitochondria
Jess_ASV_df4<-Jess_ASV_df3[!grepl('Eukaryota', Jess_ASV_df3$taxonomy),]##76 obs 171 variables variables #remove Eukaryotes
#removing taxonomy column
Jess_ASV_df_notaxo <- Jess_ASV_df4[,-171]
Jess_mapping_df<-read.table(file=here("data", "raw","Jess_mapping_file.tsv"), row.names=1, header=T, sep="\t", comment.char = "", check.names = F) #170 obs. 5 variables 
Jess_mapping_df$Pair<-as.factor(Jess_mapping_df$Pair)
