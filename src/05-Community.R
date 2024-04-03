###Bacterial community analyses
##Mutants
#Remove WT14 samples: these samples are biologically unappropriate for our analysis, as they were not cultivated and sampled in the same conditions as the other samples. Therefore, only 2 wild type samples were appropriate for the following analyses, which was statistically detrimental.
merge_noWT14 <- subset_samples(mut.phylo_filtered, !(Type%in%c("WT14.1","WT14.2")))
merge_noWT14
mut.phylo_filtered

#Extract tax and otu tables
mut.tax_table <- tax_table(merge_noWT14)
#remplacer les noms de colonnes Rank1, rank2.. par les noms de taxons
colnames(mut.tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank8", "rank9", "rank10", "rank11","rank12")
mut.OTU_table <- otu_table(merge_noWT14)

#Count the number of sequences in each sample. 
Seq<-sample_sums(merge_noWT14)
Seq
min(Seq) # gives the lowest number of sequences in a sample = 6 576
max(Seq) # gives the highest number of sequences in a sample = 49 766

#Kept ASVs that are at least present once in two samples (this step is needed, because we took out WT14 samples, some ASVs are only present in 1 or 2 samples). Two samples because, we have two WT roots and two WT soils (in case ASVs are specific to one or the other). Updated taxa and OTU table after.
mut.phylo_filtered <- filter_taxa(merge_noWT14, function(x) sum(x >= 1)>= (2), TRUE)
mut.OTU_table_phylo_filtered <- otu_table(mut.phylo_filtered)
mut.tax_table_phylo_filtered = tax_table(mut.phylo_filtered)
sample_sums(mut.phylo_filtered)

#Order mutant types
correct.order <- c( "ago1-27", "dcl1-2", "hen1-4", "RTL1", "RTL1myc","WT","Unplanted_soil")
sample_data(mut.phylo_filtered)$Type <- factor(sample_data(mut.phylo_filtered)$Type, levels = correct.order)

#Roots/rhizoplane only
mut.phylo_roots <- subset_samples(mut.phylo_filtered, Compartment=="RACINE")

#Plot abundances, top 10 phyla in all samples
mut.pseq_roots <- mut.phylo_roots %>% aggregate_taxa(level = "Phylum")
mut.keep_phyla <- names(sort(taxa_sums(mut.pseq_roots), decreasing = TRUE))[1:10]

# Transform Taxa counts to relative abundance using total number of reads - include others
mut_phylum_relabun <- data.frame(apply(mut.pseq_roots@otu_table, 1, "/", colSums(mut.pseq_roots@otu_table))*100)
rowSums(mut_phylum_relabun) #100 everywhere
mut_phylum_relabun_top10 <- mut_phylum_relabun[,colnames(mut_phylum_relabun)%in%mut.keep_phyla]
others <- 100-rowSums(mut_phylum_relabun_top10)
mut_phylum_relabun_top10 <- data.frame(mut_phylum_relabun_top10, others)
rowSums(mut_phylum_relabun_top10)#Sanity check; 100 all over

#Add mapping file info
mut.map.s <- mut.map[order(row.names(mut.map)),]#Sort mapping file
mut.map.root.s <- mut.map.s[mut.map.s$Compartment=="RACINE",]
mut.map.root.s <- mut.map.root.s[-(20:21),] #Two samples that were removed from the OTU table
mut_phylum_relabun_top10.s <- mut_phylum_relabun_top10[order(row.names(mut_phylum_relabun_top10)),] #Sort Phylum table
row.names(mut_phylum_relabun_top10.s)==row.names(mut.map.root.s)#Sanity check - all TRUE
mut.phy.map.root <- data.frame(mut.map.root.s,mut_phylum_relabun_top10.s)#Merge
rowSums(mut.phy.map.root[,3:13]) #Should be all 100s
mut.phy.map.root.long <- gather(mut.phy.map.root,Phylum,relabund,3:13) #transform in long format for ggplot

mut_StackedBarPlot_phylum_rel <- ggplot(mut.phy.map.root.long, aes(x =Type, y = relabund, fill = Phylum)) +
  geom_bar(stat = "summary", fun ="mean", position = "stack") +
  labs(y = "Relative abundance (%)", x="Genotype") +
  #facet_grid(~ Compartment, scales = "free", labeller = labeller(Compartment = comp.labs)) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_discrete(labels=c("ago1-27", "dcl1-2", "hen1-4", "RTL1", "RTL1myc","WT"))+
  theme_bw()+
  scale_fill_manual(values = c("#264653", "#F4A261", "#E56B6F", "#62929E", "#5B8E7D","#8EA8C3", "#BCB8B1","#30AB65","#FFC857","#BDD9BF","#F9844A"))

mut_StackedBarPlot_phylum_rel

#Anovas
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Acidobacteria"),])) #F=7.805, P=0.000856***
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Actinobacteria"),])) #F=3.208, P=0.363*
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Bacteroidetes"),])) #NS
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Chlorobi"),])) #NS
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Chloroflexi"),])) #F=5.678, P=0.00391**
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Gemmatimonadetes"),])) #NS
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Planctomycetes"),])) #F=5.292, P=0.00534**
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Proteobacteria"),])) #NS
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Spirochaetae"),])) #F=3.558, P=0.0255
summary(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Verrucomicrobia"),])) #F=6.367, P=0.00231**
#Tukey
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Acidobacteria"),])) 
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Actinobacteria"),]))
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Chloroflexi"),])) 
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Planctomycetes"),]))
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Spirochaetae"),])) 
TukeyHSD(aov(relabund~Type, data = mut.phy.map.root.long[(mut.phy.map.root.long$Compartment=="RACINE" & mut.phy.map.root.long$Phylum=="Verrucomicrobia"),]))

#Permanova
row.names(mut.map.root.s) == colnames(mut.phylo_roots@otu_table)#Order is ok
set.seed(22345)
adonis2(data.frame(t(mut.phylo_roots@otu_table))~Type, data = data.frame(mut.map.root.s), method = "bray") #F=3.3386, P=0.001***
#Pairwise adonis
set.seed(9456)
pairwise.adonis2(data.frame(t(mut.phylo_roots@otu_table))~Type, data = data.frame(mut.map.root.s), method = "bray")
#Ago-hen, ago-RTL1, hen-RTL1myc, RTL1-RTL1myc, (WT all nearly... but 2 reps so not significant)

##Alpha diversity
#Rarefaction roots
mut.reads_roots <- sample_sums(mut.phylo_roots)
min(mut.reads_roots) #6563
mut.phylo_rar_roots <- rarefy_even_depth(mut.phylo_roots, sample.size = 6563, rngseed = 711, replace = FALSE)
mut.alpha_roots <- estimate_richness(mut.phylo_rar_roots, measures=c("Shannon"))
row.names(mut.alpha_roots)==row.names(mut.map.root.s)#Already in the right order

# Test normality
shapiro.test(mut.alpha_roots$Shannon) #W = 0.96111, p-value = 0.5388
#Anova
summary(aov(Shannon~mut.map.root.s$Type, data = mut.alpha_roots)) #F=4.246, P=0.0132
TukeyHSD(aov(Shannon~mut.map.root.s$Type, data = mut.alpha_roots)) 
# ago dcl hen RTL RTLmyc  WT
#  a  a   ab  ab  a       b
#Tukey letters - for plot
mut.alpha.tukey <- data.frame(letters = c("a", "a", "ab", "ab", "a", "b"), x = c(1,2,3,4,5,6), y=c(5.4, 5.7, 5.4, 5.4, 5.55, 4.8))


#Plot
mut.alpha.map <- cbind(mut.map.root.s, mut.alpha_roots)
mut.div.plot <- ggplot(mut.alpha.map,aes(x =Type, y = Shannon))+
  labs(x="Genotype", y="Shannon diversity")+
  geom_boxplot() +
  geom_point() +
  geom_text(data = mut.alpha.tukey, aes(x=x, y=y, label = letters))+ #Add letters from tukey
  theme_bw() 
mut.div.plot


###miPEP - only look in roots/rhizoplane for article - data already filtered for ASV present in at least 3 samples
#Extract the OTU and taxa table
miPEP.OTU_table_phylo_filtered <- otu_table(miPEP.phylo_filtered)
miPEP.tax_table_phylo_filtered <- tax_table(miPEP.phylo_filtered)
#Replace colnames by taxonomical rank
colnames(miPEP.tax_table_phylo_filtered) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "rank8", "rank9", "rank10", "rank11","rank12")

#Before continuing, we must discard the sample cR7 from the root samples, which has too little reads.
miPEP.phylo_filtered <- subset_samples(miPEP.phylo_filtered, X.SampleID != "cR7")
miPEP.phylo_filtered #should have only 95 samples now (-cR7)

#Reorder
sample_data(miPEP.phylo_filtered)$Type <- factor(sample_data(miPEP.phylo_filtered)$Type, levels = c("miPEP-A", "miPEP-B", "miPEP-C", "Control-mi-PEP"))

miPEP.phylo_roots <- subset_samples(miPEP.phylo_filtered, Compartment=="RACINE")
miPEP.pseq_roots <- miPEP.phylo_roots %>% aggregate_taxa(level = "Phylum")
miPEP.keep_phyla <- names(sort(taxa_sums(miPEP.pseq_roots), decreasing = TRUE))[1:10]

# Transform Taxa counts to relative abundance using total number of reads - include others
miPEP_phylum_relabun <- data.frame(apply(miPEP.pseq_roots@otu_table, 1, "/", colSums(miPEP.pseq_roots@otu_table))*100)
rowSums(miPEP_phylum_relabun) #100 everywhere
miPEP_phylum_relabun_top10 <- miPEP_phylum_relabun[,colnames(miPEP_phylum_relabun)%in%miPEP.keep_phyla]
others <- 100-rowSums(miPEP_phylum_relabun_top10)
miPEP_phylum_relabun_top10 <- data.frame(miPEP_phylum_relabun_top10, others)
rowSums(miPEP_phylum_relabun_top10)#Sanity check; 100 all over

#Add mapping file info
miPEP.map.s <- miPEP.map[order(row.names(miPEP.map)),]#Sort mapping file
miPEP.map.root.s <- miPEP.map.s[miPEP.map.s$Compartment=="RACINE",]
miPEP.map.root.s <- miPEP.map.root.s[-23,] #One sample that was removed from the OTU table
miPEP_phylum_relabun_top10.s <- miPEP_phylum_relabun_top10[order(row.names(miPEP_phylum_relabun_top10)),] #Sort Phylum table
row.names(miPEP_phylum_relabun_top10.s)==row.names(miPEP.map.root.s)#Sanity check - all TRUE
miPEP.phy.map.root <- data.frame(miPEP.map.root.s,miPEP_phylum_relabun_top10.s)#Merge
rowSums(miPEP.phy.map.root[,5:15]) #Should be all 100s
miPEP.phy.map.root.long <- gather(miPEP.phy.map.root,Phylum,relabund,5:15) #transform in long format for ggplot

miPEP_StackedBarPlot_phylum_rel <- ggplot(miPEP.phy.map.root.long, aes(x =Type, y = relabund, fill = Phylum)) +
  geom_bar(stat = "summary", fun ="mean", position = "stack") +
  labs(y = "Relative abundance (%)", x="Treatment") +
  #facet_grid(~ Compartment, scales = "free", labeller = labeller(Compartment = comp.labs)) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_x_discrete(labels=c("Water", "miPEP-A", "miPEP-B", "miPEP-C"))+
  theme_bw()+
  scale_fill_manual(values = c("#264653", "#F4A261", "#E56B6F","#5B8E7D","#04395E","#25A18E","#8EA8C3" ,"#BCB8B1","#30AB65","#FFC857","#F9844A"))

miPEP_StackedBarPlot_phylum_rel

#Anovas
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Acidobacteria"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Actinobacteria"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Bacteroidetes"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Chloroflexi"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Cyanobacteria"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="FCPU426"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Gemmatimonadetes"),])) #NS
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Planctomycetes"),])) #40.78 3.64e-10 ***
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Proteobacteria"),])) #8.91 0.000287 ***
summary(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Verrucomicrobia"),])) #NS
#Tukey
TukeyHSD(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Planctomycetes"),]))
TukeyHSD(aov(relabund~Type, data = miPEP.phy.map.root.long[(miPEP.phy.map.root.long$Compartment=="RACINE" & miPEP.phy.map.root.long$Phylum=="Proteobacteria"),]))

#Permanova
row.names(miPEP.map.root.s) == colnames(miPEP.phylo_roots@otu_table)#Order is ok
set.seed(22345)
adonis2(data.frame(t(miPEP.phylo_roots@otu_table))~Type, data = data.frame(miPEP.map.root.s), method = "bray") #R2=0.16026, P=0.031
#Pairwise adonis
set.seed(9456)
pairwise.adonis2(data.frame(t(miPEP.phylo_roots@otu_table))~Type, data = data.frame(miPEP.map.root.s), method = "bray")

##Alpha diversity
#Rarefaction roots
miPEP.reads_roots <- sample_sums(miPEP.phylo_roots)
min(miPEP.reads_roots) #9385
miPEP.phylo_rar_roots <- rarefy_even_depth(miPEP.phylo_roots, sample.size = 9385, rngseed = 711, replace = FALSE)
miPEP.alpha_roots <- estimate_richness(miPEP.phylo_rar_roots, measures=c("Shannon"))
row.names(miPEP.alpha_roots)==row.names(miPEP.map.root.s)#Already in the right order

# Test normality
shapiro.test(miPEP.alpha_roots$Shannon) #W = 0.95816, p-value = 0.2606
#Anova
summary(aov(Shannon~miPEP.map.root.s$Type, data = miPEP.alpha_roots)) #F=0.953, P=0.429

#qPCR
miPEP.qPCR.16S.s <- miPEP.qPCR.16S[order(miPEP.qPCR.16S$Sample),]
miPEP.qPCR.16S.s$Sample==miPEP.map.s$X.SampleID #Sort check: 96 TRUE
miPEP.qPCR.16S.sample <- cbind(miPEP.map.s,miPEP.qPCR.16S.s)

#Test normality
shapiro.test(miPEP.qPCR.16S.sample[miPEP.qPCR.16S.sample$Compartment=="RACINE",]$X16S) #W = 0.9264, p-value = 0.03112
#Anova
summary(aov(X16S~Type, data = miPEP.qPCR.16S.sample[miPEP.qPCR.16S.sample$Compartment=="RACINE",]))#F=1.94, P=0.146

###Growth curves miRNA in AA mix
Mix17<-Dosedf3[Dosedf3$N_source=="Mix17",] #Keep only mix of amino acid for this publication
Mix17<-Mix17[Mix17$miRNA=="mix_5",] #Keep only all miRNAs for this publication
Mix17_mean<-Mix17 %>%
  group_by(Time,Treatment) %>%
  summarise(Mean = mean(OD_600),sd=sd(OD_600), se = sd / sqrt(5))

#Plot
boxMix17_mean<- ggplot(Mix17_mean, aes(x=Time, y=Mean ,fill=Treatment,group=Treatment, shape=Treatment, colour=Treatment)) +
  geom_errorbar(alpha=0.4,aes(ymin=Mean-se, ymax=Mean+se), width=.3, position=position_dodge(0.08))+
  geom_line(alpha=0.6,aes(linetype=Treatment))+
  geom_point(aes(shape=Treatment),alpha=0.6,)+
  theme_classic() + 
  scale_fill_manual(values = c("#003049","#2a7f62"))+
  scale_colour_manual(values = c("#003049","#2a7f62"))+
  labs(x="", y="")+ ylim(0,2.0)+theme_pubclean()+
  theme(legend.title = element_text(size = 24)) +
  theme(legend.title = element_text(size = 12) ,legend.text = element_text(size = 12))+
  theme(axis.text.x=element_text(size=6))+
  annotate("rect",xmin = 7.5, xmax=9.5, ymin = 0, ymax=2.0,alpha=.2, fill="steelblue4")+
  labs(x="Time (hours)", y="Optical density (600nm)")

boxMix17_mean

###Community miRNA in AA mix
#Pre- stacked bar charts steps
##Separating taxonomic groups 
taxonomy_16s<-as.data.frame(Jess_ASV_df4[,171,drop = FALSE])#6484 obs of 1 variable
colnames(taxonomy_16s) <- c("taxonomy")
tax_sep <- separate(taxonomy_16s, taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
## Eliminating taxon level, characters and spacers
tax_sep$Domain <- stringr::str_replace(tax_sep$Domain,'[k]', '')
tax_sep$Domain <- stringr::str_replace_all(tax_sep$Domain,'[__]', '')
tax_sep$Phylum <- stringr::str_replace(tax_sep$Phylum,'[p__]', '')
tax_sep$Phylum <- stringr::str_replace_all(tax_sep$Phylum,'[__]', '')
tax_sep$Class <- stringr::str_replace(tax_sep$Class,'[c__]', '')
tax_sep$Class <- stringr::str_replace_all(tax_sep$Class,'[__]', '')
tax_sep$Order <- stringr::str_replace(tax_sep$Order, '[o]', '')
tax_sep$Order <- stringr::str_replace_all(tax_sep$Order, '[__]', '')
tax_sep$Family <- stringr::str_replace(tax_sep$Family, '[f]', '')
tax_sep$Family <- stringr::str_replace_all(tax_sep$Family, '[__]', '')
tax_sep$Genus <- stringr::str_replace(tax_sep$Genus, '[g]','')
tax_sep$Genus <- stringr::str_replace_all(tax_sep$Genus, '[__]', '')

## Merging ASVs with clean taxonomy
Jess_ASV_taxo<- cbind(Jess_ASV_df4[,1:170], tax_sep)
Jess_ASV_taxo$'Family and Genus' <- paste(Jess_ASV_taxo$Family,Jess_ASV_taxo$Genus)
Just_Tax<-as.data.frame(Jess_ASV_taxo[,177,drop = FALSE])
rowstax<-row.names(Just_Tax)
Just_Tax$'ASVs'<-rowstax
rowstax==Just_Tax$'ASVs' #Sanity check
Just_Tax["Others",]<-rbind("Others") #Add a row with Others as a ASV and Taxonomy

## Preparing abundance matrix
tax_clean_16s <- Jess_ASV_taxo # working file for taxonomy
com16s_notaxo <- tax_clean_16s[,1:170] # creating an asv table without taxonomy column
tcom16s_notaxo <- t(com16s_notaxo) # transposing table, samples as rows, species as columns to match metadata
com16_rel <- (tcom16s_notaxo/rowSums(tcom16s_notaxo)) #Normalization by relative abundance
(rowSums(com16_rel)) #Sanity check

## Shaping data frame for stack bar charts
# Order mapping file
dim(Jess_mapping_df) #170 5
mapping16_sorted <-Jess_mapping_df[order(row.names(Jess_mapping_df)),] # order mappingfile
rownames(mapping16_sorted) # Check
rowsmap<-row.names(mapping16_sorted)
mapping16_sorted$'row.names'<-rowsmap
row.names(mapping16_sorted)==mapping16_sorted$'row.names' #Sanity check

# Order ASV table
com16_rel_sorted <- com16_rel[order(row.names(com16_rel)),]
dim(com16_rel_sorted)# 170  76
rownames(com16_rel_sorted)
com16_abund <- com16_rel_sorted[,colMeans(com16_rel_sorted) > 0.01]
dim(com16_abund) # 170   7
Others<-as.data.frame(1 - rowSums(com16_abund))
colnames(Others)<-"Others"
com16_abund_all<-cbind(com16_abund,Others)
rowscom<-row.names(com16_abund_all)
com16_abund_all$'ID'<-rowscom

com16_abund_melt<-reshape2::melt(com16_abund_all)
mapcom16<- cbind(com16_abund_melt, mapping16_sorted)
mapcom16_tax<- merge(mapcom16, Just_Tax, by.x ="variable", by.y = "ASVs")

##Fixing legend order and factor names
mapcom16_tax$`Family and Genus` <- factor(mapcom16_tax$`Family and Genus`, levels = c("Enterobacteriaceae Citrobacter","Enterobacteriaceae Enterobacter","Enterobacteriaceae Raoultella", "Moraxellaceae Acinetobacter","Pseudomonadaceae Pseudomonas","Others"))

##Subsetting data frame by Compartment
all5_mapcom16<-mapcom16_tax[mapcom16_tax$Dose=="10",]
all5_mapcom16_mixAA<-all5_mapcom16[all5_mapcom16$N_Source=="mix_17AA",]

#Stack bar charts: mean of relative abundances
#Calculating mean abundances 
mapcom16_mean_tax_mixAA<-all5_mapcom16_mixAA %>%
  group_by(N_Source, Treatment,variable,`Family and Genus`) %>%
  summarise(Mean = mean(value))

#Stack bar with mean abundances
stack_mean_mixAA <-ggplot(mapcom16_mean_tax_mixAA, aes(fill =`Family and Genus`, y = Mean, x=Treatment)) +
  geom_bar(stat = "identity",position = "fill") +
  ylab("Relative abundance") + 
  xlab("")+
  theme_bw() +
  scale_fill_manual(values =c("#ffc43d","#f95d6a","#d45087","#2f4b7c","#457b9d","lightgrey"), guide = guide_legend(label.theme = element_text(face = "italic", size = 11))) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_x_discrete()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),strip.background = element_rect(fill="white", linewidth = 0, linetype="solid"))+theme(plot.title = element_text(hjust = 0.5))
stack_mean_mixAA

###ASV significantly affected 
## Preparing abundance matrix
ASV_df_notaxo_t <- t(Jess_ASV_df_notaxo) # transposing table, samples as rows, species as columns to match metadata
#Merging
row.names(mapping16_sorted)==row.names(com16_rel)#Check, should be all true
com16_rel_map<- merge(mapping16_sorted[,-6],com16_rel, by = 'row.names', all = TRUE )

#Exposure to all five miRNAs
Dose10uM<-com16_rel_map[com16_rel_map$Dose=="10",]
Mix10<-Dose10uM[Dose10uM$N_Source=="mix_17AA",]

#Subsetting by ASVs of interest
ASVs_interest<- c("Row.names","N_Source","miRNA","Dose","Treatment","Pair","5","1","7","9")
ASVs_interest_df<- Mix10[,colnames(Mix10) %in% ASVs_interest ] #10 obs 10 variables

ASV5Mix<-ASVs_interest_df[(ASVs_interest_df$N_Source=="mix_17AA"), c("Dose","Treatment","N_Source","Pair","5")]
colnames(ASV5Mix)[5] <- "Relative_abundance"  
ASV5Mix$Taxonomy<- "ASV5 \n Acinetobacter"

ASV1Mix<-ASVs_interest_df[(ASVs_interest_df$N_Source=="mix_17AA"), c("Dose","Treatment","N_Source","Pair","1")]
colnames(ASV1Mix)[5] <- "Relative_abundance"  
ASV1Mix$Taxonomy<- "ASV1 \n Enterobacter"

ASV7Mix<-ASVs_interest_df[(ASVs_interest_df$N_Source=="mix_17AA"), c("Dose","Treatment","N_Source","Pair","7")]
colnames(ASV7Mix)[5] <- "Relative_abundance"  
ASV7Mix$Taxonomy<- "ASV7 \n Citrobacter"

ASV9Mix<-ASVs_interest_df[(ASVs_interest_df$N_Source=="mix_17AA"), c("Dose","Treatment","N_Source","Pair","9")]
colnames(ASV9Mix)[5] <- "Relative_abundance"  
ASV9Mix$Taxonomy<- "ASV9 \n Enterobact."

combined_ASVs <- rbind(ASV5Mix,ASV1Mix,ASV7Mix,ASV9Mix)
combined_ASVs$N_Taxonomy<-paste(combined_ASVs$N_Source,combined_ASVs$Taxonomy)
combined_ASVs$N_Source[combined_ASVs$N_Source == 'mix_17AA'] <- 'Mix of 17 amino acids'

boxall<- ggplot(combined_ASVs, aes(x=Treatment, y=Relative_abundance)) +
  geom_boxplot(outlier.shape = NA, alpha=0.4, aes(fill=Treatment, color=Treatment )) +
  geom_point(aes(fill=Treatment, color=Treatment),pch = 21, position = position_jitterdodge())+
  theme_light() +
  facet_wrap2(vars(Taxonomy), scales="free_y", nrow=1)+
  scale_x_discrete()+
  theme(axis.text.x=element_blank() ,axis.ticks.x = element_blank(),strip.background = element_rect(fill="white", size=0, linetype="solid"))+
  theme(strip.text = element_text(colour = 'grey30', size=10,face = "italic" ))+
  scale_fill_manual(values = c("#2a7f62","#003049"))+
  scale_colour_manual(values = c("#2a7f62","#003049"))+
  labs(x="", y="Relative abundance")+
  theme(legend.title = element_text(size = 14) ,legend.text = element_text(size=12))


#Verifying assmumptions
#Normality assumption
combined_ASVs %>%
  group_by(Taxonomy) %>%
  shapiro_test(Relative_abundance) #ok

#Equality of variances assumption
combined_ASVs %>%
  group_by(Taxonomy) %>%
  levene_test(Relative_abundance~Treatment) #ok

#Computing the statistical test
stat.test <- combined_ASVs %>%
  group_by(Taxonomy) %>%
  t_test(Relative_abundance ~ Treatment) %>%
  adjust_pvalue(method = "holm") %>%
  add_significance()
stat.test

#second version
#Adding the coordinates of the p-values
stat.test.pos<-stat.test %>% add_xy_position(x="Treatment", group="Taxonomy",fun = "max")
box_sig2<- boxall+stat_pvalue_manual(
  stat.test.pos, hide.ns = TRUE, 
  label = "p.adj.signif",tip.length = 0, size=6)
box_sig2
