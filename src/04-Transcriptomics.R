###Transcriptomic analyses
##Variovorax
# Convert the column from characters to factors (required by DESeq)
vario.info.2 <- cbind(vario.info,vario.info) #Order does not seem to work on a single column df? Returns a vector...
vario.info.2 <- vario.info.2[order(rownames(vario.info.2)),] 
vario.transcripto <- vario.transcripto[,order(colnames(vario.transcripto))]
rownames(vario.info.2) <- gsub('---', '...', rownames(vario.info.2 ))#Not same sample names
row.names(vario.info.2)==colnames(vario.transcripto)#All true
vario.info.2$Type <- gsub('-', '', vario.info.2$Type)#Suggested by DESeq
vario.info.2$Type <- factor(vario.info.2$Type)#Suggested by DESeq

# Set up the dds object
sample_dds.vario <- DESeqDataSetFromMatrix(countData = vario.transcripto,
                                     colData = vario.info.2,
                                     design = ~Type
)
# Run DESeq (this line of code is short, but it automatically performs all the differential expression analysis)
sample_dds.vario <- DESeq(sample_dds.vario)

# Check result. Differential expression is comparing difference between two conditions.
# Comparing difference between Type A and Type B (as defined in InfoTable)

#20 min
sample_res.20 <- results(sample_dds.vario, contrast = c("Type", "mir20", "scrb20")) #last charact vect = base level
summary(sample_res.20)
# Save the results to a table
res.20 <- data.frame(sample_res.20)
plotMA(sample_res.20, ylim=c(-6,6), colSig="#d01c8b", colLine="black", cex=1)
vario.annotations$gene_id <- gsub("_1","", vario.annotations$gene_id)#Fix extra character in annotation gene_id
res_gene.20 <- merge(x=res.20, y=vario.annotations, by.x=0, by.y="gene_id")
write.table(res.20, file = here("output", "tables", "DESeq-results-20min-vario.txt"))
# Filtering to find significant genes using FDR cutoff of 0.05
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results.20 <- res_gene.20[which(res_gene.20$padj < padj.cutoff),]
# save results using customized file_name
write.table(significant_results.20, file = here("output", "tables", "sign-pajd-20min-vario.txt"))

#120min
sample_res.120 <- results(sample_dds.vario, contrast = c("Type", "mir120", "scrb120")) 
summary(sample_res.120)
plotMA(sample_res.120, ylim=c(-6,6), colSig="#d01c8b", colLine="black", cex=1)
# Save the results to a table
res.120 <- data.frame(sample_res.120)
res_gene.120 <- merge(x=res.120, y=vario.annotations, by.x=0, by.y="gene_id")
write.table(res.120, file = here("output", "tables", "DESeq-results-120min-vario.txt"))
# Filtering to find significant genes using FDR cutoff of 0.05
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results.120 <- res_gene.120[which(res_gene.120$padj < padj.cutoff),]
# save results using customized file_name
write.table(significant_results.120, file = here("output", "tables", "sign-pajd-120min-vario.txt"))

#Look at overlap
intersect(significant_results.20$Row.names, significant_results.120$Row.names) #one gene: gene-VARPA_RS01000

#Volcano plots
vario.volcano.20 <- EnhancedVolcano(res.20, lab = row.names(res.20), selectLab = "",
                x = 'log2FoldChange', y = 'padj',
                pCutoff = 0.05, FCcutoff = 1.5, 
                pointSize =3, labSize = 4, legendPosition = "bottom",
                title = "", subtitle = "", ylim = c( -log10(6.5)))
vario.volcano.20

vario.volcano.120 <- EnhancedVolcano(res.120, lab = row.names(res.120), selectLab = "",
                                    x = 'log2FoldChange', y = 'padj',
                                    pCutoff = 0.05, FCcutoff = 1.5, 
                                    pointSize =3, labSize = 4, legendPosition = "bottom",
                                    title = "", subtitle = "", ylim = c( -log10(4)))
vario.volcano.120


##Bacillus
# Convert the column from characters to factors (required by DESeq)
bacillus.info.2 <- cbind(bacillus.info,bacillus.info) #Order does not seem to work on a single column df? Returns a vector...
bacillus.info.2 <- bacillus.info.2[order(rownames(bacillus.info.2)),] 
bacillus.transcripto <- bacillus.transcripto[,order(colnames(bacillus.transcripto))]
rownames(bacillus.info.2) <- gsub('---', '...', rownames(bacillus.info.2 ))#Not same sample names
row.names(bacillus.info.2)==colnames(bacillus.transcripto)#All true
bacillus.info.2$Type <- gsub('-', '', bacillus.info.2$Type)#Suggested by DESeq
bacillus.info.2$Type <- factor(bacillus.info.2$Type)#Suggested by DESeq

# Set up the dds object
sample_dds.bacillus <- DESeqDataSetFromMatrix(countData = bacillus.transcripto,
                                     colData = bacillus.info.2,
                                     design = ~Type
)
# Run DESeq (this line of code is short, but it automatically performs all the differential expression analysis)
sample_dds.bacillus <- DESeq(sample_dds.bacillus)
# Check result. Differential expression is comparing difference between two conditions.
#20 min
sample_res.20.bacillus <- results(sample_dds.bacillus, contrast = c("Type", "mir20", "scrb20")) #last charact vect = base level
summary(sample_res.20.bacillus) #No significantly diff genes at padj<0.1
#120min
sample_res.120.bacillus <- results(sample_dds.bacillus, contrast = c("Type", "mir120", "scrb120")) 
summary(sample_res.120.bacillus)#No significantly diff genes at padj<0.1

###miPEP and qPCR
#Identify outliers
#pri.miR159c
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "miPEP",], variable="pri.miR159c") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "scrambled",], variable="pri.miR159c") #4.176437, row 18
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "water",], variable="pri.miR159c")# None

#Phosphatidate.cytydylytransferase
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "miPEP",], variable="Phosphatidate.cytidylyltransferase") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "scrambled",], variable="Phosphatidate.cytidylyltransferase") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "water",], variable="Phosphatidate.cytidylyltransferase")# None

#alpha.2.macroglobulin
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "miPEP",], variable="alpha.2.macroglobulin") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "scrambled",], variable="alpha.2.macroglobulin") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "water",], variable="alpha.2.macroglobulin")# None

#LysR
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "miPEP",], variable="LysR") #0.646076, row 1
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "scrambled",], variable="LysR") #None
identify_outliers(vario.qpcr.miPEP[vario.qpcr.miPEP$Treatment == "water",], variable="LysR")# None

#Remove outliers
vario.qpcr.miPEP[18,2] <- NA
vario.qpcr.miPEP[1,5] <- NA

#Make long for plotting
qpcr.long<- gather(vario.qpcr.miPEP, Gene, Ratio, 2:5)
#Order the different genes
qpcr.long$Gene <- factor(qpcr.long$Gene, c("pri.miR159c", "Phosphatidate.cytidylyltransferase", "alpha.2.macroglobulin", "LysR"))
qpcr.long$Treatment <- factor(qpcr.long$Treatment, c("water", "miPEP", "scrambled"))

#To have correct label names
gene.labs <- c("pri-miR159c", "Alpha-2-macroglobulin","CdsA", "LysR")
names(gene.labs) <- c("pri.miR159c",  "alpha.2.macroglobulin", "Phosphatidate.cytidylyltransferase","LysR")

qpcr.plot <- ggplot(qpcr.long[qpcr.long$Treatment!="scrambled",],aes(x =Treatment, y = Ratio))+
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(strip.text = element_text(size = 8))+
  ylab("Relative expression ratio")+
  facet_wrap(~Gene, nrow = 2, ncol =2, scales = "free_y", labeller = labeller(Gene = gene.labs))

qpcr.plot

###### Statistical analysis ######
shapiro.test(qpcr.long[qpcr.long$Gene=="pri.miR159c",3]) #P=0.2269
shapiro.test(log(qpcr.long[qpcr.long$Gene=="LysR",3])) #P=0.1229
shapiro.test(qpcr.long[qpcr.long$Gene=="alpha.2.macroglobulin",3]) #P=0.9153
shapiro.test(log(qpcr.long[qpcr.long$Gene=="Phosphatidate.cytidylyltransferase",3])) #P=0.7232

#T-test between water and miPEP
t.test(qpcr.long[qpcr.long$Gene=="pri.miR159c" & qpcr.long$Treatment=="water",3], 
       qpcr.long[qpcr.long$Gene=="pri.miR159c" & qpcr.long$Treatment=="miPEP",3]) #t=-3.30551, P=0.0145
t.test(log(qpcr.long[qpcr.long$Gene=="Phosphatidate.cytidylyltransferase" & qpcr.long$Treatment=="water",3]), 
       log(qpcr.long[qpcr.long$Gene=="Phosphatidate.cytidylyltransferase" & qpcr.long$Treatment=="miPEP",3])) #t=3.4297, P=0.006474
t.test(qpcr.long[qpcr.long$Gene=="alpha.2.macroglobulin" & qpcr.long$Treatment=="water",3], 
       qpcr.long[qpcr.long$Gene=="alpha.2.macroglobulin" & qpcr.long$Treatment=="miPEP",3]) #t=-2.92, P=0.04621
t.test(log(qpcr.long[qpcr.long$Gene=="LysR" & qpcr.long$Treatment=="water",3]), 
       log(qpcr.long[qpcr.long$Gene=="LysR" & qpcr.long$Treatment=="miPEP",3])) #t=4.2681, P=0.005365

#T-test between water and scrambled
t.test(qpcr.long[qpcr.long$Gene=="pri.miR159c" & qpcr.long$Treatment=="water",3], 
       qpcr.long[qpcr.long$Gene=="pri.miR159c" & qpcr.long$Treatment=="scrambled",3]) #t=-0.37597, P=0.7167
t.test(log(qpcr.long[qpcr.long$Gene=="Phosphatidate.cytidylyltransferase" & qpcr.long$Treatment=="water",3]), 
       log(qpcr.long[qpcr.long$Gene=="Phosphatidate.cytidylyltransferase" & qpcr.long$Treatment=="scrambled",3])) #t=0.66799, P=0.5246
t.test(qpcr.long[qpcr.long$Gene=="alpha.2.macroglobulin" & qpcr.long$Treatment=="water",3], 
       qpcr.long[qpcr.long$Gene=="alpha.2.macroglobulin" & qpcr.long$Treatment=="scrambled",3]) #t=0.17628, P=0.8634
t.test(log(qpcr.long[qpcr.long$Gene=="LysR" & qpcr.long$Treatment=="water",3]), 
       log(qpcr.long[qpcr.long$Gene=="LysR" & qpcr.long$Treatment=="scrambled",3])) #t=0.17075, P=0.8676
