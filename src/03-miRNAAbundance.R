###Arrange miRNA abundance tables, keep only abundant miRNA with >10 reads in all reps and absent in unplanted (1 read or less)

##Arabidopsis-FR-Rhizosphere
FR_ara_prelim <- FR_ara[,colnames(FR_ara)%in%row.names(FR_ara_map[FR_ara_map$Other=="Preliminary.exp",])]
FR_ara_prelim_colsums <- colSums(FR_ara_prelim)#calculate colsums for normalization later
FR_ara_prelim <- FR_ara_prelim[rowSums(FR_ara_prelim)!=0,] #Remove zero only lines
FR_ara_prelim <- FR_ara_prelim[order(row.names(FR_ara_prelim)),]#Order
miRNA.df <- miRNA.df[order(miRNA.df$Identity),]#Order
miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_prelim),]$Identity == row.names(FR_ara_prelim)#Check - all True
FR_ara_prelim$Sequence <- miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_prelim), "Sequence"] #Add sequence to table
FR_ara_prelim$miRNA <- miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_prelim), "miRNA"] #Add sequence to table
FR_ara_prelim <- FR_ara_prelim %>% #Concatenate identical mature miRNAs
  group_by(Sequence) %>%
  mutate(miRNAs = paste0(miRNA, collapse = ";")) %>%
  group_by(miRNAs,Sequence)%>%
  summarise(across(1:6, sum))
colSums(FR_ara_prelim[,3:8])==FR_ara_prelim_colsums#Check
FR_ara_prelim <- FR_ara_prelim[(FR_ara_prelim[,5]>10 & FR_ara_prelim[,6]>10 & FR_ara_prelim[,7]>10), ] #Keep only >10 reads
FR_ara_prelim_norm <- data.frame(apply(FR_ara_prelim[,3:8], 1, "/", FR_ara_prelim_colsums))#normalize
rowSums(FR_ara_prelim_norm)#sanity check
colnames(FR_ara_prelim_norm) <- FR_ara_prelim$miRNAs

#Prepare for plot
FR_ara_prelim_norm <- FR_ara_prelim_norm[3:5,]
row.names(FR_ara_prelim_norm) <- c("Rep1", "Rep2", "Rep3")
FR_ara_prelim_long <- gather(FR_ara_prelim_norm)#long format for plotting
colnames(FR_ara_prelim_long) <- c("miRNA", "RelativeAbundance")
FR_ara_prelim_long <- FR_ara_prelim_long[order(FR_ara_prelim_long$miRNA),]
FR_ara_prelim_long$miRNA <- gsub("miR162a-3p;miR162b-3p", "miR162a,b", FR_ara_prelim_long$miRNA)
FR_ara_prelim_long$miRNA <- gsub("miR166a-3p;miR166b-3p;miR166c;miR166d;miR166e-3p;miR166f;miR166g", "miR166a,b,c,d,e,f,g", FR_ara_prelim_long$miRNA)
FR_ara_prelim_long$miRNA <- gsub("miR165a-3p;miR165b", "miR165a,b", FR_ara_prelim_long$miRNA)

#Plot
ara_miRNA_box <- ggplot(FR_ara_prelim_long, aes(x=reorder(miRNA, -RelativeAbundance), y=100*RelativeAbundance))+
              geom_boxplot()+
              theme_bw()+
              ylab("Proportion of Arabidopsis reads (%)")+
              xlab("")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.y =  element_text(size=8))
ara_miRNA_box

##Brachypodium-FR-Rhizosphere
FR_bra_prelim <- FR_bra[,colnames(FR_bra)%in%row.names(FR_bra_map[FR_bra_map$Experiment=="Preliminary.exp",])]
FR_bra_prelim_colsums <- colSums(FR_bra_prelim)#calculate colsums for normalization later
FR_bra_prelim <- FR_bra_prelim[rowSums(FR_bra_prelim)!=0,] #Remove zero only lines
FR_bra_prelim <- FR_bra_prelim[order(row.names(FR_bra_prelim)),]#Order
miRNA.df <- miRNA.df[order(miRNA.df$Identity),]#Order
miRNA.df[miRNA.df$Identity %in% row.names(FR_bra_prelim),]$Identity == row.names(FR_bra_prelim)#Check - all True
FR_bra_prelim$Sequence <- miRNA.df[miRNA.df$Identity %in% row.names(FR_bra_prelim), "Sequence"] #Add sequence to table
FR_bra_prelim$miRNA <- miRNA.df[miRNA.df$Identity %in% row.names(FR_bra_prelim), "miRNA"] #Add sequence to table
FR_bra_prelim <- FR_bra_prelim %>% #Concatenate identical mature miRNAs
  group_by(Sequence) %>%
  mutate(miRNAs = paste0(miRNA, collapse = ";")) %>%
  group_by(miRNAs,Sequence)%>%
  summarise(across(1:5, sum))
colSums(FR_bra_prelim[,3:7])==FR_bra_prelim_colsums#Check
FR_bra_prelim <- FR_bra_prelim[(FR_bra_prelim[,3]>10 & FR_bra_prelim[,4]>10), ] #Keep only >10 reads
FR_bra_prelim_norm <- data.frame(apply(FR_bra_prelim[,3:7], 1, "/", FR_bra_prelim_colsums))#normalize
rowSums(FR_bra_prelim_norm)#sanity check
colnames(FR_bra_prelim_norm) <- FR_bra_prelim$miRNAs

#Prepare for plot
FR_bra_prelim_norm <- FR_bra_prelim_norm[1:2,]
row.names(FR_bra_prelim_norm) <- c("Rep1", "Rep3")
FR_bra_prelim_long <- gather(FR_bra_prelim_norm)#long format for plotting
colnames(FR_bra_prelim_long) <- c("miRNA", "RelativeAbundance")
FR_bra_prelim_long$miRNA <- gsub("miR156b-5p;miR156c;miR156d-5p;miR156e-5p;miR156f-5p;miR156g-5p;miR156h-5p;miR156i-5p", "miR156b,c,d,e,f,g,h,i", FR_bra_prelim_long$miRNA)
FR_bra_prelim_long$miRNA <- gsub("miR166a-3p;miR166b-3p;miR166c-3p;miR166d-3p;miR166i-3p", "miR166a,b,c,d,i", FR_bra_prelim_long$miRNA)
FR_bra_prelim_long$miRNA <- gsub("miR167a;miR167b;miR167f", "miR167a,b,f", FR_bra_prelim_long$miRNA)
FR_bra_prelim_long$miRNA <- gsub("miR167c-5p;miR167d-5p;miR167e-5p;miR167g", "miR167c,d,e,g", FR_bra_prelim_long$miRNA)
FR_bra_prelim_long$miRNA <- gsub("miR396a-5p;miR396b-5p", "miR396a,b", FR_bra_prelim_long$miRNA)
FR_bra_prelim_long <- FR_bra_prelim_long[order(FR_bra_prelim_long$miRNA),]

#Plot
bra_miRNA_box <- ggplot(FR_bra_prelim_long, aes(x=reorder(miRNA, -RelativeAbundance), y=100*RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Proportion of Brachypodium reads (%)")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.y =  element_text(size=8))
bra_miRNA_box

#Look for overlap
intersect(colnames(FR_ara_prelim_norm),colnames(FR_bra_prelim_norm))

##Arabidospis root Canada
CAN_ara_nitro <- CAN_ara[,colnames(CAN_ara)%in%row.names(CAN_map[(CAN_map$Compartment=="Endosphere" & CAN_map$Nitrogen=="No.Added.Nitrogen" & CAN_map$Plant=="Arabidopsis"),])]
CAN_ara_nitro_colsums <- colSums(CAN_ara_nitro)#calculate colsums for normalization later
CAN_ara_nitro <- CAN_ara_nitro[rowSums(CAN_ara_nitro)!=0,] #Remove zero only lines
CAN_ara_nitro <- CAN_ara_nitro[order(row.names(CAN_ara_nitro)),]#Order
miRNA.df <- miRNA.df[order(miRNA.df$Identity),]#Order
miRNA.df[miRNA.df$Identity %in% row.names(CAN_ara_nitro),]$Identity == row.names(CAN_ara_nitro)#Check - all True
CAN_ara_nitro$Sequence <- miRNA.df[miRNA.df$Identity %in% row.names(CAN_ara_nitro), "Sequence"] #Add sequence to table
CAN_ara_nitro$miRNA <- miRNA.df[miRNA.df$Identity %in% row.names(CAN_ara_nitro), "miRNA"] #Add sequence to table
CAN_ara_nitro <- CAN_ara_nitro %>% #Concatenate identical mature miRNAs
  group_by(Sequence) %>%
  mutate(miRNAs = paste0(miRNA, collapse = ";")) %>%
  group_by(miRNAs,Sequence)%>%
  summarise(across(1:5, sum))
colSums(CAN_ara_nitro[,3:7])==CAN_ara_nitro_colsums#Check
CAN_ara_nitro <- CAN_ara_nitro[(CAN_ara_nitro[,3]>10 & CAN_ara_nitro[,4]>10 & CAN_ara_nitro[,5]>10& CAN_ara_nitro[,6]>10& CAN_ara_nitro[,7]>10), ] #Keep only >10 reads
CAN_ara_nitro_norm <- data.frame(apply(CAN_ara_nitro[,3:7], 1, "/", CAN_ara_nitro_colsums))#normalize
rowSums(CAN_ara_nitro_norm)#sanity check
colnames(CAN_ara_nitro_norm) <- CAN_ara_nitro$miRNAs

#Prepare for plot
CAN_ara_nitro_long <- gather(CAN_ara_nitro_norm)#long format for plotting
colnames(CAN_ara_nitro_long) <- c("miRNA", "RelativeAbundance")
CAN_ara_nitro_long$miRNA <- gsub("miR160a-5p;miR160b;miR160c-5p", "miR160a,b,c", CAN_ara_nitro_long$miRNA)
CAN_ara_nitro_long$miRNA <- gsub("miR162a-3p;miR162b-3p", "miR162a,b", CAN_ara_nitro_long$miRNA)
CAN_ara_nitro_long$miRNA <- gsub("miR165a-3p;miR165b", "miR165a,b", CAN_ara_nitro_long$miRNA)
CAN_ara_nitro_long$miRNA <- gsub("miR166a-3p;miR166b-3p;miR166c;miR166d;miR166e-3p;miR166f;miR166g", "miR166a,b,c,d,e,f,g", CAN_ara_nitro_long$miRNA)
CAN_ara_nitro_long$miRNA <- gsub("miR168a-5p;miR168b-5p", "miR168a,b", CAN_ara_nitro_long$miRNA)
CAN_ara_nitro_long <- CAN_ara_nitro_long[order(CAN_ara_nitro_long$miRNA),]

#Plot 
root_ara_miRNA_box <- ggplot(CAN_ara_nitro_long, aes(x=reorder(miRNA, -RelativeAbundance), y=100*RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Proportion of Arabidopsis reads (%)")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.title.y =  element_text(size=8))
root_ara_miRNA_box

#Look for overlap
intersect(colnames(FR_ara_prelim_norm), colnames(CAN_ara_nitro_norm))

##Bacteria from soil and rhizosphere - France, use only 020320 - other dates have some weird samples
FR_ara_bact <- FR_ara[,colnames(FR_ara)%in%row.names(FR_ara_map[FR_ara_map$Date=="20320",])]
FR_ara_bact <- FR_ara_bact[rowSums(FR_ara_bact)!=0,] #Remove zero only lines
FR_ara_bact_colsums <- colSums(FR_ara_bact)#calculate colsums for normalization later
FR_ara_bact <- FR_ara_bact[order(row.names(FR_ara_bact)),]#Order
miRNA.df <- miRNA.df[order(miRNA.df$Identity),]#Order
miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_bact),]$Identity == row.names(FR_ara_bact)#Check - all True
FR_ara_bact$Sequence <- miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_bact), "Sequence"] #Add sequence to table
FR_ara_bact$miRNA <- miRNA.df[miRNA.df$Identity %in% row.names(FR_ara_bact), "miRNA"] #Add sequence to table
FR_ara_bact <- FR_ara_bact %>% #Concatenate identical mature miRNAs
  group_by(Sequence) %>%
  mutate(miRNAs = paste0(miRNA, collapse = ";")) %>%
  group_by(miRNAs,Sequence)%>%
  summarise(across(1:6, sum))
colSums(FR_ara_bact[,3:8])==FR_ara_bact_colsums#Check

FR_ara_bact <- FR_ara_bact[(FR_ara_bact[,4]>5 & FR_ara_bact[,5]>5 & FR_ara_bact[,6]>5), ] #Keep only >5 reads
FR_ara_bact_norm <- data.frame(apply(FR_ara_bact[,3:8], 1, "/", FR_ara_bact_colsums))#normalize
rowSums(FR_ara_bact_norm)#sanity check
colnames(FR_ara_bact_norm) <- FR_ara_bact$miRNAs

#Prepare for plot
FR_ara_bact_norm <- FR_ara_bact_norm[2:4,]
FR_ara_bact_long <- gather(FR_ara_bact_norm)#long format for plotting
colnames(FR_ara_bact_long) <- c("miRNA", "RelativeAbundance")
FR_ara_bact_long$miRNA <- gsub("miR162a-3p;miR162b-3p", "miR162a,b", FR_ara_bact_long$miRNA)
FR_ara_bact_long <- FR_ara_bact_long[order(FR_ara_bact_long$miRNA),]

#Plot
ara_bact_miRNA_box <- ggplot(FR_ara_bact_long, aes(x=reorder(miRNA, -RelativeAbundance), y=100*RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Proportion of Arabidopsis reads (%)")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.y =  element_text(size=8))
ara_bact_miRNA_box
