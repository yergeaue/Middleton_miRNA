###Arrange miRNA abundance tables, keep only abundant miRNA with >10 reads in all reps and absent in unplanted (1 read or less)

##Arabidopsis-FR-Rhizosphere
FR_ara_prelim <- FR_ara[,colnames(FR_ara)%in%row.names(FR_ara_map[FR_ara_map$Other=="Preliminary.exp",])]
FR_ara_prelim_colsums <- colSums(FR_ara_prelim)#calculate colsums for normalization later
FR_ara_prelim <- FR_ara_prelim[rowSums(FR_ara_prelim)!=0,] #Remove zero only lines
FR_ara_prelim <- FR_ara_prelim[(FR_ara_prelim[,3]>10 & FR_ara_prelim[,4]>10 & FR_ara_prelim[,5]>10), ] #Keep only >10 reads
FR_ara_prelim_norm <- data.frame(apply(FR_ara_prelim, 1, "/", FR_ara_prelim_colsums))#normalize
rowSums(FR_ara_prelim_norm)#sanity check

#Prepare for plot
colnames(FR_ara_prelim_norm) <- gsub("_.*", "", colnames(FR_ara_prelim_norm))
FR_ara_prelim_norm <- FR_ara_prelim_norm[3:5,]
row.names(FR_ara_prelim_norm) <- c("Rep1", "Rep2", "Rep3")
FR_ara_prelim_norm$ath.miR165a.b <- FR_ara_prelim_norm$ath.miR165a.3p+FR_ara_prelim_norm$ath.miR165a #Merge: have same mature miRNA
FR_ara_prelim_norm$ath.miR166a.b.c.d.e.f.g <- FR_ara_prelim_norm$ath.miR166f+FR_ara_prelim_norm$ath.miR166g #Merge: have same mature miRNA
FR_ara_prelim_norm <- FR_ara_prelim_norm[,-c(1,11,15,16)]
FR_ara_prelim_long <- gather(FR_ara_prelim_norm)#long format for plotting
colnames(FR_ara_prelim_long) <- c("miRNA", "RelativeAbundance")
FR_ara_prelim_long <- FR_ara_prelim_long[order(FR_ara_prelim_long$miRNA),]

#Plot
ara_miRNA_box <- ggplot(FR_ara_prelim_long, aes(x=reorder(miRNA, -RelativeAbundance), y=RelativeAbundance))+
              geom_boxplot()+
              theme_bw()+
              ylab("Relative abundance")+
              xlab("")+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ara_miRNA_box

##Brachypodium-FR-Rhizosphere
FR_bra_prelim <- FR_bra[,colnames(FR_bra)%in%row.names(FR_bra_map[FR_bra_map$Experiment=="Preliminary.exp",])]
FR_bra_prelim_colsums <- colSums(FR_bra_prelim)#calculate colsums for normalization later
FR_bra_prelim <- FR_bra_prelim[rowSums(FR_bra_prelim)!=0,] #Remove zero only lines
FR_bra_prelim <- FR_bra_prelim[(FR_bra_prelim[,1]>10 & FR_bra_prelim[,1]>10), ] #Keep only >10 reads
FR_bra_prelim_norm <- data.frame(apply(FR_bra_prelim, 1, "/", FR_bra_prelim_colsums))#normalize
rowSums(FR_bra_prelim_norm)#sanity check

#Prepare for plot
colnames(FR_bra_prelim_norm) <- gsub("_.*", "", colnames(FR_bra_prelim_norm))
FR_bra_prelim_norm <- FR_bra_prelim_norm[1:2,]
row.names(FR_bra_prelim_norm) <- c("Rep1", "Rep3")
FR_bra_prelim_norm$ath.miR165a.b <- FR_ara_prelim_norm$ath.miR165a.3p+FR_ara_prelim_norm$ath.miR165a #Merge: have same mature miRNA
FR_bra_prelim_norm$bdi.miR166a.b.c.d.i <- FR_bra_prelim_norm[,1]+FR_bra_prelim_norm[,2]+FR_bra_prelim_norm[,4]+FR_bra_prelim_norm[,12]+FR_bra_prelim_norm[,13] #Merge: have same mature miRNA
FR_bra_prelim_norm$bdi.miR167c.d.e.g <- FR_bra_prelim_norm[,3]+FR_bra_prelim_norm[,5]+FR_bra_prelim_norm[,11]+FR_bra_prelim_norm[,16] #Merge: have same mature miRNA
FR_bra_prelim_norm$bdi.miR159b.c.d.e.f.g.h.1 <- FR_bra_prelim_norm[,8]
FR_bra_prelim_norm$bdi.miR396a.b.c.d <- FR_bra_prelim_norm[,15]+FR_bra_prelim_norm[,17]
FR_bra_prelim_norm <- FR_bra_prelim_norm[,-c(1,2,3,4,5,8,11,12,13,15,16,17)]
FR_bra_prelim_long <- gather(FR_bra_prelim_norm)#long format for plotting
colnames(FR_bra_prelim_long) <- c("miRNA", "RelativeAbundance")
FR_bra_prelim_long <- FR_bra_prelim_long[order(FR_bra_prelim_long$miRNA),]

#Plot
bra_miRNA_box <- ggplot(FR_bra_prelim_long, aes(x=reorder(miRNA, -RelativeAbundance), y=RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
bra_miRNA_box

#Look for overlap
intersect(gsub("\\..*","",gsub("ath.","", colnames(FR_ara_prelim_norm))), gsub("\\..*","",gsub("bdi.","", colnames(FR_bra_prelim_norm))))

##Arabidospis root Canada
CAN_ara_nitro <- CAN_ara[,colnames(CAN_ara)%in%row.names(CAN_map[(CAN_map$Compartment=="Endosphere" & CAN_map$Nitrogen=="No.Added.Nitrogen" & CAN_map$Plant=="Arabidopsis"),])]
CAN_ara_nitro_colsums <- colSums(CAN_ara_nitro)#calculate colsums for normalization later
CAN_ara_nitro <- CAN_ara_nitro[rowSums(CAN_ara_nitro)!=0,] #Remove zero only lines
CAN_ara_nitro <- CAN_ara_nitro[(CAN_ara_nitro[,1]>10 & CAN_ara_nitro[,2]>10 & CAN_ara_nitro[,3]>10& CAN_ara_nitro[,4]>10& CAN_ara_nitro[,5]>10), ] #Keep only >10 reads
CAN_ara_nitro_norm <- data.frame(apply(CAN_ara_nitro, 1, "/", CAN_ara_nitro_colsums))#normalize
rowSums(CAN_ara_nitro_norm)#sanity check

#Prepare for plot
colnames(CAN_ara_nitro_norm) <- gsub("_.*", "", colnames(CAN_ara_nitro_norm))
#CAN_ara_nitro_norm <- CAN_ara_nitro_norm[3:5,]
#row.names(CAN_ara_nitro_norm) <- c("Rep1", "Rep2", "Rep3")
CAN_ara_nitro_long <- gather(CAN_ara_nitro_norm)#long format for plotting
colnames(CAN_ara_nitro_long) <- c("miRNA", "RelativeAbundance")
CAN_ara_nitro_long <- CAN_ara_nitro_long[order(CAN_ara_nitro_long$miRNA),]

#Plot 
root_ara_miRNA_box <- ggplot(CAN_ara_nitro_long, aes(x=reorder(miRNA, -RelativeAbundance), y=RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
root_ara_miRNA_box

#Look for overlap
intersect(colnames(FR_ara_prelim_norm), colnames(CAN_ara_nitro_norm))

##Bacteria from soil and rhizosphere - France, use only 020320 - other dates have some weird samples
FR_ara_bact <- FR_ara[,colnames(FR_ara)%in%row.names(FR_ara_map[FR_ara_map$Date=="20320",])]
FR_ara_bact <- FR_ara_bact[rowSums(FR_ara_bact)!=0,] #Remove zero only lines
FR_ara_bact_colsums <- colSums(FR_ara_bact)#calculate colsums for normalization later
FR_ara_bact <- FR_ara_bact[(FR_ara_bact[,2]>5 & FR_ara_bact[,3]>5 & FR_ara_bact[,4]>5), ] #Keep only >5 reads
FR_ara_bact_norm <- data.frame(apply(FR_ara_bact, 1, "/", FR_ara_bact_colsums))#normalize
rowSums(FR_ara_bact_norm)#sanity check


#Prepare for plot
colnames(FR_ara_bact_norm) <- gsub("_.*", "", colnames(FR_ara_bact_norm))
FR_ara_bact_norm <- FR_ara_bact_norm[2:4,]
FR_ara_bact_long <- gather(FR_ara_bact_norm)#long format for plotting
colnames(FR_ara_bact_long) <- c("miRNA", "RelativeAbundance")
FR_ara_bact_long <- FR_ara_bact_long[order(FR_ara_bact_long$miRNA),]

#Plot
ara_bact_miRNA_box <- ggplot(FR_ara_bact_long, aes(x=reorder(miRNA, -RelativeAbundance), y=RelativeAbundance))+
  geom_boxplot()+
  theme_bw()+
  ylab("Relative abundance")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ara_bact_miRNA_box
