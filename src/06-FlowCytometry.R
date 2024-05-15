#Calculating mean abundances 
FACS_df_mean<-FACS_df %>%
  group_by(Treatment,Day, Strain) %>%
  summarise_if(is.numeric, mean) #27 obs. 15 variables

#Boxplots: Blue_Green_Red
boxall<- ggplot(FACS_df_mean, aes(x=Treatment, y=RBG_P3)) +
  geom_boxplot(outlier.shape = NA, alpha=0.6, aes(fill=Strain, color=Strain )) +
  geom_point(aes(fill=Strain, color=Strain, shape=Strain), position = position_jitterdodge())+theme_minimal()+
  theme(strip.text = element_text(colour = 'grey30', size=10,face = "italic" ))+
  scale_fill_manual(values = c("grey50","grey10"))+
  scale_colour_manual(values = c("grey50","grey10"))
boxall
boxall_2<-boxall +labs(x="", y="Percentage of active Cy5-positive \n bacteria (%)")+ theme(legend.title = element_text(size = 12) ,legend.text = element_text(size=12))+scale_x_discrete(labels=c('pCp-Cy5', 'Plant miRNA', 'Scrambled miRNA'))
boxall_2

boxP3_R_mean<- ggplot(FACS_df_mean, aes(x=Treatment, y=P3_R_Median)) +
  geom_boxplot(outlier.shape = NA, alpha=0.6, aes(fill=Strain, color=Strain )) +
  geom_point(aes(fill=Strain, color=Strain,shape=Strain), position = position_jitterdodge())+theme_minimal()+
  theme(strip.text = element_text(colour = 'grey30', size=10,face = "italic" ))+
  scale_fill_manual(values = c("#f28482","#d90429"))+
  scale_colour_manual(values = c("#f28482","#d90429"))
boxP3_R_mean
boxP3_R_mean_2<-boxP3_R_mean +labs(x="", y="Median of Cy5 intensity")+ theme(legend.title = element_text(size = 12) ,legend.text = element_text(size=12))+scale_x_discrete(labels=c('pCp-Cy5', 'Plant miRNA', 'Scrambled miRNA'))
boxP3_R_mean_2

#Computing the statistical test
stat.test.Pop <- FACS_df_mean %>%
  group_by(Strain) %>%
  dunn_test(RBG_P3 ~ Treatment, p.adjust.method = "BH") %>%
  add_significance()
stat.test.Pop


#Adding the coordinates of the p-values
stat.test.pos<-stat.test.Pop %>% add_xy_position(x="Treatment", group="Strain",fun = "max")
box_sig.Pop<- boxall_2+stat_pvalue_manual(
  stat.test.pos, hide.ns = TRUE, 
  label = "p.adj.signif",tip.length = 0, size=6)
box_sig.Pop


#Computing the statistical test
stat.test.MFI <- FACS_df_mean %>%
  group_by(Strain) %>%
  dunn_test(P3_R_Median ~ Treatment, p.adjust.method = "BH") %>%
  add_significance()
stat.test.MFI

#Adding the coordinates of the p-values
stat.test.pos_2<-stat.test.MFI%>% add_xy_position(x="Treatment", group="Strain",fun = "max")
box_sig.MFI<- boxP3_R_mean_2+stat_pvalue_manual(
  stat.test.pos_2, hide.ns = TRUE, 
  label='p.adj.signif',tip.length = 0, size=6)
box_sig.MFI
