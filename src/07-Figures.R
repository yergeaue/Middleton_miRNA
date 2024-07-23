###Figures

#Figure 1 - miRNA are present in the rhizosphere and associate with bacteria
fig1 <- ggarrange(ara_miRNA_box,bra_miRNA_box,root_ara_miRNA_box,ara_bact_miRNA_box, labels = c("A","B", "C", "D"), common.legend = F, nrow=2, ncol=2)
fig1
ggsave(fig1, filename = here("output", "figs", "fig1.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 7, width = 7, units = "in")
ggsave(fig1, filename = here("output", "figs", "fig1.pdf"), dpi = 600, device = "pdf", height = 7, width = 7, units = "in")


#Figure 2 - transcriptomic and microscropy
#fig2 <- ggarrange(ggarrange(vario.volcano.20,vario.volcano.120, labels = c("A", "B"), common.legend = TRUE, legend = "bottom", nrow=1, ncol=2,font.label = list(size = 18) ),vario.g, bacillus.g, labels = c("", "C","D"), ncol = 1, font.label = list(size = 18), nrow=3 )
fig2 <- ggarrange(vario.volcano.20,vario.volcano.120, vario.g, bacillus.g, nrow=2, ncol=2,font.label = list(size = 18) ,labels = c("A", "B", "C","D") )
fig2
ggsave(fig2, filename = here("output", "figs", "fig2.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")
ggsave(fig2, filename = here("output", "figs", "fig2.pdf"), dpi = 600, device = "pdf", height = 14, width = 14, units = "in")


#Figure 3 - qPCR bacterial transcriptome
fig3 <- qpcr.plot
fig3
ggsave(fig3, filename = here("output", "figs", "fig3.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 3.5, width = 3.5, units = "in")
ggsave(fig3, filename = here("output", "figs", "fig3.pdf"), dpi = 600, device = "pdf", height = 3.5, width = 3.5, units = "in")

#Figure 4 - miRNA affect the bacterial community
fig4 <- ggarrange(mut_StackedBarPlot_phylum_rel, mut.div.plot, 
                  miPEP_StackedBarPlot_phylum_rel, boxMix17_mean, stack_mean_mixAA, 
                  box_sig2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
fig4
ggsave(fig4, filename = here("output", "figs", "fig4.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")
ggsave(fig4, filename = here("output", "figs", "fig4.pdf"), dpi = 600, device = "pdf", height = 14, width = 14, units = "in")


#Fig S2 - Flow cytometry results
figS2<-ggarrange(box_sig.Pop, box_sig.MFI,nrow=2, labels=(c('A','B')))
ggsave(figS2, file=here("output","figs", "figS2.tiff"),units="cm", width=14, height=18, compression='lzw' )
