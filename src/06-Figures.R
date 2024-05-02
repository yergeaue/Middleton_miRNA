###Figures

#Figure 1 - miRNA are present in the rhizosphere
fig1 <- ggarrange(ara_miRNA_box,bra_miRNA_box,root_ara_miRNA_box, labels = c("A","B", "C"), common.legend = F, nrow=3, ncol=1)
fig1
ggsave(fig1, filename = here("output", "figs", "fig1.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 7, units = "in")


#Figure 2 - miRNA are taken up by rhizosphere bacteria
vario <- readTIFF(source =  here("data", "raw", "Variovorax_zoom.tiff")) #Include microscopy images
vario.g <- rasterGrob(vario, interpolate = TRUE)
bacillus <- readTIFF(source =  here("data", "raw", "Bacillus_zoom.tiff")) #Include microscopy images
bacillus.g <- rasterGrob(bacillus, interpolate = TRUE)
fig2 <- ggarrange(ara_bact_miRNA_box, vario.g, bacillus.g, labels = c("A", "B", "C"), ncol = 1, font.label = list(size = 18), nrow=3 )
fig2
ggsave(fig2, filename = here("output", "figs", "fig2.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 7, units = "in")

#Figure 3 - miRNA change the bacterial transcriptome
fig3 <- ggarrange(ggarrange(vario.volcano.20,vario.volcano.120, ncol=2, nrow=1, 
                            labels = c("A","B"), common.legend = TRUE, legend = "bottom"), 
                  qpcr.plot, ncol=1, nrow = 2, labels = c("","C"))
fig3
ggsave(fig3, filename = here("output", "figs", "fig3.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")

#Figure 4 - miRNA affect the bacterial community
fig4 <- ggarrange(mut_StackedBarPlot_phylum_rel, mut.div.plot, 
                  miPEP_StackedBarPlot_phylum_rel, boxMix17_mean, stack_mean_mixAA, 
                  box_sig2, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))
fig4
ggsave(fig4, filename = here("output", "figs", "fig4.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")
