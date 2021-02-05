
realone <- ggplot(melted_cormat, aes("", "", fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="spearman\nCorrelation") +
  coord_fixed()

realone



ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="SparCC\nCorrelation") +
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
        
# Print the heatmap
print(ggheatmap)







ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="SparCC\nCorrelation") +
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        panel.background = element_blank())

# Print the heatmap
print(ggheatmap)





labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", oral_sparcc_Modules),
               yLabels = paste(" ", oral_clr1_Modules),
               colorLabels = TRUE,
               xSymbols = paste("SparCC ", oral_sparcc_Modules, ": ", oral_sparcc_ModTotals, sep=""),
               ySymbols = paste("CLR1 ", oral_clr1_Modules, ": ", oral_clr1_ModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)




plot(oral_trimmed_sparcc, 
     vertex.label ="",
     vertex.size= V(oral_net_sparcc)$strength * 10,
     vertex.label.color="black",
     vertex.color = adjustcolor(V(oral_net_sparcc)$color, alpha.f=0.6)) # make nodes slightly opaque to improve readability






plot(oral_trimmed_clr1, 
     vertex.size= V(oral_net_clr1)$strength * 3,
     vertex.label="",
     vertex.label.color="black",
     vertex.color = adjustcolor(V(oral_net_clr1)$color, alpha.f=0.6)) # make nodes slightly opaque to improve readability
     



plot(oral_trimmed_comp, 
     vertex.size= V(oral_net_comp)$strength * 3,
     vertex.label="",
     vertex.label.color="black",
     vertex.color = adjustcolor(V(oral_net_comp)$color, alpha.f=0.6)) # make nodes slightly opaque to improve readability

# Individual modules

plot(brown_net_sparcc, 
     vertex.label.color="black",
     vertex.color = adjustcolor(module, alpha.f=0.6)) # make nodes slightly opaque to improve readability)





plot(blue_net_clr1,
     vertex.label.color="black",
     vertex.color = adjustcolor("blue", alpha.f=0.5)) # make nodes slightly opaque to improve readability)



plot(net_clr1, 
     vertex.label.color="black",
     vertex.color = adjustcolor("brown", alpha.f=0.6)) # make nodes slightly opaque to improve readability)




plot(net_comp, 
     vertex.label.color="black",
     vertex.color = adjustcolor("pink", alpha.f=0.6)) # make nodes slightly opaque to improve readability)





module="red"

probes = colnames(oral_dat)


inModule = is.finite(match(oral_moduleColors_sparcc, module));
modProbes= probes[inModule]

now <- toString(modProbes)
length(modProbes)
now


oral_comp_Modules



module="turquoise"

probes = colnames(oral_dat)


inModule = is.finite(match(oral_moduleColors_sparcc, module));
modProbes= probes[inModule]

now <- toString(modProbes)
length(modProbes)
now


oral_sparcc_Modules




