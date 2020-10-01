#import libraries
library(factoextra)
library(ggcorrplot)

#read csv
bzs2 <- read.csv('BZS2_transformation.csv')
bzs2 <- as.data.frame(bzs2)

#compute PCA
df_PCA.variables <- bzs2[c('BZTalanine','glycosylatedBZT',
                               'BZTacetylalanine','aniline','amino_3_phenol',
                               'phenazine','methylBZT','methoxyBZT')]

res.pca <- prcomp(df_PCA.variables, scale = TRUE)

#biplot, grouped by genotype
biplot.genotype <- fviz_pca_biplot(res.pca, label="var", habillage=bzs2$Genotype, col.var = "#000000",
                                          addEllipses=TRUE, legend.title = "Genotype", ellipse.type = "confidence", 
                                          pointsize = 5, labelsize = 7, arrowsize=1.5) + 
  theme_classic() + theme(text = element_text(size = 20))

#correlation plot, principal components and location descriptors
cor_variables <- cbind(res.pca$x, km = bzs2$km, density = bzs2$density)
cor_plot <- ggcorrplot(cor(cor_variables), p.mat = cor_pmat(cor_variables), hc.order=FALSE, type='lower')