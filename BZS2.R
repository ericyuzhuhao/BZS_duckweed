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

#linear models for each response variable
#responses are scaled to avoid error in MCMCglmm ('Mixed model equations singular: use a (stronger) prior') 
library(MCMCglmm)

#linear model for percent decrease in benzotriazole concentration
bzs_model <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density +
                        Salt:Microbes + Salt:density + Microbes:density +
                        Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp1 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density +
                        Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp2 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density +
                              Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp3 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density +
                              Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp4 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp5 <- MCMCglmm(scale(BZT_percent_d)~Salt + density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model.temp6 <- MCMCglmm(scale(BZT_percent_d)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole alanine
BZTalanine_model <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density +
                        Salt:Microbes + Salt:density + Microbes:density +
                        Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp1 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density +
                               Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp2 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density +
                                     Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp3 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density +
                                     Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp4 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp5 <- MCMCglmm(scale(BZTalanine)~Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model.temp6 <- MCMCglmm(scale(BZTalanine)~density, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of glycosylated benzotriazole 
glycosylatedBZT_model <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density +
                               Salt:Microbes + Salt:density + Microbes:density +
                               Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp1 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density +
                                    Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp2 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density +
                                          Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp3 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density +
                                          Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp4 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp5 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model.temp6 <- MCMCglmm(scale(glycosylatedBZT)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole acetyl-alanine
BZTacetylalanine_model <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density +
                                    Salt:Microbes + Salt:density + Microbes:density +
                                    Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of aniline
aniline_model <- MCMCglmm(scale(aniline)~Salt + Microbes + density +
                                     Salt:Microbes + Salt:density + Microbes:density +
                                     Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp1 <- MCMCglmm(scale(aniline)~Salt + Microbes + density +
                            Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp2 <- MCMCglmm(scale(aniline)~Salt + Microbes + density +
                                  Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp3 <- MCMCglmm(scale(aniline)~Salt + Microbes + density +
                                  Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp4 <- MCMCglmm(scale(aniline)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp5 <- MCMCglmm(scale(aniline)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model.temp6 <- MCMCglmm(scale(aniline)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of amino_3_phenol 
amino_3_phenol_model <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density +
                            Salt:Microbes + Salt:density + Microbes:density +
                            Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp1 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density +
                                   Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp2 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density +
                                         Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp3 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density +
                                         Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp4 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp5 <- MCMCglmm(scale(amino_3_phenol)~Salt + density, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model.temp6 <- MCMCglmm(scale(amino_3_phenol)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of phenazine
phenazine_model <- MCMCglmm(scale(phenazine)~Salt + Microbes + density +
                                   Salt:Microbes + Salt:density + Microbes:density +
                                   Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp1 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density +
                              Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp2 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density +
                                    Salt:Microbes + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp3 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density +
                                    Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp4 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp5 <- MCMCglmm(scale(phenazine)~Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model.temp6 <- MCMCglmm(scale(phenazine)~density, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of methylbenzotriazole
methylBZT_model <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density +
                              Salt:Microbes + Salt:density + Microbes:density +
                              Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp1 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density +
                              Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp2 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density +
                                    Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp3 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density +
                                    Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp4 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp5 <- MCMCglmm(scale(methylBZT)~Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model.temp6 <- MCMCglmm(scale(methylBZT)~Microbes, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of methoxybenzotriazole
methoxyBZT_model <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density +
                              Salt:Microbes + Salt:density + Microbes:density +
                              Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model.temp1 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density +
                               Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model.temp2 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density +
                                     Salt:Microbes + Microbes:density,random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model.temp3 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density +
                                     Microbes:density,random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model.temp4 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model.temp5 <- MCMCglmm(scale(methoxyBZT)~Salt + density, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for frond number on day 10
X10FN_model <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density +
                               Salt:Microbes + Salt:density + Microbes:density +
                               Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp1 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density +
                          Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp2 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density +
                                Salt:Microbes + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp3 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density +
                                Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp4 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp5 <- MCMCglmm(scale(X10.FN)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model.temp6 <- MCMCglmm(scale(X10.FN)~Microbes, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10
od600_model <- MCMCglmm(scale(log(od600))~Salt + Microbes + density +
                          Salt:Microbes + Salt:density + Microbes:density +
                          Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp1 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density +
                          Salt:Microbes + Salt:density + Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp2 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density +
                                Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp3 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density +
                                Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp4 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp5 <- MCMCglmm(scale(log(od600))~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model.temp6 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10, restricting to inoculated wells
od600i_model <- MCMCglmm(scale(log(od600))~Salt + density + Salt:density, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model.temp1 <- MCMCglmm(scale(log(od600))~Salt + density, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model.temp2 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

#linear model for percent decrease in benzotriazole concentration, with distance to city center as the location descriptor
bzs_model_km <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + km +
                        Salt:Microbes + Salt:km + Microbes:km +
                        Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_km.temp1 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + km +
                           Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_km.temp2 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + km +
                                 Salt:Microbes + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

#salt microbes interactions significant when km is included
bzs_model_km.temp3 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + km +
                                 Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_km.temp4 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes +
                                 Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_km.temp5 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

#salt is significant
bzs_model_km.temp6 <- MCMCglmm(scale(BZT_percent_d)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole alanine, with distance to city center as the location descriptor
BZTalanine_model_km <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + km +
                               Salt:Microbes + Salt:km + Microbes:km +
                               Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp1 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + km +
                                  Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp2 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + km +
                                        Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp3 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + km + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp4 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp5 <- MCMCglmm(scale(BZTalanine)~Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_km.temp6 <- MCMCglmm(scale(BZTalanine)~km, random = ~ Genotype , data=bzs2,verbose=F)


#linear model for amount of glycosylated benzotriazole, with distance to city center as the location descriptor 
glycosylatedBZT_model_km <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km +
                                    Salt:Microbes + Salt:km + Microbes:km +
                                    Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp1 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km +
                                       Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp2 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km +
                                             Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp3 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km  + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp4 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp5 <- MCMCglmm(scale(glycosylatedBZT)~Salt + km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_km.temp6 <- MCMCglmm(scale(glycosylatedBZT)~km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole acetyl-alanine, with distance to city center as the location descriptor 
BZTacetylalanine_model_km <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + km +
                                     Salt:Microbes + Salt:km + Microbes:km +
                                     Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of aniline, with distance to city center as the location descriptor 
aniline_model_km <- MCMCglmm(scale(aniline)~Salt + Microbes + km +
                            Salt:Microbes + Salt:km + Microbes:km +
                            Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp1 <- MCMCglmm(scale(aniline)~Salt + Microbes + km +
                               Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp2 <- MCMCglmm(scale(aniline)~Salt + Microbes + km + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp3 <- MCMCglmm(scale(aniline)~Salt + Microbes + km + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp4 <- MCMCglmm(scale(aniline)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp5 <- MCMCglmm(scale(aniline)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_km.temp6 <- MCMCglmm(scale(aniline)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of amino_3_phenol, with distance to city center as the location descriptor  
amino_3_phenol_model_km <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + km +
                                   Salt:Microbes + Salt:km + Microbes:km +
                                   Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_km.temp1 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + km +
                                      Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_km.temp2 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + km + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_km.temp3 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_km.temp4 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_km.temp5 <- MCMCglmm(scale(amino_3_phenol)~Salt + km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of phenazine, with distance to city center as the location descriptor 
phenazine_model_km <- MCMCglmm(scale(phenazine)~Salt + Microbes + km +
                              Salt:Microbes + Salt:km + Microbes:km +
                              Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp1 <- MCMCglmm(scale(phenazine)~Salt + Microbes + km +
                                 Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp2 <- MCMCglmm(scale(phenazine)~Salt + Microbes + km +
                                       Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp3 <- MCMCglmm(scale(phenazine)~Salt + Microbes + km + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp4 <- MCMCglmm(scale(phenazine)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp5 <- MCMCglmm(scale(phenazine)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_km.temp6 <- MCMCglmm(scale(phenazine)~Microbes, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of methylbenzotriazole, with distance to city center as the location descriptor 
methylBZT_model_km <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km +
                              Salt:Microbes + Salt:km + Microbes:km +
                              Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp1 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km +
                                 Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp2 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km +
                                       Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp3 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp4 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp5 <- MCMCglmm(scale(methylBZT)~Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_km.temp6 <- MCMCglmm(scale(methylBZT)~km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of methoxybenzotriazole, with distance to city center as the location descriptor 
methoxyBZT_model_km <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km +
                               Salt:Microbes + Salt:km + Microbes:km +
                               Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_km.temp1 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km +
                                  Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_km.temp2 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km +
                                        Salt:Microbes + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_km.temp3 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_km.temp4 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_km.temp5 <- MCMCglmm(scale(methoxyBZT)~Salt + km, random = ~ Genotype , data=bzs2,verbose=F)


#linear model for frond number on day 10, with distance to city center as the location descriptor 
X10FN_model_km <- MCMCglmm(scale(X10.FN)~Salt + Microbes + km +
                          Salt:Microbes + Salt:km + Microbes:km +
                          Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp1 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + km +
                             Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp2 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + km +
                                   Salt:Microbes + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp3 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp4 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp5 <- MCMCglmm(scale(X10.FN)~Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_km.temp6 <- MCMCglmm(scale(X10.FN)~Microbes, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10, with distance to city center as the location descriptor 
od600_model_km <- MCMCglmm(scale(log(od600))~Salt + Microbes + km +
                          Salt:Microbes + Salt:km + Microbes:km +
                          Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp1 <- MCMCglmm(scale(log(od600))~Salt + Microbes + km +
                             Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp2 <- MCMCglmm(scale(log(od600))~Salt + Microbes + km +
                                   Salt:Microbes + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp3 <- MCMCglmm(scale(log(od600))~Salt + Microbes + km +
                                   Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp4 <- MCMCglmm(scale(log(od600))~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp5 <- MCMCglmm(scale(log(od600))~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_km.temp6 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10, restricting to inoculated wells, with distance to city center as the location descriptor 
od600i_model_km <- MCMCglmm(scale(log(od600))~Salt + km + Salt:km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_km.temp1 <- MCMCglmm(scale(log(od600))~Salt + km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_km.temp2 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)



#linear model for percent decrease in benzotriazole concentration, both location descriptors
bzs_model_both <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                             Salt:Microbes + Salt:density + Salt:km +
                             Microbes:density + Microbes:km +
                             density:km +
                             Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                             Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp1 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                             Salt:Microbes + Salt:density + Salt:km +
                             Microbes:density + Microbes:km +
                             density:km +
                             Salt:Microbes:density + Salt:Microbes:km + 
                             Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp2 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density + Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp3 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp4 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp5 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + 
                                   Microbes:density + Microbes:km +
                                   density:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp6 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + 
                                   Microbes:density + Microbes:km +
                                   density:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp7 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + 
                                   Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp8 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp9 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + km +
                                   Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp10 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density + 
                                   Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp11 <- MCMCglmm(scale(BZT_percent_d)~Salt + Microbes + density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp12 <- MCMCglmm(scale(BZT_percent_d)~Salt + density, random = ~ Genotype , data=bzs2,verbose=F)

bzs_model_both.temp13 <- MCMCglmm(scale(BZT_percent_d)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole alanine, both location descriptors
BZTalanine_model_both <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                    Salt:Microbes + Salt:density + Salt:km +
                                    Microbes:density + Microbes:km +
                                    density:km +
                                    Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                    Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_both.temp1 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                    Salt:Microbes + Salt:density + Salt:km +
                                    Microbes:density + Microbes:km +
                                    density:km +
                                    Salt:Microbes:density + Salt:Microbes:km +
                                    Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_both.temp2 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Salt:Microbes:density + 
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_both.temp3 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_both.temp4 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                          Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTalanine_model_both.temp5 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                          Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)
BZTalanine_model_both.temp6 <- MCMCglmm(scale(BZTalanine)~Salt + Microbes + density + km +
                                          Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of glycosylated benzotriazole, both location descriptors
glycosylatedBZT_model_both <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                       Salt:Microbes + Salt:density + Salt:km +
                                       Microbes:density + Microbes:km +
                                       density:km +
                                       Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                       Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp1 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         density:km +
                                         Salt:Microbes:density + Salt:density:km +
                                         Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp2 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:density + Salt:km +
                                               Microbes:density + Microbes:km +
                                               density:km +
                                               Salt:Microbes:density + Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp3 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:density + Salt:km +
                                               Microbes:density + Microbes:km +
                                               density:km +
                                               Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp4 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:density + Salt:km +
                                               Microbes:density + Microbes:km +
                                               density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp5 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:km +
                                               Microbes:density + Microbes:km +
                                               density:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp6 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:km +
                                               Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp7 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:km + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp8 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp9 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km +
                                               Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp10 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + density + km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp11 <- MCMCglmm(scale(glycosylatedBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp12 <- MCMCglmm(scale(glycosylatedBZT)~Salt + km, random = ~ Genotype , data=bzs2,verbose=F)

glycosylatedBZT_model_both.temp13 <- MCMCglmm(scale(glycosylatedBZT)~km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of benzotriazole acetyl-alanine, both location descriptors
BZTacetylalanine_model_both <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp1 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Salt:Microbes:density + Salt:density:km +
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp2 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                                Salt:Microbes + Salt:density + Salt:km +
                                                Microbes:density + Microbes:km +
                                                density:km +
                                                Salt:Microbes:density + 
                                                Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp3 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                                Salt:Microbes + Salt:density + Salt:km +
                                                Microbes:density + Microbes:km +
                                                density:km +
                                                Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp4 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                                Salt:Microbes + Salt:density + Salt:km +
                                                Microbes:density + 
                                                density:km +
                                                Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp5 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                                Salt:Microbes + Salt:density + Salt:km +
                                                Microbes:density + 
                                                Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

BZTacetylalanine_model_both.temp6 <- MCMCglmm(scale(BZTacetylalanine)~Salt + Microbes + density + km +
                                                Salt:Microbes + Salt:density + 
                                                Microbes:density + 
                                                Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of aniline, both location descriptors
aniline_model_both <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                 Salt:Microbes + Salt:density + Salt:km +
                                 Microbes:density + Microbes:km +
                                 density:km +
                                 Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                 Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp1 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                 Salt:Microbes + Salt:density + Salt:km +
                                 Microbes:density + Microbes:km +
                                 density:km +
                                 Salt:Microbes:km + Salt:density:km +
                                 Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp2 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:Microbes + Salt:density + Salt:km +
                                       Microbes:density + Microbes:km +
                                       density:km +
                                       Salt:Microbes:km + Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp3 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:Microbes + Salt:density + Salt:km +
                                       Microbes:density + Microbes:km +
                                       density:km +
                                       Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp4 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:Microbes + Salt:density + Salt:km +
                                       Microbes:density + Microbes:km +
                                       density:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp5 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:Microbes + Salt:density + Salt:km +
                                       Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp6 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:density + Salt:km +
                                       Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp7 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Salt:km +
                                       Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp8 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp9 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km +
                                       Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp10 <- MCMCglmm(scale(aniline)~Salt + Microbes + density + km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp11 <- MCMCglmm(scale(aniline)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp12 <- MCMCglmm(scale(aniline)~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

aniline_model_both.temp13 <- MCMCglmm(scale(aniline)~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of amino_3_phenol, both location descriptors
amino_3_phenol_model_both <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density + km +
                                        Salt:Microbes + Salt:density + Salt:km +
                                        Microbes:density + Microbes:km +
                                        density:km +
                                        Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                        Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_both.temp1 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density + km +
                                        Salt:Microbes + Salt:density + Salt:km +
                                        Microbes:density + Microbes:km +
                                        density:km +
                                        Salt:Microbes:density + Salt:Microbes:km + 
                                        Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_both.temp2 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density + km +
                                              Salt:Microbes + Salt:density + Salt:km +
                                              Microbes:density + Microbes:km +
                                              density:km +
                                              Salt:Microbes:density + Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

amino_3_phenol_model_both.temp3 <- MCMCglmm(scale(amino_3_phenol)~Salt + Microbes + density + km +
                                              Salt:Microbes + Salt:density + Salt:km +
                                              Microbes:density + Microbes:km +
                                              Salt:Microbes:density + Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for amount of phenazine, both location descriptors
phenazine_model_both <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                   Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp1 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density + Salt:density:km +
                                   Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp2 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         density:km +
                                         Salt:density:km +
                                         Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp3 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         density:km +
                                         Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp4 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                         Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         density:km +
                                         Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp5 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                         Salt:density + Salt:km +
                                         Microbes:density +
                                         density:km +
                                         Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp6 <- MCMCglmm(scale(phenazine)~Salt + Microbes + density + km +
                                         Salt:density + Salt:km +
                                         density:km +
                                         Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp7 <- MCMCglmm(scale(phenazine)~Salt + density + km +
                                         Salt:density + Salt:km +
                                         density:km +
                                         Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp8 <- MCMCglmm(scale(phenazine)~Salt + density + km +
                                         Salt:density + Salt:km +
                                         density:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp9 <- MCMCglmm(scale(phenazine)~Salt + density + km +
                                         Salt:density + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp10 <- MCMCglmm(scale(phenazine)~Salt + density + km + Salt:km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp11 <- MCMCglmm(scale(phenazine)~Salt + density + km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp12 <- MCMCglmm(scale(phenazine)~density + km, random = ~ Genotype , data=bzs2,verbose=F)

phenazine_model_both.temp13 <- MCMCglmm(scale(phenazine)~density, random = ~ Genotype , data=bzs2,verbose=F)


#linear model for amount of methylbenzotriazole, both location descriptors
methylBZT_model_both <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                   Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp1 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                   Salt:Microbes + Salt:density + Salt:km +
                                   Microbes:density + Microbes:km +
                                   density:km +
                                   Salt:Microbes:density + Salt:Microbes:km + 
                                   Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp2 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         density:km +
                                         Salt:Microbes:density + Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp3 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         Salt:Microbes:density + Salt:Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp4 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km +
                                         Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp5 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:Microbes + Salt:density + Salt:km +
                                         Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp6 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:density + Salt:km +
                                         Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp7 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:density + Salt:km +
                                         Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp8 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Salt:km +
                                         Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp9 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km +
                                         Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp10 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + density + km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp11 <- MCMCglmm(scale(methylBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp12 <- MCMCglmm(scale(methylBZT)~Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methylBZT_model_both.temp13 <- MCMCglmm(scale(methylBZT)~Microbes, random = ~ Genotype , data=bzs2,verbose=F)


#linear model for amount of methoxybenzotriazole, both location descriptors
methoxyBZT_model_both <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                    Salt:Microbes + Salt:density + Salt:km +
                                    Microbes:density + Microbes:km +
                                    density:km +
                                    Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                                    Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp1 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                    Salt:Microbes + Salt:density + Salt:km +
                                    Microbes:density + Microbes:km +
                                    density:km +
                                    Salt:Microbes:density + Salt:density:km +
                                    Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp2 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Salt:Microbes:density + 
                                          Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp3 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km +
                                          Salt:Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)


methoxyBZT_model_both.temp4 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:density + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp5 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + Salt:km +
                                          Microbes:density + Microbes:km +
                                          density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp6 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + 
                                          Microbes:density + Microbes:km +
                                          density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp7 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + 
                                          Microbes:density + 
                                          density:km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp8 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Salt:Microbes + 
                                          Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp9 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km +
                                          Microbes:density, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp10 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + density + km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp11 <- MCMCglmm(scale(methoxyBZT)~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

methoxyBZT_model_both.temp12 <- MCMCglmm(scale(methoxyBZT)~Salt + km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for frond number on day 10, both location descriptors
X10FN_model_both <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                               Salt:Microbes + Salt:density + Salt:km +
                               Microbes:density + Microbes:km +
                               density:km +
                               Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                               Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp1 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                               Salt:Microbes + Salt:density + Salt:km +
                               Microbes:density + Microbes:km +
                               density:km +
                               Salt:Microbes:density + Salt:density:km +
                               Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp2 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Salt:Microbes:density + 
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp3 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp4 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + 
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp5 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                                     Salt:Microbes +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp6 <- MCMCglmm(scale(X10.FN)~Salt + Microbes + density + km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

X10FN_model_both.temp7 <- MCMCglmm(scale(X10.FN)~Microbes + density + km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10, both location descriptors 
od600_model_both <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                               Salt:Microbes + Salt:density + Salt:km +
                               Microbes:density + Microbes:km +
                               density:km +
                               Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
                               Microbes:density:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp1 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                               Salt:Microbes + Salt:density + Salt:km +
                               Microbes:density + Microbes:km +
                               density:km +
                               Salt:Microbes:density + Salt:Microbes:km + Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp2 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Salt:Microbes:km + Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp3 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km +
                                     density:km +
                                     Salt:density:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp4 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km +
                                     density:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp5 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:density + Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp6 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + Salt:km +
                                     Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp7 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + Salt:density + 
                                     Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp8 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes + 
                                     Microbes:km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp9 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km +
                                     Salt:Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp10 <- MCMCglmm(scale(log(od600))~Salt + Microbes + density + km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp11 <- MCMCglmm(scale(log(od600))~Salt + Microbes + km, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp12 <- MCMCglmm(scale(log(od600))~Salt + Microbes, random = ~ Genotype , data=bzs2,verbose=F)

od600_model_both.temp13 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2,verbose=F)

#linear model for od600 on day 10, restricting to inoculated wells, both location descriptors
od600i_model_both <- MCMCglmm(scale(log(od600))~Salt + density + km +
                                Salt:density + Salt:km +
                                density:km +
                                Salt:density:km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp1 <- MCMCglmm(scale(log(od600))~Salt + density + km +
                                Salt:density + Salt:km +
                                density:km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp2 <- MCMCglmm(scale(log(od600))~Salt + density + km +
                                      Salt:density + 
                                      density:km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp3 <- MCMCglmm(scale(log(od600))~Salt + density + km +
                                      Salt:density, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp4 <- MCMCglmm(scale(log(od600))~Salt + density + km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp5 <- MCMCglmm(scale(log(od600))~Salt + km, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)

od600i_model_both.temp6 <- MCMCglmm(scale(log(od600))~Salt, random = ~ Genotype , data=bzs2[bzs2$Microbes == 'Yes',],verbose=F)