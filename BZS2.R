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

