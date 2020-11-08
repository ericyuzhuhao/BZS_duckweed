#import libraries
library(readxl)
library(MCMCglmm)

#read data
bzs1 <- read_xlsx('BZS1_transformation.xlsx')
bzs1 <- as.data.frame(bzs1)
area <- pi*25
bzs1['density'] <- bzs1['road_length']/area

#linear models for each response variable
#responses are scaled to avoid error in MCMCglmm ('Mixed model equations singular: use a (stronger) prior') 

#linear model for percent decrease in benzotriazole concentration
bzs_model <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                        BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                        Salt:Microbes + Salt:density +
                        Microbes:density + 
                        BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                        Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp1 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                        BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                        Salt:Microbes + Salt:density +
                        Microbes:density + 
                        BZT_init:Salt:Microbes + BZT_init:Salt:density + 
                        Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp2 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density + 
                              BZT_init:Salt:Microbes + BZT_init:Salt:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp3 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density + 
                              BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp4 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp5 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + 
                              Salt:Microbes + Salt:density +
                              Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp6 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + 
                              Salt:Microbes + Salt:density +
                              Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp7 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              Salt:Microbes + Salt:density +
                              Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp8 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              Salt:density +
                              Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model.temp9 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
                              Salt:density, random = ~ Genotype , data=bzs1,verbose=F)


#linear model for amount of benzotriazole alanine
BZTalanine_model <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density + 
                              BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                              Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp1 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Microbes:density + 
                               Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp2 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp3 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp4 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + 
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp5 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + 
                                     Salt:Microbes + 
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp6 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + 
                                     Salt:Microbes + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp7 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes + 
                                     Salt:Microbes + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model.temp8 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes + 
                                     Salt:Microbes + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of glycosylated benzotriazole 
glycosylatedBZT_model <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                    Salt:Microbes + Salt:density +
                                    Microbes:density + 
                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                                    Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp1 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                    Salt:Microbes + Salt:density +
                                    Microbes:density + 
                                    BZT_init:Salt:density + BZT_init:Microbes:density + 
                                    Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp2 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                          Salt:Microbes + Salt:density +
                                          Microbes:density + 
                                          BZT_init:Salt:density + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp3 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                          Salt:Microbes + Salt:density +
                                          Microbes:density + 
                                          BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp4 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                          Salt:Microbes + Salt:density +
                                          Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp5 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                          Salt:density +
                                          Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp6 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:density +
                                          Salt:density +
                                          Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp7 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + BZT_init:density +
                                          Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp8 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt + 
                                          Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp9 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + density + 
                                          BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp10 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + density + 
                                          BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model.temp11 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + 
                                           BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of benzotriazole acetyl-alanine
BZTacetylalanine_model <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                                     Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp1 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes + BZT_init:Microbes:density + 
                                     Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp2 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                           Salt:Microbes + Salt:density +
                                           Microbes:density + 
                                           BZT_init:Salt:Microbes + 
                                           Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp3 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                           Salt:Microbes + Salt:density +
                                           Microbes:density + 
                                           BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp4 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                           Salt:Microbes + Salt:density +
                                           Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp5 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + 
                                           Salt:Microbes + Salt:density +
                                           Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp6 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + 
                                           Salt:Microbes + 
                                           Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp7 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes + 
                                           Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp8 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp9 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + density + 
                                           BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp10 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                           BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model.temp11 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                            BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of aniline
aniline_model <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                            BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                            Salt:Microbes + Salt:density +
                            Microbes:density + 
                            BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                            Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model.temp1 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                            BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                            Salt:Microbes + Salt:density +
                            Microbes:density + 
                            BZT_init:Salt:Microbes + BZT_init:Salt:density + 
                            Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model.temp2 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                  Salt:Microbes + Salt:density +
                                  Microbes:density + 
                                  BZT_init:Salt:Microbes + 
                                  Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model.temp3 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                  Salt:Microbes + Salt:density +
                                  Microbes:density + 
                                  BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model.temp4 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                  Salt:Microbes + Salt:density +
                                  BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model.temp5 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + density + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                  Salt:Microbes + 
                                  BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of methylbenzotriazole
methylBZT_model <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density + 
                              BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                              Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model.temp1 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + density + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                              Salt:Microbes + Salt:density +
                              Microbes:density + 
                              BZT_init:Salt:density + BZT_init:Microbes:density + 
                              Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model.temp2 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + density + 
                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                    Salt:Microbes + Salt:density +
                                    Microbes:density + 
                                    BZT_init:Salt:density + BZT_init:Microbes:density + 
                                    Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of methoxybenzotriazole
methoxyBZT_model <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                               Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model.temp1 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of pthalic acid
pthalic_acid_model <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                               Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model.temp1 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + density + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                 Salt:Microbes + Salt:density +
                                 Microbes:density + 
                                 BZT_init:Salt:Microbes + BZT_init:Salt:density + 
                                 Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model.temp2 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + density + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                       Salt:Microbes + Salt:density +
                                       Microbes:density + 
                                       BZT_init:Salt:density + 
                                       Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model.temp3 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + density + 
                                       BZT_init:Salt + BZT_init:density +
                                       Salt:Microbes + Salt:density +
                                       Microbes:density + 
                                       BZT_init:Salt:density + 
                                       Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model.temp4 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + density + 
                                       BZT_init:Salt + BZT_init:density +
                                       Salt:Microbes + Salt:density +
                                       Microbes:density + 
                                       BZT_init:Salt:density + 
                                       Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of hydroxyBZT
hydroxyBZT_model <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                 Salt:Microbes + Salt:density +
                                 Microbes:density + 
                                 BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                                 Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp1 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:density + 
                               Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp2 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes +
                                     Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp3 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     Microbes:density + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp4 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density +
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp5 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp6 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                                     Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp7 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes +
                                     Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp8 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp9 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + density + 
                                     BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp10 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                     BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model.temp11 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Microbes + 
                                      BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for duckweed pixel area
px_model <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                               Salt:Microbes + Salt:density +
                               Microbes:density + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                               Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp1 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                       BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                       Salt:Microbes + Salt:density +
                       Microbes:density + 
                       BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp2 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density + 
                             BZT_init:Salt:Microbes + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp3 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density + 
                             BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp4 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp5 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp6 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp7 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp8 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Microbes + BZT_init:density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp9 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp10 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + density, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp11 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp12 <- MCMCglmm(scale(px.mn)~Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model.temp13 <- MCMCglmm(scale(px.mn)~Salt, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for optical density
od_model <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                       BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                       Salt:Microbes + Salt:density +
                       Microbes:density + 
                       BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Microbes:density + 
                       Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp1 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                       BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                       Salt:Microbes + Salt:density +
                       Microbes:density + 
                       BZT_init:Salt:Microbes + BZT_init:Microbes:density + 
                       Salt:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp2 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density + 
                             BZT_init:Salt:Microbes + BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp3 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density + 
                             BZT_init:Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp4 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:Microbes + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp5 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:density +
                             Salt:Microbes + Salt:density +
                             Microbes:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp6 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:density +
                             Salt:Microbes + Salt:density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp7 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + BZT_init:density +
                             Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp8 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             BZT_init:Salt + 
                             Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp9 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density + 
                             Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp10 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + density, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp11 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

od_model.temp12 <- MCMCglmm(scale(od.mn)~Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)


#linear model for optical density, restricting to inoculated wells
odi_model <- MCMCglmm(scale(od.mn)~BZT_init + Salt + density +
                       BZT_init:Salt + BZT_init:density +
                       Salt:density +
                       BZT_init:Salt:density, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp1 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + density +
                        BZT_init:Salt + BZT_init:density +
                        Salt:density, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp2 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + density +
                              BZT_init:density +
                              Salt:density, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp3 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + density +
                              Salt:density, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp4 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + density, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp5 <- MCMCglmm(scale(od.mn)~BZT_init + Salt, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model.temp6 <- MCMCglmm(scale(od.mn)~Salt, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

#linear model for percent decrease in benzotriazole concentration, with distance to city center as the location descriptor 
bzs_model_km <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                        Salt:Microbes + Salt:km +
                        Microbes:km + 
                        BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                        Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp1 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                           BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                           Salt:Microbes + Salt:km +
                           Microbes:km + 
                           BZT_init:Salt:Microbes + BZT_init:Salt:km + 
                           Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp2 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km + 
                                 BZT_init:Salt:Microbes + BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp3 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km + 
                                 BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp4 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp5 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + 
                                 Salt:Microbes + Salt:km +
                                 Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp6 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp7 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp8 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km +
                                 Salt:km +
                                 Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_km.temp9 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + km +
                                 Salt:km, random = ~ Genotype , data=bzs1,verbose=F)


#linear model for amount of benzotriazole alanine, with distance to city center as the location descriptor 
BZTalanine_model_km <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                               Salt:Microbes + Salt:km +
                               Microbes:km + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                               Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp1 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                  Salt:Microbes + Salt:km +
                                  Microbes:km + 
                                  BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp2 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp3 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp4 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + 
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp5 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + 
                                        Salt:Microbes + 
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp6 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + 
                                        Salt:Microbes +
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTalanine_model_km.temp7 <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                        BZT_init:Salt + BZT_init:Microbes + 
                                        Salt:Microbes +
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of glycosylated benzotriazole, with distance to city center as the location descriptor  
glycosylatedBZT_model_km <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                    Salt:Microbes + Salt:km +
                                    Microbes:km + 
                                    BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                                    Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp1 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                       Salt:Microbes + Salt:km +
                                       Microbes:km + 
                                       BZT_init:Salt:Microbes + BZT_init:Microbes:km + 
                                       Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp2 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                             Salt:Microbes + Salt:km +
                                             Microbes:km + 
                                             BZT_init:Salt:Microbes + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp3 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                             Salt:Microbes + Salt:km +
                                             Microbes:km + 
                                             BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp4 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                             Salt:Microbes + Salt:km +
                                             Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp5 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                             Salt:Microbes + 
                                             Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp6 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + 
                                             Salt:Microbes + 
                                             Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp7 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + BZT_init:Microbes + 
                                             Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp8 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt + 
                                             Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp9 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + km + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp10 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_km.temp11 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt +
                                              BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of benzotriazole acetyl-alanine, with distance to city center as the location descriptor 
BZTacetylalanine_model_km <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                     Salt:Microbes + Salt:km +
                                     Microbes:km + 
                                     BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                                     Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp1 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes + BZT_init:Microbes:km + 
                                        Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp2 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                              Salt:Microbes + Salt:km +
                                              Microbes:km + 
                                              BZT_init:Salt:Microbes + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp3 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                              Salt:Microbes + Salt:km +
                                              Microbes:km + 
                                              BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp4 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                              Salt:Microbes + Salt:km +
                                              Microbes:km, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp5 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + 
                                              Salt:Microbes + Salt:km +
                                              Microbes:km, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp6 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + 
                                              Salt:Microbes + 
                                              Microbes:km, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp7 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes + 
                                              Salt:Microbes, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp8 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp9 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + km + 
                                              BZT_init:Salt, random = ~ Genotype, data=bzs1,verbose=F)

BZTacetylalanine_model_km.temp10 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes +  
                                              BZT_init:Salt, random = ~ Genotype, data=bzs1,verbose=F)


#linear model for amount of aniline, with distance to city center as the location descriptor 
aniline_model_km <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + km + 
                            BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                            Salt:Microbes + Salt:km +
                            Microbes:km + 
                            BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                            Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model_km.temp1 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                               Salt:Microbes + Salt:km +
                               Microbes:km + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:km +  
                               Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)


aniline_model_km.temp2 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + km + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                     Salt:Microbes + Salt:km +
                                     Microbes:km + 
                                     BZT_init:Salt:Microbes + BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model_km.temp3 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + km + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                     Salt:Microbes + Salt:km +
                                     Microbes:km + 
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

aniline_model_km.temp4 <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + km + 
                                     BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                     Salt:Microbes + Salt:km +
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)


#linear model for amount of methylbenzotriazole, with distance to city center as the location descriptor 
methylBZT_model_km <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                              BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                              Salt:Microbes + Salt:km +
                              Microbes:km + 
                              BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                              Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp1 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km + 
                                 BZT_init:Salt:km + BZT_init:Microbes:km + 
                                 Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp2 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                       Salt:Microbes + Salt:km +
                                       Microbes:km + 
                                       BZT_init:Microbes:km + 
                                       Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)


methylBZT_model_km.temp3 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                       Salt:Microbes + Salt:km +
                                       Microbes:km + 
                                       BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp4 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                       Salt:Microbes + Salt:km +
                                       Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp5 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Microbes + BZT_init:km +
                                       Salt:Microbes + Salt:km +
                                       Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp6 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Microbes + BZT_init:km +
                                       Salt:km +
                                       Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp7 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:km +
                                       Salt:km +
                                       Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp8 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       Salt:km +
                                       Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp9 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km + 
                                       Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp10 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp11 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp12 <- MCMCglmm(scale(methylBZT)~Salt + km, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_km.temp13 <- MCMCglmm(scale(methylBZT)~km, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of methoxybenzotriazole, with distance to city center as the location descriptor 
methoxyBZT_model_km <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                               Salt:Microbes + Salt:km +
                               Microbes:km + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                               Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_km.temp1 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + km + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                  Salt:Microbes + Salt:km +
                                  Microbes:km + 
                                  BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of pthalic acid, with distance to city center as the location descriptor 
pthalic_acid_model_km <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                 BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                 Salt:Microbes + Salt:km +
                                 Microbes:km + 
                                 BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                                 Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp1 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                    Salt:Microbes + Salt:km +
                                    Microbes:km + 
                                    BZT_init:Salt:Microbes + BZT_init:Salt:km +
                                    Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp2 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                          BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                          Salt:Microbes + Salt:km +
                                          Microbes:km + 
                                          BZT_init:Salt:km +
                                          Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp3 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                          BZT_init:Salt + BZT_init:km +
                                          Salt:Microbes + Salt:km +
                                          Microbes:km + 
                                          BZT_init:Salt:km +
                                          Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp4 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                          BZT_init:Salt + BZT_init:km +
                                          Salt:Microbes + Salt:km +
                                          Microbes:km + 
                                          BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp5 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                          BZT_init:Salt + BZT_init:km +
                                          Salt:km +
                                          Microbes:km + 
                                          BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp6 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + km + 
                                          BZT_init:Salt + BZT_init:km +
                                          Salt:km +
                                          BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_km.temp7 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + km + 
                                          BZT_init:Salt + BZT_init:km +
                                          Salt:km +
                                          BZT_init:Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for amount of hydroxyBZT, with distance to city center as the location descriptor 
hydroxyBZT_model_km <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                               Salt:Microbes + Salt:km +
                               Microbes:km + 
                               BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                               Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp1 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                  BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                  Salt:Microbes + Salt:km +
                                  Microbes:km + 
                                  BZT_init:Salt:Microbes + BZT_init:Salt:km + 
                                  Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp2 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes + 
                                        Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp3 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km + 
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)


hydroxyBZT_model_km.temp4 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km +
                                        Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp5 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp6<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                        BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                        Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp7<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp8<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp9<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + km + 
                                       BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp10<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                       BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_km.temp11<- MCMCglmm(scale(hydroxyBZT)~BZT_init + Microbes + 
                                        BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for duckweed pixel area, with distance to city center as the location descriptor 
px_model_km <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                       Salt:Microbes + Salt:km +
                       Microbes:km + 
                       BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                       Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp1 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                          BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                          Salt:Microbes + Salt:km +
                          Microbes:km + 
                          BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp2 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km + 
                                BZT_init:Salt:Microbes + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp3 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km + 
                                BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp4 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp5 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp6 <- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp7<- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp8<- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Microbes + BZT_init:km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp9<- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km + 
                               BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp10<- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes + km, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp11<- MCMCglmm(scale(px.mn)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp12<- MCMCglmm(scale(px.mn)~Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

px_model_km.temp13<- MCMCglmm(scale(px.mn)~Salt, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for optical km, with distance to city center as the location descriptor 
od_model_km <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                       BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                       Salt:Microbes + Salt:km +
                       Microbes:km + 
                       BZT_init:Salt:Microbes + BZT_init:Salt:km + BZT_init:Microbes:km + 
                       Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp1 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                          BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                          Salt:Microbes + Salt:km +
                          Microbes:km + 
                          BZT_init:Salt:Microbes + BZT_init:Microbes:km + 
                          Salt:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp2 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km + 
                                BZT_init:Salt:Microbes + BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp3 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km + 
                                BZT_init:Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp4 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:Microbes + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp5 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:km +
                                Salt:Microbes + Salt:km +
                                Microbes:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp6 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:km +
                                Salt:Microbes + Salt:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp7 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:km +
                                Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp8 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:Salt + BZT_init:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp9 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km + 
                                BZT_init:km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp10 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + Microbes + km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp11 <- MCMCglmm(scale(od.mn)~Salt + Microbes + km, random = ~ Genotype , data=bzs1,verbose=F)

od_model_km.temp12 <- MCMCglmm(scale(od.mn)~Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#linear model for optical km, restricting to inoculated wells, with distance to city center as the location descriptor 
odi_model_km <- MCMCglmm(scale(od.mn)~BZT_init + Salt + km +
                        BZT_init:Salt + BZT_init:km +
                        Salt:km +
                        BZT_init:Salt:km, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp1 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + km +
                           BZT_init:Salt + BZT_init:km +
                           Salt:km, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp2 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + km +
                                 BZT_init:km +
                                 Salt:km, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp3 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + km +
                                 Salt:km, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp4 <- MCMCglmm(scale(od.mn)~BZT_init + Salt + km, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp5 <- MCMCglmm(scale(od.mn)~BZT_init + Salt, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

odi_model_km.temp6 <- MCMCglmm(scale(od.mn)~Salt, random = ~ Genotype , data=bzs1[bzs1$Microbes == 'Yes',],verbose=F)

#linear model for percent decrease in benzotriazole concentration, both location descriptors
#removed scaling to avoid error message
#some effects not estimable
# bzs_model_both <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                              Salt:Microbes + Salt:density + Salt:km +
#                              Microbes:density + Microbes:km + 
#                              density:km +
#                              BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                              BZT_init:Microbes:density + BZT_init:Microbes:km + BZT_init:density:km +
#                              Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
#                              Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp1 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                              Salt:Microbes + Salt:density + Salt:km +
#                              Microbes:density + Microbes:km + 
#                              density:km +
#                              BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                              BZT_init:Microbes:density + BZT_init:Microbes:km + BZT_init:density:km +
#                              Salt:Microbes:km + Salt:density:km +
#                              Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp2 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    Microbes:density + Microbes:km + 
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:Microbes:density + BZT_init:Microbes:km + BZT_init:density:km +
#                                    Salt:density:km +
#                                    Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp3 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    Microbes:density + Microbes:km + 
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:Microbes:density + BZT_init:density:km +
#                                    Salt:density:km +
#                                    Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp4 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    Microbes:density + Microbes:km + 
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km +
#                                    Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp5 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    Microbes:density +  
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km +
#                                    Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp6 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    Microbes:density +  
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp7 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    density:km +
#                                    BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp8 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    density:km +
#                                    BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp9 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:density + BZT_init:km +
#                                    Salt:Microbes + Salt:density + Salt:km +
#                                    density:km +
#                                    BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp10 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                    BZT_init:Salt + BZT_init:density + BZT_init:km +
#                                    Salt:density + Salt:km +
#                                    density:km +
#                                    BZT_init:Salt:density + BZT_init:Salt:km +
#                                    BZT_init:density:km +
#                                    Salt:density:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# bzs_model_both.temp11 <- MCMCglmm(BZT_percent_d~BZT_init + Salt + Microbes + density + km +
#                                     BZT_init:Salt + BZT_init:density + BZT_init:km +
#                                     Salt:density + Salt:km +
#                                     BZT_init:Salt:density + BZT_init:Salt:km, random = ~ Genotype, data=bzs1,verbose=F)
# 
# #linear model for amount of BZTalanine, both location descriptors
# #some effects not estimable
# BZTalanine_model_both <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + density + km +
#                              BZT_init:Salt + BZT_init:Microbes + BZT_init:density + BZT_init:km +
#                              Salt:Microbes + Salt:density + Salt:km +
#                              Microbes:density + Microbes:km + 
#                              density:km +
#                              BZT_init:Salt:Microbes + BZT_init:Salt:density + BZT_init:Salt:km +
#                              BZT_init:Microbes:density + BZT_init:Microbes:km + BZT_init:density:km +
#                              Salt:Microbes:density + Salt:Microbes:km + Salt:density:km +
#                              Microbes:density:km, random = ~ Genotype, data=bzs1,verbose=F)

#linear model for percent decrease in benzotriazole concentration, without location variable
bzs_model_noloc <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + 
                              BZT_init:Salt + BZT_init:Microbes +
                              Salt:Microbes +
                              BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_noloc.temp1 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + 
                              BZT_init:Salt + BZT_init:Microbes +
                              Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_noloc.temp2 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + 
                                    BZT_init:Salt +
                                    Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_noloc.temp3 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + 
                                    Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

bzs_model_noloc.temp4 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)


#significance of random effects, bzs model without location variable
bzs_model_DIC <- 0
bzs_model_norand_DIC <- 0
for(i in 1:10) {
  bzs_model_DIC[i] <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  bzs_model_norand_DIC[i] <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes, data=bzs1,verbose=F)$DIC
}

#linear model for amount of glycoslyated bzt, without location variable
glycosylatedBZT_model_noloc <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + 
                              BZT_init:Salt + BZT_init:Microbes +
                              Salt:Microbes +
                              BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_noloc.temp1 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + 
                                          BZT_init:Salt + BZT_init:Microbes +
                                          Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_noloc.temp2 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + 
                                                BZT_init:Salt + 
                                                Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_noloc.temp3 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + Microbes + 
                                                BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

glycosylatedBZT_model_noloc.temp4 <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + 
                                                BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, glycosylatedBZT model without location variable
glycosylatedBZT_model_DIC <- 0
glycosylatedBZT_model_norand_DIC <- 0
for(i in 1:10) {
  glycosylatedBZT_model_DIC[i] <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  glycosylatedBZT_model_norand_DIC[i] <- MCMCglmm(scale(glycosylatedBZT)~BZT_init + Salt + 
                                                    BZT_init:Salt,data=bzs1,verbose=F)$DIC
}

#linear model for amount of BZTalanine, without location variable
BZTalanine_model_noloc <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                          BZT_init:Salt + BZT_init:Microbes +
                                          Salt:Microbes +
                                          BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, BZTalanine model without location variable
BZTalanine_model_DIC <- 0
BZTalanine_model_norand_DIC <- 0
for(i in 1:10) {
  BZTalanine_model_DIC[i] <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                        BZT_init:Salt + BZT_init:Microbes +
                                        Salt:Microbes +
                                        BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  BZTalanine_model_norand_DIC[i] <- MCMCglmm(scale(BZTalanine)~BZT_init + Salt + Microbes + 
                                         BZT_init:Salt + BZT_init:Microbes +
                                         Salt:Microbes +
                                         BZT_init:Salt:Microbes, data=bzs1,verbose=F)$DIC
}

#linear model for amount of BZTacetylalanine, without location variable
BZTacetylalanine_model_noloc <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                          BZT_init:Salt + BZT_init:Microbes +
                                          Salt:Microbes +
                                          BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_noloc.temp1 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                           BZT_init:Salt + BZT_init:Microbes +
                                           Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_noloc.temp2 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                                 BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_noloc.temp3 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt + Microbes + 
                                                 BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

BZTacetylalanine_model_noloc.temp4 <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt +  
                                                 BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, BZTacetylalanine model without location variable
BZTacetylalanine_model_DIC <- 0
BZTacetylalanine_model_norand_DIC <- 0
for(i in 1:10) {
  BZTacetylalanine_model_DIC[i] <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt +  
                                              BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  BZTacetylalanine_model_norand_DIC[i] <- MCMCglmm(scale(BZTacetylalanine)~BZT_init + Salt +  
                                                     BZT_init:Salt, data=bzs1,verbose=F)$DIC
}


#linear model for amount of methylBZT, without location variable
methylBZT_model_noloc <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + 
                                           BZT_init:Salt + BZT_init:Microbes +
                                           Salt:Microbes +
                                           BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp1 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + 
                                    BZT_init:Salt + BZT_init:Microbes +
                                    Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp2 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + 
                                          BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp3 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes + 
                                          BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp4 <- MCMCglmm(scale(methylBZT)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp5 <- MCMCglmm(scale(methylBZT)~Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methylBZT_model_noloc.temp6 <- MCMCglmm(scale(methylBZT)~Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, methylBZT model without location variable
methylBZT_model_DIC <- 0
methylBZT_model_norand_DIC <- 0
for(i in 1:10) {
  methylBZT_model_DIC[i] <- MCMCglmm(scale(methylBZT)~Microbes, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  methylBZT_model_norand_DIC[i] <- MCMCglmm(scale(methylBZT)~Microbes, data=bzs1,verbose=F)$DIC
}

#linear model for amount of methoxyBZT, without location variable
methoxyBZT_model_noloc <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + 
                                    BZT_init:Salt + BZT_init:Microbes +
                                    Salt:Microbes +
                                    BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp1 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes +
                                     Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp2 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes + 
                                           BZT_init:Salt +
                                           Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp3 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes +
                                           Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp4 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt + Microbes, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp5 <- MCMCglmm(scale(methoxyBZT)~BZT_init + Salt, random = ~ Genotype , data=bzs1,verbose=F)

methoxyBZT_model_noloc.temp6 <- MCMCglmm(scale(methoxyBZT)~BZT_init, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, methoxyBZT model without location variable
methoxyBZT_model_DIC <- 0
methoxyBZT_model_norand_DIC <- 0
for(i in 1:10) {
  methoxyBZT_model_DIC[i] <- MCMCglmm(scale(methoxyBZT)~BZT_init, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  methoxyBZT_model_norand_DIC[i] <- MCMCglmm(scale(methoxyBZT)~BZT_init, data=bzs1,verbose=F)$DIC
}

#linear model for amount of aniline, without location variable
aniline_model_noloc <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes +
                                     Salt:Microbes +
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, aniline model without location variable
aniline_model_DIC <- 0
aniline_model_norand_DIC <- 0
for(i in 1:10) {
  aniline_model_DIC[i] <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes +
                                     Salt:Microbes +
                                     BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  aniline_model_norand_DIC[i] <- MCMCglmm(scale(aniline)~BZT_init + Salt + Microbes + 
                                            BZT_init:Salt + BZT_init:Microbes +
                                            Salt:Microbes +
                                            BZT_init:Salt:Microbes, data=bzs1,verbose=F)$DIC
}

#linear model for amount of pthalic_acid, without location variable
pthalic_acid_model_noloc <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + 
                                  BZT_init:Salt + BZT_init:Microbes +
                                  Salt:Microbes +
                                  BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_noloc.temp1 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + 
                                       BZT_init:Salt + BZT_init:Microbes +
                                       Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_noloc.temp2 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + 
                                             BZT_init:Salt + BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_noloc.temp3 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + Microbes + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_noloc.temp4 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

pthalic_acid_model_noloc.temp5 <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + 
                                             BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, pthalic_acid model without location variable
pthalic_acid_model_DIC <- 0
pthalic_acid_model_norand_DIC <- 0
for(i in 1:10) {
  pthalic_acid_model_DIC[i] <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + 
                                          BZT_init:Salt, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  pthalic_acid_model_norand_DIC[i] <- MCMCglmm(scale(pthalic_acid)~BZT_init + Salt + 
                                                 BZT_init:Salt, data=bzs1,verbose=F)$DIC
}


#linear model for amount of hydroxyBZT, without location variable
hydroxyBZT_model_noloc <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                       BZT_init:Salt + BZT_init:Microbes +
                                       Salt:Microbes +
                                       BZT_init:Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_noloc.temp1 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                     BZT_init:Salt + BZT_init:Microbes +
                                     Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_noloc.temp2 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                           BZT_init:Microbes +
                                           Salt:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_noloc.temp3 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt + Microbes + 
                                           BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_noloc.temp4 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + Salt +  
                                           BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

hydroxyBZT_model_noloc.temp5 <- MCMCglmm(scale(hydroxyBZT)~BZT_init + 
                                           BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)

#significance of random effects, hydroxyBZT model without location variable
hydroxyBZT_model_DIC <- 0
hydroxyBZT_model_norand_DIC <- 0
for(i in 1:10) {
  hydroxyBZT_model_DIC[i] <- MCMCglmm(scale(hydroxyBZT)~BZT_init + 
                                        BZT_init:Microbes, random = ~ Genotype , data=bzs1,verbose=F)$DIC
  hydroxyBZT_model_norand_DIC[i] <- MCMCglmm(scale(hydroxyBZT)~BZT_init + 
                                               BZT_init:Microbes, data=bzs1,verbose=F)$DIC
}

#MANOVA
res.man <- manova(cbind(BZT_percent_d, hydroxyBZT, BZTalanine, BZTacetylalanine, 
                        glycosylatedBZT, pthalic_acid, methoxyBZT, methylBZT, aniline) ~ BZT_init*Salt*Microbes, data = bzs1)
