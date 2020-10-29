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

bzs_model.temp10 <- MCMCglmm(scale(BZT_percent_d)~BZT_init + Salt + Microbes + density + 
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
