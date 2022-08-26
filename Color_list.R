##############################
#
# Color list for heatmaps
#
# Written by Beáta Szeitz
# Project: SCLC proteomics
#
##############################

colorlist <- list("Subtype" = 
                    c("SCLC-A"="#0000FE",
                      "SCLC-Y"="#FB02FE",
                      "SCLC-N"="#21FE06",
                      "SCLC-P"="#FEFE0A"),
                  "ASCL1 qPCR" = 
                    circlize::colorRamp2(c(0,0.06,1), 
                                         c("white","#0098fe","#0000FE")),
                  "NEUROD1 qPCR" = 
                    circlize::colorRamp2(c(0,0.11,1), 
                                         c("white","#97fe06","#21FE06")),
                  "POU2F3 qPCR" = 
                    circlize::colorRamp2(c(0,0.03,1), 
                                         c("white","#FEFE0A","#fedd0a")),
                  "YAP1 qPCR" = 
                    circlize::colorRamp2(c(0,0.20,1), 
                                         c("white","#FB02FE","#d402fe")),
                  "ASCL1.transcript" = 
                    circlize::colorRamp2(c(0,0.06,1), 
                                         c("white","#0098fe","#0000FE")),
                  "NEUROD1.transcript" = 
                    circlize::colorRamp2(c(0,0.11,1), 
                                         c("white","#97fe06","#21FE06")),
                  "POU2F3.transcript" = 
                    circlize::colorRamp2(c(0,0.03,1), 
                                         c("white","#FEFE0A","#fedd0a")),
                  "YAP1.transcript" = 
                    circlize::colorRamp2(c(0,0.20,1), 
                                         c("white","#FB02FE","#d402fe")),
                  "Cell line origin" =
                    c("lung"="#3f8098",
                      "pleural.effusion"="lightblue",
                      "metastatic"="#ffc3a0"),
                  "Chemo" = 
                    c("chemo-naive" ="white",
                      "post-chemo" = "brown"),
                  "Culture type" = 
                    c("Adherent"="#23b866",
                      "Semi-adherent"="#ef5675",
                      "Suspension"="#7a5195"),
                  "MSbatch" = 
                    c("1" = "white",
                      "2" = "black"),
                  "Data.Acquisiton" = 
                    c("DDA" = "white",
                      "DIA" = "black"),
                  "TP53 mutation status" = 
                    c("WT" = "white",
                      "Mut (non-deleterious)" = "#66bca4",
                      "Mut (deleterious)" = "#b66f6f"),
                  "RB1 mutation status" = 
                    c("WT" = "white",
                      "Mut (non-deleterious)" = "#66bca4",
                      "Mut (deleterious)" = "#b66f6f"))
