# load libraries
library("ggpubr")
library(ComplexUpset)
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(plyr)
library(scales)
library(forcats)
library(rmarkdown)
library(BiocParallel)
library(dplyr)
library(edgeR)
library(limma)
library(ggrepel)
library(ggplot2)
library(gplots)
library(grDevices)
require(philentropy)
library(rtracklayer)
library(stringr)
require(variancePartition)
library(reshape)
library(Glimma)
library(plyr)
library(corrplot)
library(ggpubr)
library(caret)
library(glmnet)
library(vroom)
library(matrixStats)
library("data.table")
library(DESeq2)
library(dittoSeq)
library(Hmisc)
library(tidyr)
library(gridExtra)
library(grid)
require(openxlsx)
library(mvIC)
library(RColorBrewer)
library(devtools)
library(reshape2)
library(edgeR)
library(tidyr)
library(limma)
library(pheatmap)
library(tidyverse, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")
library(tidyr, lib.loc = "/usr/local/biotools/rpackages/R-4.1.2-2021-11-11")
library(vctrs, lib.loc = "/usr/local/biotools/rpackages/R-4.2.2-2023-02-01")

# paths, colors, shapes and more
#LBD <- "LBD"
#AD <- "AD"
#PA <- "PA"
#CONTROL <- "CONTROL"
#control_color <- "#4682B0"
#AD_color <- "#B4464B"
#PA_color <- "#B4AF46"
#LBD_color <- "gray35"
#control_shape <- c(15) # square
#AD_shape <- c(16) # circle
#PA_shape <- c(17) # triangle
#LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B0", "gray35")
ATSColors <- c("#4682B0",  "gray35", "gray65", "gray", "gray85")

SexColors <- c("#490092", "#D55E00")
colorbindColors <- dittoColors()
correlationColors <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

input_path = c("/research/labs/neurology/fryer/m239830/LBD_CWOW/proteomics/LBD_CWOW_proteomics/counts/")

saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}


metadata <-
  read.delim(
    "/research/labs/neurology/fryer/projects/LBD_CWOW/proteomics/counts/CWOW_metadata_430.tsv"
  )
metadata$APOE.1 <- NULL
metadata$Note <- NULL
metadata$Cing.LB[is.na(metadata$Cing.LB)] <- 0
metadata$Braak.NFT[is.na(metadata$Braak.NFT)] <- 0
metadata$Thal.amyloid[is.na(metadata$Thal.amyloid)] <- 0

# set factor levels
metadata$TYPE <-
  factor(metadata$TYPE, levels = c("CONTROL", "LBD"))

df <- as.data.frame(metadata)
df$group <- df$TYPE

df$group.redefined <-
  ifelse(
    df$group == "LBD" &
      df$Braak.NFT <= 3 & df$Thal.amyloid < 2,
    "LBD_S",
    ifelse(
      df$group == "LBD" &
        df$Braak.NFT > 3 & df$Thal.amyloid >= 2,
      "LBD_ATS",
      ifelse(
        df$group == "LBD" &
          df$Braak.NFT <= 3 & df$Thal.amyloid >= 2,
        "LBD_AS",
        ifelse(
          df$group == "LBD" &
            df$Braak.NFT > 3 & df$Thal.amyloid < 2,
          "LBD_TS",
          ifelse(
            df$group == "CONTROL" &
              df$Braak.NFT <= 3 & df$Thal.amyloid < 2,
            "CONTROL",
            ifelse(
              df$group == "AD" & df$Braak.NFT > 3 & df$Thal.amyloid >= 3,
              "AD",
              ifelse(
                df$group == "PA" &
                  df$Braak.NFT <= 3 & df$Thal.amyloid >= 2,
                "PA",
                "OTHER"
              )
            )
          )
        )
      )
    )
  )

# check
table(df$group)
table(df$group.redefined)


df <- df %>%
  filter(!(group.redefined %in% c("OTHER", "NA")))

rm(metadata)
metadata <- df
rm(df)

table(metadata$group.redefined)
metadata$ATS  <- metadata$group.redefined
metadata$ATS <-
  factor(metadata$ATS, levels = c("CONTROL", "LBD_ATS", "LBD_AS", "LBD_TS", "LBD_S"))
