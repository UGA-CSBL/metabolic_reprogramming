library(tidyverse)
library(ggrepel)

# projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
#               "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
projects <- c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")

plot.avg.exp <- function(df, alpha.level = 0.9, suffix) {
  ymax <- ceiling(max(df$MeanNormCounts))
  p <- df %>%
    ggplot(aes(x = tumor_stage, y = MeanNormCounts,
               group = Pathway, color = Pathway, label = Pathway))+
    geom_point()+
    geom_line()+
    ggtitle(proj)+
    labs(
      x = "Tumor Stage",
      y = "log(Normalized Counts)"
    )+
    # Plot labels in the right margin
    coord_cartesian(clip = 'off')+
    geom_label_repel(
      aes(
        label = ifelse(tumor_stage == "IV", Pathway, ''),
        fill = Pathway, color = Pathway, alpha = alpha.level
      ),
      fontface = 'bold', color = 'white',
      segment.color = 'grey50',
      direction = "y",
      xlim = c(4.2, 20),
      ylim = c(NA, ymax)
    )+
    theme_minimal()+
    theme(
      legend.position = "none",
      plot.margin = unit(c(1,7,0,0), "in")
    )

  ggsave(p, file = str_glue(
    "~/storage/data/metabolic_reprogramming/average_pathway_exp/plot/{proj}_{suffix}.png"),
    device = "png", width = 20, height = 20,
    units = "in", dpi = "retina")
}


for (proj in projects) {
  df <- data.table::fread(str_glue(
    "~/storage/data/metabolic_reprogramming/average_pathway_exp/{proj}.csv"
  )) %>%
    mutate(MeanNormCounts = log(MeanNormCounts))
  pathways <- sort(unique(df$Pathway))
  middle.index <- ceiling(length(pathways) / 2)

  pw1 <- pathways[1:middle.index]
  pw2 <- pathways[(middle.index+1):length(pathways)]

  plot.avg.exp(df[df$Pathway %in% pw1, ], suffix = "1", alpha.level = 1)
  plot.avg.exp(df[df$Pathway %in% pw2, ], suffix = "2", alpha.level = 1)
}

# sessionInfo()
# R version 3.5.2 (2018-12-20)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
#
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
#
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggrepel_0.8.0   forcats_0.4.0   stringr_1.4.0   dplyr_0.8.0.1   purrr_0.3.1     readr_1.3.1     tidyr_0.8.3     tibble_2.0.1    ggplot2_3.1.0   tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.0        cellranger_1.1.0  pillar_1.3.1      compiler_3.5.2    plyr_1.8.4        tools_3.5.2       jsonlite_1.6      lubridate_1.7.4   gtable_0.2.0      nlme_3.1-137
# [11] lattice_0.20-38   pkgconfig_2.0.2   rlang_0.3.1       cli_1.0.1         rstudioapi_0.9.0  yaml_2.2.0        haven_2.1.0       withr_2.1.2       xml2_1.2.0        httr_1.4.0
# [21] generics_0.0.2    hms_0.4.2         grid_3.5.2        tidyselect_0.2.5  glue_1.3.1        R6_2.4.0          readxl_1.3.0      sessioninfo_1.1.1 modelr_0.1.4      magrittr_1.5
# [31] backports_1.1.3   scales_1.0.0      rvest_0.3.2       assertthat_0.2.0  colorspace_1.4-0  stringi_1.4.3     lazyeval_0.2.1    munsell_0.5.0     broom_0.5.1       crayon_1.3.4
