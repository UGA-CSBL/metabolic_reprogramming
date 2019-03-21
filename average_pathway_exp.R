library(tidyverse)
library(DESeq2)

library(BiocParallel)
register(MulticoreParam(20))

# projects <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
#               "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA")
projects <- c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")

gene.list <- data.table::fread("~/storage/data/fenton/cleaned_gene_list.csv") %>%
  group_by(Pathway) %>%
  nest(Ensembl, Symbol, .key = Genes) %>%
  filter(str_detect(Pathway, "Fenton", negate = T))

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv") %>%
  mutate(tumor_stage = case_when(
    sample_type == "Solid Tissue Normal" ~ "Normal",
    str_detect(tumor_stage, "(^i$)|(\\si[abc]?$)|(1)") ~ "I",
    str_detect(tumor_stage, "(^ii$)|(\\si{2}[abc]?$)|(2)") ~ "II",
    str_detect(tumor_stage, "(^iii$)|(\\si{3}[abc]?$)|(3)") ~ "III",
    str_detect(tumor_stage, "(^iv$)|(\\siv[abc]?$)|(4)") ~ "IV",
    T ~ "Unknown"
  )) %>%
  filter(project %in% projects) %>%
  filter(tumor_stage != "Unknown") %>%
  select(-c(fileID, filename)) %>%
  as_tibble()

for (proj in projects) {
  print(str_glue("Working on project: {proj}"))
  design.mat <- clinical %>%
    filter(project == proj) %>%
    filter(tumor_stage != "Normal") %>%
    select(barcode, tumor_stage) %>%
    column_to_rownames("barcode") %>%
    mutate_if(is.character, as.factor)

  df <- data.table::fread(str_glue(
    "~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
    column_to_rownames("Ensembl") %>%
    head(-5) %>%
    select(one_of(rownames(design.mat)))

  dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = design.mat,
                                design= ~ tumor_stage)
  print(str_glue("\tEstimating size factors for {proj}"))
  dds <- estimateSizeFactors(dds)
  df.norm <- counts(dds, normalized = T)
  df.norm <- df.norm[rowSums(df.norm) >= 10, ]

  dfn <- df.norm %>%
    as.data.frame() %>%
    rownames_to_column("Ensembl") %>%
    mutate(Ensembl = str_extract(Ensembl, "^[^.]+"))

  # Get the most connected genes in each pathway across all stages
  mcgs <- tibble(Pathway = gene.list$Pathway, mcg = NA)
  for (i in seq(nrow(gene.list))) {
    pw.genes <- gene.list$Genes[[i]]$Ensembl
    pw.norm.counts <- dfn %>%
      filter(Ensembl %in% pw.genes) %>%
      column_to_rownames("Ensembl")

    if(nrow(pw.norm.counts) <= 2) {
      mcgs$mcg[i] <- rownames(pw.norm.counts)[1]
    } else {
      mcgs$mcg[i] <- cor(t(pw.norm.counts)) %>%
        as.data.frame() %>%
        rownames_to_column("Ensembl1") %>%
        gather(key = "Ensembl2", value = "Correlation", -Ensembl1) %>%
        filter(Ensembl1 != Ensembl2) %>%
        group_by(Ensembl1) %>%
        summarise(Connectivity = sum(Correlation)) %>%
        top_n(1, Connectivity) %>%
        .$Ensembl1
    }
  }

  clin.proj <- design.mat %>%
    rownames_to_column("Barcode")

  dfn <- dfn %>%
    filter(Ensembl %in% mcgs$mcg) %>%
    gather(key = "Barcode", value = "NormCounts", -Ensembl) %>%
    left_join(clin.proj, by = c("Barcode" = "Barcode")) %>%
    group_by(Ensembl, tumor_stage) %>%
    summarise(MeanNormCounts = mean(NormCounts)) %>%
    left_join(mcgs, by = c("Ensembl" = "mcg")) %>%
    select(Pathway, everything())
  data.table::fwrite(dfn, file = str_glue(
    "~/storage/data/metabolic_reprogramming/average_pathway_exp/{proj}.csv"))
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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] DESeq2_1.22.2               SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0          Biobase_2.42.0
# [7] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0         forcats_0.4.0
# [13] stringr_1.4.0               dplyr_0.8.0.1               purrr_0.3.1                 readr_1.3.1                 tidyr_0.8.3                 tibble_2.0.1
# [19] ggplot2_3.1.0               tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-137           bitops_1.0-6           lubridate_1.7.4        bit64_0.9-7            RColorBrewer_1.1-2     httr_1.4.0             tools_3.5.2            backports_1.1.3
# [9] R6_2.4.0               rpart_4.1-13           Hmisc_4.2-0            DBI_1.0.0              lazyeval_0.2.1         colorspace_1.4-0       nnet_7.3-12            withr_2.1.2
# [17] tidyselect_0.2.5       gridExtra_2.3          bit_1.1-14             compiler_3.5.2         cli_1.0.1              rvest_0.3.2            htmlTable_1.13.1       xml2_1.2.0
# [25] scales_1.0.0           checkmate_1.9.1        genefilter_1.64.0      digest_0.6.18          foreign_0.8-71         XVector_0.22.0         base64enc_0.1-3        pkgconfig_2.0.2
# [33] htmltools_0.3.6        htmlwidgets_1.3        rlang_0.3.1            readxl_1.3.0           rstudioapi_0.9.0       RSQLite_2.1.1          generics_0.0.2         jsonlite_1.6
# [41] acepack_1.4.1          RCurl_1.95-4.12        magrittr_1.5           GenomeInfoDbData_1.2.0 Formula_1.2-3          Matrix_1.2-16          Rcpp_1.0.0             munsell_0.5.0
# [49] stringi_1.4.3          yaml_2.2.0             zlibbioc_1.28.0        plyr_1.8.4             blob_1.1.1             grid_3.5.2             crayon_1.3.4           lattice_0.20-38
# [57] haven_2.1.0            splines_3.5.2          annotate_1.60.1        hms_0.4.2              locfit_1.5-9.1         knitr_1.22             pillar_1.3.1           geneplotter_1.60.0
# [65] XML_3.98-1.19          glue_1.3.1             latticeExtra_0.6-28    data.table_1.12.0      modelr_0.1.4           cellranger_1.1.0       gtable_0.2.0           assertthat_0.2.0
# [73] xfun_0.5               xtable_1.8-3           broom_0.5.1            survival_2.43-3        memoise_1.1.0          AnnotationDbi_1.44.0   cluster_2.0.7-1
