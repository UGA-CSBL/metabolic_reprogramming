library(tidyverse)
library(scales)
library(ggrepel)

# Cancer types, pathway gene lists and tumor staging information -------------
projects <- c(
  c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA"),
  c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")
)

# proteasomes <- data.table::fread(
#   "~/storage/data/metabolic_reprogramming/proteasome.csv") %>%
#   .$`Approved symbol`
#
# proteasomes <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
#   filter(external_gene_name %in% proteasomes) %>%
#   mutate(Pathway = "Proteasome") %>%
#   select(Pathway, Symbol = external_gene_name, Ensembl = ensembl_gene_id)

annot <- data.table::fread(
  "~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv"
) %>%
  select(Ensembl = ensembl_gene_id, GeneName = external_gene_name)

gangliosides <- tibble(
  Pathway = "Ganglioside biosynthesis",
  Symbol = c("B4GALT6", "UGCG", "ST6GAL1", "ST3GAL5", "ST3GAL2",
             "B4GALNT1", "B3GALT4", "ST3GAL3", "ST8SIA1")
) %>%
  left_join(annot, by = c("Symbol" = "GeneName"))

sialic.acids <- tibble(
  Pathway = "CMP-sialic acid biosynthesis",
  Symbol = c("GNE", "NANS", "NANP", "CMAS")
) %>%
  left_join(annot, by = c("Symbol" = "GeneName"))

gene.list <- data.table::fread("~/storage/data/fenton/cleaned_gene_list.csv") %>%
  filter(str_detect(Pathway, "Fenton", negate = T)) %>%
  select(-PathwayType) %>%
  # bind_rows(proteasomes) %>%
  bind_rows(gangliosides, sialic.acids) %>%
  # Fix multiple mappings in symbol -> Ensembl
  left_join(annot, by = "Ensembl") %>%
  drop_na(GeneName) %>%
  distinct() %>%
  mutate(
    Pathway = str_replace_all(Pathway, "[\\s/]", "_"),
    Symbol = case_when(
      Symbol != GeneName ~ GeneName,
      T ~ Symbol
    )) %>%
  select(-GeneName)

gene.mapping <- gene.list %>%
  select(-Pathway) %>%
  distinct()

# https://seer.cancer.gov/csr/1975_2015/browse_csr.php?sectionSEL=1&pageSEL=sect_01_table.04
# https://www.nature.com/articles/s41416-018-0140-8
rel.survival <- data.table::fread(
  "~/storage/data/metabolic_reprogramming/tcga_survival.csv"
)
colnames(rel.survival) <- c("Project", "TCGA", "SEER")
rel.survival <- rel.survival %>%
  mutate(
    Project = str_glue("TCGA-{Project}"),
    TCGA = str_extract(TCGA, "\\([\\d.]+"),
    SEER = str_extract(SEER, "\\([\\d.]+")
  ) %>%
  mutate(
    TCGA = as.numeric(str_remove(TCGA, "\\(")),
    SEER = as.numeric(str_remove(SEER, "\\("))
  ) %>%
  filter(Project %in% projects) %>%
  arrange(TCGA, SEER)

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
  # filter(!tumor_stage %in% c("Unknown", "Normal")) %>%
  filter(tumor_stage != "Unknown") %>%
  select(project, tumor_stage, barcode) %>%
  as_tibble()

FPKMtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}

# Get DEA results for each cancer type and stage -----------------------------
tumor.stages <- c("I", "II", "III", "IV")
res.dea <- tibble()
for (proj in projects) {
  for (t.stage in tumor.stages) {
    try(df <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/csv/{proj}_{t.stage}_vs_N.csv"
    )), silent = T)
    try(df <- data.table::fread(str_glue(
      "~/storage/data/TCGA/DEA/csv/all_stages_proj/{proj}_{t.stage}_vs_N.csv"
    )), silent = T)
    df <- df %>%
      rename(Ensembl = V1) %>%
      mutate(Ensembl = str_extract(Ensembl, "^[^.]+")) %>%
      inner_join(gene.mapping, by = "Ensembl") %>%
      filter(padj <= 0.01) %>%
      arrange(-abs(log2FoldChange)) %>%
      mutate(Project = proj, Stage = t.stage) %>%
      select(Project, Stage, Symbol, log2FoldChange)
    res.dea <- bind_rows(res.dea, df)
  }
}
# Get TPMs -------------------------------------------------------------------
res.exp <- tibble()
for (proj in projects) {
  message(str_glue("Reading project {proj}."))
  df <- data.table::fread(
    str_glue("/home/yi/CSBL_shared/RNASeq/TCGA/FPKM/{proj}.FPKM.csv")
  ) %>%
    mutate_if(is.numeric, FPKMtoTPM) %>%
    # Ensembl ID to gene symbol
    mutate(Ensembl = gsub("\\.\\d+$", "", Ensembl)) %>%
    filter(Ensembl %in% gene.list$Ensembl) %>%
    select(one_of(c("Ensembl", clinical$barcode[clinical$project == proj]))) %>%
    gather(key = "Barcode", value = "TPM", -Ensembl) %>%
    left_join(clinical, by = c("Barcode" = "barcode"))

  sample.num <- df %>%
    group_by(Ensembl) %>%
    tally()
  sample.num <- sample.num$n[1]
  # filter out lowly-expressed genes: <0.5 in more than half the samples
  genes.to.keep <- df %>%
    group_by(Ensembl) %>%
    add_tally(TPM >= 0.5) %>%
    filter(n > (sample.num / 2)) %>%
    distinct(Ensembl) %>%
    .$Ensembl

  df <- df %>%
    filter(Ensembl %in% genes.to.keep) %>%
    group_by(Ensembl, tumor_stage) %>%
    summarise(MeanTPM = mean(TPM)) %>%
    left_join(gene.list, by = "Ensembl") %>%
    mutate(Project = proj)
  res.exp <- bind_rows(res.exp, df)
}

res <- res.exp %>%
  left_join(res.dea, by = c("Project" = "Project",
                            "tumor_stage" = "Stage",
                            "Symbol" = "Symbol")) %>%
  mutate(
    tumor_stage = factor(tumor_stage,
                         levels = c("Normal", "I", "II", "III", "IV")),
    Project = factor(Project, levels = rel.survival$Project)
  )
rm(annot, res.exp, res.dea, df, genes.to.keep, t.stage, tumor.stages)
pw.list <- unique(gene.list$Pathway)

# Proteasome differential expression -----------------------------------------
# pw.list <- pw.list[pw.list != "Proteasome"]
# Proteasome.SigDiff <- res %>%
#   filter(Pathway == "Proteasome") %>%
#   filter(!is.na(log2FoldChange)) %>%
#   group_by(Project, Ensembl) %>%
#   summarise(DiffExp = max(log2FoldChange) - min(log2FoldChange)) %>%
#   ungroup() %>%
#   filter(DiffExp >= quantile(DiffExp)[4]) %>%
#   left_join(gene.mapping, by = "Ensembl") %>%
#   mutate(ProjSymbol = paste0(Project, "_", Symbol)) %>%
#   distinct(ProjSymbol, .keep_all = T)
# Proteasome.p <- res %>%
#   filter(Pathway == "Proteasome") %>%
#   mutate(
#     ProjSymbol = paste0(Project, "_", Symbol),
#     LineColor = case_when(
#       ProjSymbol %in% Proteasome.SigDiff$ProjSymbol ~ "blue",
#       T ~ "grey90"
#     ))

# Plot for each pathway ------------------------------------------------------
for (pw in pw.list) {
  # Filter for genes that have harmonic expression across stages
  SigDiff <- res %>%
    filter(Pathway == pw) %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(Project, Ensembl) %>%
    summarise(DiffExp = max(log2FoldChange) - min(log2FoldChange)) %>%
    ungroup() %>%
    filter(DiffExp >= quantile(DiffExp)[4]) %>%
    left_join(gene.mapping, by = "Ensembl") %>%
    mutate(ProjSymbol = paste0(Project, "_", Symbol)) %>%
    distinct(ProjSymbol, .keep_all = T)

  p <- res %>%
    filter(Pathway == pw) %>%
    mutate(
      ProjSymbol = paste0(Project, "_", Symbol),
      LineColor = case_when(
        ProjSymbol %in% SigDiff$ProjSymbol ~ "red",
        T ~ "grey"
      ),
      PointColor = case_when(
        is.na(log2FoldChange) ~ "grey",
        T ~ "red"
      )) %>%
    # bind_rows(Proteasome.p) %>%
    # mutate(LineColor = factor(LineColor,
    #                           levels = c("grey", "red", "blue", "grey90"))) %>%
    mutate(
      LineColor = factor(LineColor,
                         levels = c("grey", "red")),
      PointColor= factor(PointColor,
                         levels = c("grey", "red"))) %>%
    ggplot(aes(x = tumor_stage, y = MeanTPM))+
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x, n = 10),
                       labels = trans_format("log2", math_format(2^.x)))+
    geom_point(aes(group = Project, color = PointColor))+
    geom_line(aes(group = interaction(Project, Symbol), color = LineColor))+
    facet_wrap(~Project, nrow = 3, ncol = 5)+
    # Label genes in the right margin
    geom_label_repel(
      aes(
        label = ifelse(
          (tumor_stage == "IV") & (str_detect(LineColor, "grey", negate = T)),
          Symbol, ''),
        fill = Symbol
      ),
      fontface = 'bold', color = 'white',
      segment.color = 'grey50',
      direction = "y",
      xlim = c(5, 15)
    )+
    # Label log2 Fold Change for red lines
    geom_text_repel(
      aes(label = ifelse((LineColor == "red") & (!is.na(log2FoldChange)),
                         round(log2FoldChange, 2), "")),
      # segment.color = 'grey50',
      size = 2,
      colour = "#0392ce",
      max.iter = 1000
    )+
    ggtitle(pw)+
    labs(
      x = "Tumor Stage",
      y = "Average TPM"
    )+
    # Plot labels in the margins
    coord_cartesian(xlim = c(1.5, 4.5), clip = 'off')+
    theme_minimal()+
    theme(
      legend.position = "none",
      panel.spacing = unit(1, "in"),
      plot.margin = unit(c(0, 1, 0, 0), units = "in")
    )+
    # scale_colour_manual(values=c("gray60", "red", "blue", "gray90"))
    scale_colour_manual(values=c("gray80", "red"))

  ggsave(p, file = str_glue(
    "~/storage/data/metabolic_reprogramming/average_pathway_exp/plot/log2FC/{pw}.png"),
    device = "png", width = 25, height = 15,
    units = "in", dpi = "retina")
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] ggrepel_0.8.0   scales_1.0.0    forcats_0.4.0   stringr_1.4.0   dplyr_0.8.0.1   purrr_0.3.1     readr_1.3.1     tidyr_0.8.3     tibble_2.0.1    ggplot2_3.1.0   tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.0        cellranger_1.1.0  pillar_1.3.1      compiler_3.5.2    plyr_1.8.4        tools_3.5.2       jsonlite_1.6      lubridate_1.7.4   gtable_0.2.0      nlme_3.1-137      lattice_0.20-38   pkgconfig_2.0.2   rlang_0.3.1
# [14] cli_1.0.1         rstudioapi_0.9.0  yaml_2.2.0        haven_2.1.0       withr_2.1.2       xml2_1.2.0        httr_1.4.0        generics_0.0.2    hms_0.4.2         grid_3.5.2        tidyselect_0.2.5  glue_1.3.1        data.table_1.12.0
# [27] R6_2.4.0          readxl_1.3.0      modelr_0.1.4      magrittr_1.5      backports_1.1.3   rvest_0.3.2       assertthat_0.2.0  colorspace_1.4-0  stringi_1.4.3     lazyeval_0.2.1    munsell_0.5.0     broom_0.5.1       crayon_1.3.4
