library(TCGAbiolinks)
library(xml2)
library(tidyverse)

setwd("~/storage/data/TCGA/clinical/TCGAbiolinks")

# Download xml clinical data files -------------------------------------------
projects <- c(
  c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA"),
  c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")
)
for (proj in projects) {
  query <- GDCquery(project = proj,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
}

# Parse xml files for each cancer type ---------------------------------------
projects <- list.files("./GDCdata/")

for (proj in projects.downloaded) {
  patients <- list.files(
    stringr::str_glue("./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/"),
    pattern = "xml$", recursive = T)
  filenames <- stringr::str_glue(
    "./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/{patients}"
  )
  patients <- map_chr(str_split(patients, "/"), function(x) x[[1]])

  res <- NULL
  for (i in seq_along(filenames)) {
    df <- read_xml(filenames[i])
    barcode <- df %>%
      xml_find_first("//shared:bcr_patient_barcode") %>%
      xml_text()

    # In project:patient, look for all fields with text
    df <- xml_find_all(df, "//*[local-name(.) = 'patient']//*[text()]")
    df <- tibble(
      fileID = patients[i],
      field = map_chr(df, xml_path),
      value = map_chr(df, xml_text)
    ) %>%
      mutate(field = str_extract(field, "[^:]+$"))
    res <- bind_rows(res, df)
  }
  fields <- c("fileID", unique(res$field))
  res <- res %>%
    distinct() %>%
    group_by(fileID, field) %>%
    summarise(value = paste(value, collapse = "|")) %>%
    ungroup() %>%
    spread(key = field, value = value) %>%
    select(one_of(fields))
  data.table::fwrite(res, str_glue("./{proj}.csv"))
}
