## code to prepare sample metadata goes here
library(magrittr)

url_sample_attributes <-
  paste0("https://storage.googleapis.com/gtex_analysis_v8/annotations",
         "/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

download.file(url = url_sample_attributes,
              destfile = file.path(
                "data-raw",
                "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
              quiet = T)

sampleMetadata <- as.data.frame(
  readr::read_tsv(
    file = file.path(
      "data-raw",
      "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
    col_names = T, show_col_types = F) %>%
    dplyr::select(SAMPID,
                  SMTS,
                  SMTSD,
                  SMAFRZE,
                  SMATSSCR) %>%
    dplyr::rename(sample_id = SAMPID,
                  tissue_type = SMTS,
                  tissue_type_detailed = SMTSD,
                  in_analysis_freeze = SMAFRZE,
                  autolysis_score = SMATSSCR) %>%
    dplyr::filter(in_analysis_freeze == "RNASEQ") %>%
    dplyr::mutate(tissue_site_detail_id =
                    stringr::str_replace(
                      stringr::str_replace_all(
                        tissue_type_detailed," - |\\)| \\(| ","_"),
                      "_$","")
    ) %>%
    dplyr::mutate(tcga_study_match = dplyr::case_when(
      tissue_type_detailed == "Skin - Sun Exposed (Lower leg)" ~ "SKCM",
      tissue_type_detailed == "Pancreas" ~ "PAAD",
      tissue_type_detailed == "Prostate" ~ "PRAD",
      tissue_type_detailed == "Adrenal Gland" ~ "ACC",
      tissue_type_detailed == "Thyroid" ~ "THCA",
      tissue_type_detailed == "Ovary" ~ "OV",
      tissue_type_detailed == "Lung" ~ "LUAD|LUSC",
      tissue_type_detailed == "Stomach" ~ "STAD",
      tissue_type_detailed == "Testis" ~ "TGCT",
      tissue_type_detailed == "Breast - Mammary Tissue" ~ "BRCA",
      tissue_type_detailed == "Liver" ~ "LIHC",
      tissue_type_detailed == "Bladder" ~ "BLCA",
      tissue_type_detailed == "Cervix - Ectocervix" ~ "CESC",
      tissue_type_detailed == "Cervix - Endocervix" ~ "CESC",
      tissue_type_detailed == "Esophagus - Muscularis" ~ "ESCA",
      tissue_type_detailed == "Esophagus - Mucosa" ~ "ESCA",
      tissue_type_detailed == "Esophagus - Gastroesophageal Junction" ~ "ESCA",
      tissue_type_detailed == "Colon - Transverse" ~ "COAD",
      tissue_type_detailed == "Colon - Sigmoid" ~ "COAD",
      TRUE ~ as.character(NA)
    ))
)

num_samples_per_tissuetype <- as.data.frame(
  sampleMetadata %>%
  dplyr::group_by(tissue_site_detail_id) %>%
  dplyr::summarise(num_tissue_samples = dplyr::n())
)

sampleMetadata <- sampleMetadata %>%
  dplyr::left_join(num_samples_per_tissuetype) %>%
  dplyr::mutate(sample_type = "normal/healthy")


usethis::use_data(sampleMetadata, overwrite = T)
