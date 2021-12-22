#' Download compressed GTEx expression dataset (.gct.gz)
#'
#' @param download_dir target local directory for download file
#' @param overwrite do not overwrite (FALSE) if expression file exist, or overwrite (TRUE)
#' @param gtex_version version of gtex database
#'
#' @keywords internal
#'
#'
download_raw_data <- function(
  download_dir = NULL,
  overwrite = F,
  gtex_version = "v8"){

  logger <- log4r::logger(
    threshold = "INFO",
    appenders = log4r::console_appender(log4r_layout))


  if(!is.null(download_dir)){
    val <- assertthat::validate_that(
      dir.exists(download_dir)
    )
    if(!is.logical(val)){
      log4r_info(paste0("ERROR: ",val))
      return()
    }

  }

  url_rnaseq_gene_tpm <-
    paste0("https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data",
           "/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

  log4r_info(logger, paste0("Starting download of compressed GTEx expression data"))
  log4r_info(logger, paste0("Remote URL: ",
                            url_rnaseq_gene_tpm))


  if(!RCurl::url.exists(url_rnaseq_gene_tpm)){
    log4r_info(paste0("Download error: ",url_rnaseq_gene_tpm," does not exist!"))
    return()
  }


  if(!file.exists(file.path(
    download_dir,
    paste0("gtex_", gtex_version,
           "_gene_rnaseq_tpm.gct.gz")
  )) | overwrite == T){

    utils::download.file(
      url = url_rnaseq_gene_tpm,
      destfile = file.path(
        download_dir,
        paste0("gtex_", gtex_version,
               "_gene_rnaseq_tpm.gct.gz")),
      quiet = T
    )

    log4r_info(logger, paste0("Completing download of compressed GTEx expression data"))


  }else{
    log4r_info(logger, paste0("Compressed GTEx expression data already exist - no downloading",
                              " (set 'overwrite = T' to force download)"))

  }

  return()

}



#' Function that retrieves GTex RNA-seq sample identifiers
#' for a given tissue type
#'
#' @param tissue Type of tissue. Either of the following: 'Blood','Brain',
#' 'Adipose Tissue','Muscle','Blood Vessel','Heart','Ovary','Uterus','Vagina',
#' 'Breast','Skin','Salivary Gland','Adrenal Gland','Thyroid','Lung','Spleen',
#' 'Pancreas','Espohagus','Stomach','Colon','Small Intestine','Prostate',
#' 'Testis','Nerve','Pituitary','Liver','Kidney','Cervix Uteri','Fallopian Tube',
#' 'Bladder','Bone Marrow'
#' @param analysis RNASEQ
#' @param max_autolysis_score Set maximum allowed autolysis score for samples
#' retrieved (1 (None) to 3 (Severe))
#'
#' @keywords internal
#'
get_sample_identifiers <- function(
  tissue = "Adrenal Gland",
  analysis = "RNASEQ",
  max_autolysis_score = 3){

  logger <- log4r::logger(
    threshold = "INFO",
    appenders = log4r::console_appender(log4r_layout))


  tissue_types_gtex <-
    c("Blood", "Brain", "Adipose Tissue",
      "Muscle", "Blood Vessel", "Heart",
      "Ovary", "Uterus", "Vagina",
      "Breast", "Skin", "Salivary Gland",
      "Adrenal Gland", "Thyroid",
      "Lung", "Spleen", "Pancreas",
      "Esophagus", "Stomach",
      "Colon", "Small Intestine",
      "Prostate", "Testis", "Nerve",
      "Pituitary", "Liver",
      "Kidney", "Cervix Uteri",
      "Fallopian Tube", "Bladder",
      "Bone Marrow")

  val <- assertthat::validate_that(
    tissue %in%
      tissue_types_gtex)

  if(!is.logical(val)){
    log4r_info(
      logger,
      paste0(
        "ERROR: tissue type '", tissue,
        "' not among the allowed tissue types in",
        " GTEx: ",
        paste(tissue_types_gtex, collapse=", ")))
    return()
  }

  stopifnot(analysis %in% c("RNASEQ"))

  sample_identifiers <- as.data.frame(
    gtexRNA::sampleMetadata %>%
    dplyr::filter(.data$tissue_type == tissue &
                  (!is.na(.data$in_analysis_freeze) &
                    .data$in_analysis_freeze == analysis)) %>%
    dplyr::filter(!is.na(.data$autolysis_score) &
                    .data$autolysis_score <= max_autolysis_score)
  )

  return(sample_identifiers)

}

# Function that retrieves gene expression data for a given set of
# genes from GTEx. A list of GTEx sample identifiers must be
# provided, for which gene expression estimates are provided (TPM, log2(TPM + 1))
#
# param sample_identifiers data frame with sample identifiers
# param sample_tissue_type tissue type
# param genes character vector with primary gene symbols
# param data_dir Directory with complete TPM matrix from GTEX (gtex_tpm.gct.gz)
#
# #export
#
# get_tpm_data <- function(
#   sample_identifiers = NULL,
#   sample_tissue_type = "Adrenal Gland",
#   data_dir = NULL,
#   genes = NULL){
#
#
#
#   stopifnot(!is.null(genes))
#   stopifnot(dir.exists(data_dir))
#
#   logger <- log4r::logger(
#     threshold = "INFO",
#     appenders = log4r::console_appender(log4r_layout))
#
#   val <- assertthat::validate_that(is.data.frame(sample_identifiers))
#   if(!is.logical(val)){
#     log4r_info(logger,
#                paste0("ERROR: Argument 'sample_identifiers' must be",
#                       "of type 'data.frame' - NOT ",
#                       class(sample_identifiers))
#                )
#   }
#
#   assertable::assert_colnames(
#     sample_identifiers,
#     c('sample_id','tissue_type','tissue_type_detailed',
#       'in_analysis_freeze','autolysis_score'),
#     only_colnames = F)
#
#
#
#   gene_regex <- paste(genes, collapse="|")
#   sample_gene_regex <- paste0("(^Name)|(",gene_regex,")")
#
#   ## random number for file with hits
#   rn <- round(stats::runif(1, min = 0, max = 10000000000))
#
#   ## rough filtering of full GTex RNAseq dataset by regex of genesymbols
#   system(paste0(
#     "gzip -dc ",
#     file.path(data_dir,"gtex_tpm.gz"),
#     "| egrep '",
#     sample_gene_regex,
#     "' > ",
#     file.path(data_dir, paste0("tmp.",rn,"gtex_tpm.txt"))))
#
#   tpm_data <- NULL
#   select_cols <- c('Description','Name')
#   if(!is.null(sample_identifiers)){
#     select_cols <- c(select_cols, sample_identifiers)
#     tpm_data <- data.table::fread(
#       file = file.path(data_dir, paste0("tmp.",rn,"gtex_tpm.txt")),
#
#       ## only choose desired samples
#       select = select_cols,
#       header = T, stringsAsFactors = F, verbose = F) %>%
#       dplyr::rename(symbol = .data$Description,
#                     ensembl_gene_id = .data$Name) %>%
#       ## more exact filtering against gene symbols
#       dplyr::filter(.data$symbol %in% genes)
#   }else{
#     tpm_data <- data.table::fread(
#       file = file.path(data_dir, paste0("tmp.",rn,"gtex_tpm.txt")),
#       header = T, stringsAsFactors = F, verbose = F) %>%
#       dplyr::rename(symbol = .data$Description,
#                     ensembl_gene_id = .data$Name) %>%
#       ## more exact filtering against gene symbols
#       dplyr::filter(.data$symbol %in% genes)
#   }
#
#   tpm_long <- data.frame()
#
#   if(nrow(tpm_data) > 0){
#     tpm_long <- as.data.frame(
#       tpm_data %>%
#         tidyr::pivot_longer(
#           cols = dplyr::starts_with("GTEX"),
#           names_to = "sample_id",
#           values_to = "tpm",
#           values_drop_na = TRUE
#         ) %>%
#         dplyr::mutate(tissue_type = .data$sample_tissue_type,
#                       sample_type = "normal",
#                       gtex_version = "gtex_v8") %>%
#         dplyr::mutate(
#           log2_tpm_plus_1 = log2(.data$tpm + 1))
#     )
#   }
#
#   system(paste0(
#     "rm -f ",
#     file.path(data_dir, paste0("tmp.",rn,"gtex_tpm.txt")))
#   )
#
#
#   return(tpm_long)
#
# }

#
# all_tpm_data <- data.frame()
# for(t in c('Ovary','Pancreas','Prostate',
#            'Skin_Sun_Exposed_Lower_leg',
#            'Thyroid','Stomach',
#            'Esophagus_Mucosa',
#            'Bladder',
#            'Breast_Mammary_Tissue',
#            'Colon_Sigmoid','Colon_Transverse',
#            'Lung','Testis','Uterus')){
#   df <- gtexRNA::get_tpm_data_api(tissue_type = t,
#                             genes = c("SYTL5","RAB27A"))
#
#   df <- df %>% dplyr::left_join(
#     dplyr::select(gtexRNA::sampleMetadata,
#                   tissue_site_detail_id,
#                   num_samples),
#     by = c("tissue_type" = "tissue_site_detail_id")
#   ) %>%
#     dplyr::distinct()
#
#   all_tpm_data <- all_tpm_data %>%
#     dplyr::bind_rows(df)
#
# }


log4r_layout <-
  function(level, ...) {
    paste0(format(Sys.time()), " - ",
           level, " - ", ..., "\n",
           collapse = "")
  }

log4r_info <- function(log4r_logger, msg) {
  log4r::info(log4r_logger, msg)
}

log4r_debug <- function(log4r_logger, msg) {
  log4r::debug(log4r_logger, msg)
}

log4r_warn <- function(log4r_logger, msg) {
  log4r::warn(log4r_logger, msg)
}

log4r_err <- function(log4r_logger, msg) {
  log4r::error(log4r_logger, msg)
}

file_is_writable <- function(path) {
  assertthat::is.string(path) &&
    file.exists(path) &&
    assertthat::is.writeable(path)
}


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
NULL

#' RCurl url.exists
#'
#' @name url_exists
#' @rdname url_exists
#' @keywords internal
#' @importFrom RCurl url.exists
NULL


#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL

utils::globalVariables(c("."))


#' Retrieve sample gene expression data from GTEx
#'
#' Retrieves sample gene expression data from GTEx, using the public API.
#' A tissue type must be provided, and also a character vector with genes of interest.
#' Gene expression estimates are provided for all samples of that tissue type,
#' with two different expression units provided (TPM, log2(TPM + 1))
#'
#' @param tissue_type tissue type
#' @param genes character vector with primary gene symbols
#' @param gtex_version version of GTEx
#'
#' @export
#'
get_tpm_data <- function(
  tissue_type = "Breast_Mammary_Tissue",
  genes = NULL,
  gtex_version = "v8"){

  logger <- log4r::logger(
    threshold = "INFO",
    appenders = log4r::console_appender(log4r_layout))

  if(is.null(genes)){
    log4r_info(logger, msg = paste0(
      "ERROR: argument 'genes' is NULL - ",
      "a character vector of gene symbols is mandatory"))
    return()
  }

  val1 <- assertthat::validate_that(
    tissue_type %in% unique(gtexRNA::sampleMetadata$tissue_site_detail_id)
  )
  if(!is.logical(val1)){
    log4r_info(logger, msg = paste0(
      "ERROR: tissue type ('",
      tissue_type,
      "') not allowed - use any of:\n ",
      paste(sort(unique(gtexRNA::sampleMetadata$tissue_site_detail_id)),
            collapse = "'\n - '"))
    )
    return()
  }

  tcga_study_match <-
    unique(gtexRNA::sampleMetadata[
      gtexRNA::sampleMetadata$tissue_site_detail_id == tissue_type,]$tcga_study_match)

  sample_type <-
    unique(gtexRNA::sampleMetadata[
      gtexRNA::sampleMetadata$tissue_site_detail_id == tissue_type,]$sample_type)

  log4r_info(logger, msg = paste0(
    "Using tissue_type: ", tissue_type,
    " - as provided with argument 'tissue_type'"))

  gene_symbol_df <-
    data.frame('symbol' = genes, stringsAsFactors = F) %>%
    dplyr::inner_join(gtexRNA::geneXref, by = "symbol")

  ensembl_gene_query <- paste(
    gene_symbol_df$ensembl_gene_id,
    collapse = "%2C")


  queryURL <- paste0(
    "https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v8",
    "&gencodeId=", ensembl_gene_query,
    "&tissueSiteDetailId=", tissue_type,
    "&format=json")

  results <- jsonlite::fromJSON(queryURL)
  tpm_long <- data.frame()
  if(length(results$geneExpression$data) == nrow(results$geneExpression)){
    for(i in 1:nrow(results$geneExpression)){
      df <- data.frame('tpm' = results$geneExpression$data[[i]],
                       'symbol' = results$geneExpression[i, ]$geneSymbol,
                       'ensembl_gene_id' =
                         results$geneExpression[i, ]$gencodeId,
                       'tissue_type' = tissue_type,
                       'sample_type' = sample_type,
                       'unit' = results$geneExpression[i, ]$unit,
                       'tcga_study_match' = tcga_study_match,
                       stringsAsFactors = F) %>%
        dplyr::mutate(
          log2_tpm_plus_1 = log2(.data$tpm + 1),
          gtex_version = gtex_version)

      tpm_long <-
        dplyr::bind_rows(tpm_long, df)
    }
  }

  return(tpm_long)

}
