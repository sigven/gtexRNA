#' Sample properties of GTEx RNAseq samples (v8)
#'
#' @format A data.frame with seven columns, and 17,832 rows:
#' \itemize{
#'   \item \emph{sample_id} - GTEx sample identifier
#'   \item \emph{tissue_type} - Tissue type annotation
#'   \item \emph{tissue_type_detailed} - Detailed annotation of tissue type
#'   \item \emph{in_analysis_freeze} - Type of analysis (RNASEQ)
#'   \item \emph{autolysis_score} - Autolysis score (1 to 3)
#'   \item \emph{tissue_site_detail_id} - Slight modification of 'tissue_type_detailed', used for querying API
#'   \item \emph{num_tissue_samples} - Number of samples for a given 'tissue_site_detail_id'
#' }
#'
"sampleMetadata"

#' Cross-reference between gene symbols and Ensembl gene identifiers
#'
#' @format A data.frame with two columns, and 56,200 rows:
#' \itemize{
#'   \item \emph{symbol} - Primary gene symbol
#'   \item \emph{ensembl_gene_id} - Ensembl gene identifier
#' }
#'
"geneXref"
