### gtexRNA - Quick retrieval of tissue-specific expression data from GTEx

#### Overview

This R package provides a simple wrapper function to query [GTEx (v8)](https://gtexportal.org/home/) for gene expression levels in healthy human tissues, using a set of human gene symbols and a particular tissue type as its main arguments. The function utilizes the [GTEx API](https://gtexportal.org/home/api-docs/index.html) to retrieve the data. A brief example:

`devtools::install_github('sigven/gtexRNA')` 

`exp_dist_gtex <- gtexRNA::get_tpm_data(genes = c('KRAS','BRAF'), tissue_type = 'Bladder')`

`exp_dist_gtex` is a data.frame with all gene expression levels pr. sample, provided with TPM and log2(TPM + 1).

To see the available types of tissues that can be queried, explore the metadata file:

`unique(gtexRNA::sampleMetadata$tissue_site_detail_id)`

**PS**: We have not yet found a way to retrieve the actual sample identifiers from the [API method](https://gtexportal.org/home/api-docs/index.html#!/expression/geneExpression), 
so `get_tpm_data()` will currently thus return sample records without any sample identifiers.

### Contact

sigven AT ifi.uio.no

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/sigven/gtexRNA/blob/main/.github/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
