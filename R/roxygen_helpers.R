#' @param vcf_annotation_tool A character scalar: which variant annotation tool was used:
#' - `vep`: Variant Effect Predictor
#' - `snpeff`: snpEff
#' - `common`: a generic VCF
#'
#' @name vcf_annotation_tool_param
NULL

#' @param vc_tool A character scalar: variant calling tool used to produce the VCF files.
#' This influences behaviour [calculate_variant_stats()].
#'
#' @name vc_tool_param
NULL

#' @param filter_variants A logical scalar: if `TRUE`, apply [filter_vcf()]
#'
#' @name filter_variants_param
NULL
