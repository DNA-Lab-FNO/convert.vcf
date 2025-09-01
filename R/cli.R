.get_convert_vcf_parsed_args <- function(args) {
  parser <- argparse::ArgumentParser(description = "Aggregate and convert VCF files to user-friendly format, including FinalistDX Excel format")
  parser$add_argument(
    "--vc_tool",
    default = "haplotypecaller",
    type = "character",
    choices = SUPPORTED_VC_TOOLS,
    help = "Name of the variant calling tool that originally produced the VCF files [default: %(default)s]"
  )
  parser$add_argument(
    "--vcf_annotation_tool",
    default = "vep",
    type = "character",
    choices = SUPPORTED_VCF_ANNOTATION_TOOLS,
    help = "Name of the variant annotation tool that was used for the original VCF files ('common' is for generic VCFs) [default: %(default)s]"
  )
  parser$add_argument(
    "--no_variant_filtering",
    action="store_true",
    help = "Skip variant filtering by FILTER column in VCF file"
  )
  parser$add_argument(
    "--output_file",
    required = TRUE,
    type = "character",
    help = "Path to output CSV (.csv), Excel (.xlsx), or R binary object (.rds) file (will be determined based on the extension)"
  )
  parser$add_argument(
    "vcf_files",
    type = "character",
    nargs = "+",
    help = "Paths to VCF files"
  )

  parser$parse_args(args = args)
}

#' @export
run_convert_vcf_cli <- function(args = commandArgs(TRUE)) {
  args <- .get_convert_vcf_parsed_args(args)

  cli::cli_alert_info("Input files ({length(args$vcf_files)}):")
  cli::cli_li(items = glue("{{.file {args$vcf_files}}}"))

  agg_df_finalist <- convert_vcf_files_to_finalist(
    args$vcf_files,
    filter_variants = !args$no_variant_filtering,
    vc_tool = args$vc_tool,
    vcf_annotation_tool = args$vcf_annotation_tool
  )

  write_output_file(agg_df_finalist, args$output_file)

  agg_df_parsed <- convert_vcf_files_to_parsed(
    args$vcf_files,
    filter_variants = !args$no_variant_filtering,
    vc_tool = args$vc_tool,
    vcf_annotation_tool = args$vcf_annotation_tool
  )

  output_vcf_parsed_file <- get_parsed_vcf_file_name(args$output_file)
  write_output_file(agg_df_parsed, output_vcf_parsed_file)
}

get_parsed_vcf_file_name <- function(output_file) {
  glue("{fs::path_ext_remove(output_file)}_original.csv.gz")
}
