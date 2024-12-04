.get_vcf_convert_parsed_args <- function(args) {
  parser <- argparse::ArgumentParser(description = "Aggregate and convert VCF files to user-friendly format, including FinalistDX Excel format")
  parser$add_argument(
    "--vc_tool",
    default = "haplotypecaller",
    type = "character",
    choices = c("haplotypecaller", "deepvariant", "strelka"),
    help = "Name of the variant calling tool that originally produced the VCF files [default: %(default)s]"
  )
  parser$add_argument(
    "--vcf_annotation_tool",
    default = "vep",
    type = "character",
    choices = c("vep", "common", "snpeff"),
    help = "Name of the variant annotation tool that was used for the original VCF files ('common' is for generic VCFs) [default: %(default)s]"
  )
  parser$add_argument(
    "--no_variant_filtering",
    action="store_true",
    help = "Skip variant filtering by FILTER column in VCF file"
  )
  # parser$add_argument(
  #   "--output_type",
  #   default = "finalist",
  #   type = "character",
  #   choices = c("original", "finalist"),
  #   help = "Output type: 'finalist' is the FinalistDX Excel format, 'original' is just user-friendly parsed VCF"
  # )
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
run_vcf_convert_cli <- function(args = commandArgs(TRUE)) {
  args <- .get_vcf_convert_parsed_args(args)

  cli::cli_alert_info("Input files ({length(args$vcf_files)}):")
  cli::cli_li(items = glue("{{.file {args$vcf_files}}}"))

  agg_df <- convert_vcf_files_to_finalist(args$vcf_files, filter_variants = !args$no_variant_filtering, vc_tool = args$vc_tool, vcf_annotation_tool = args$vcf_annotation_tool)

  cli::cli_alert_info("Writing result to {.file {args$output_file}}")
  write_output_file(agg_df, args$output_file)
}
