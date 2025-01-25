.split_by_char <- function(x, char = ",", as_matrix = TRUE, convert_to = "integer", NA_fill = 0) {
  if (as_matrix) {
    res <- stringr::str_split(x, stringr::fixed(char), simplify = TRUE)
    mode(res) <- convert_to
    res[is.na(res)] <- NA_fill
  } else {
    res <- stringr::str_split(x, stringr::fixed(char)) %>%
      purrr::map(methods::as, convert_to)
  }

  res
}

#' @export
write_output_file <- function(df, output_file) {
  output_file_ext <- fs::path_ext(output_file)

  if (output_file_ext == "csv") {
    readr::write_csv(df, output_file)
  } else if (output_file_ext == "xlsx") {
    writexl::write_xlsx(df, output_file)
  } else if (output_file_ext == "rds") {
    saveRDS(df, output_file)
  } else {
    cli::cli_abort("Cannot write output file: unknown extension '{output_file_ext}' for output file {.file {output_file}}")
  }
}

#' @export
install_cli <- function(out_dir, script = c("convert_vcf.R"), overwrite = TRUE) {
  rlang::arg_match(script)
  out_dir <- fs::path(out_dir)

  script_file <- fs::path(system.file("scripts", package = "convert.vcf", mustWork = TRUE)) / script
  cli::cli_alert_info("Copying {.file {script_file}} to {.file {out_dir}}")
  fs::file_copy(script_file, out_dir, overwrite = overwrite)
  fs::file_chmod(out_dir / script, "+x")
  cli::cli_alert_success("Done")
}
