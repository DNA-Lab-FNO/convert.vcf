.split_by_char <- function(x, char = ",", as_matrix = TRUE, convert_to = "integer", NA_fill = 0, regex = FALSE) {
  if (as_matrix) {
    res <- stringr::str_split(x, pattern = if (regex) char else stringr::fixed(char), simplify = TRUE)
    mode(res) <- convert_to
    res[is.na(res)] <- NA_fill
  } else {
    res <- stringr::str_split(x, pattern = if (regex) char else stringr::fixed(char)) %>%
      purrr::map(methods::as, convert_to)
  }

  res
}

#' @export
write_output_file <- function(df, output_file, xlsx_max_rows_per_file = 999999) {
  output_file_ext <- fs::path_ext(output_file)

  if (output_file_ext == "csv" || stringr::str_ends(output_file, stringr::fixed(".csv.gz"))) {
    cli::cli_alert_info("Writing to {.file {output_file}}")
    readr::write_csv(df, output_file, quote = "all")
  } else if (output_file_ext == "xlsx") {
    .write_xlsx_file(df, output_file, xlsx_max_rows_per_file)
  } else if (output_file_ext == "rds") {
    cli::cli_alert_info("Writing to {.file {output_file}}")
    saveRDS(df, output_file)
  } else {
    cli::cli_abort("Cannot write output file: unknown extension '{output_file_ext}' for output file {.file {output_file}}")
  }
}

.write_xlsx_file <- function(df, output_file, xlsx_max_rows_per_file) {
  max_rows <- 999999
  n_rows <- nrow(df)
  output_base_name <- fs::path_ext_remove(output_file)

  if (n_rows <= xlsx_max_rows_per_file) {
    cli::cli_alert_info("Writing to {.file {output_file}}")
    writexl::write_xlsx(df, output_file)
  } else {
    n_chunks <- ceiling(n_rows / xlsx_max_rows_per_file)
    cli::cli_alert_warning("Dataframe contains more than {xlsx_max_rows_per_file} rows - splitting output file to {n_chunks} chunks")
    for (i in seq_len(n_chunks)) {
      start <- (i - 1) * xlsx_max_rows_per_file + 1
      end <- min(i * xlsx_max_rows_per_file, n_rows)

      output_chunk_file <- glue("{output_base_name}_00{i}.xlsx")

      cli::cli_alert_info("Writing chunk {i}/{n_chunks} to {.file {output_chunk_file}}")
      writexl::write_xlsx(df[start:end, , drop = FALSE], output_chunk_file)
    }
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
