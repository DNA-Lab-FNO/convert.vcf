TEST_VCF_FILES_DIR <- fs::path("test_data/vcf")

HAPLOTYPECALLER_VCF_GLOB <- "*.haplotypecaller.filtered_VEP.ann.vcf.gz"
DEEPVARIANT_VCF_GLOB <- "*.deepvariant_VEP.ann.vcf.gz"
STRELKA_VCF_GLOB <- "*.strelka.variants_VEP.ann.vcf.gz"
MUTECT2_VCF_GLOB <- "*.mutect2.filtered_VEP.ann.vcf.gz"

EXPECTED_FINALIST_FILES_DIR <- fs::path("expected_data/finalist")

.get_vcf_files <- function(run_name, glob) {
  fs::dir_ls(TEST_VCF_FILES_DIR / run_name, glob = glob)
}

.load_expected_finalist_df <- function(vcf_files, run_name) {
  purrr::map_dfr(vcf_files, function(vcf_file) {
    base_filename <- fs::path_file(vcf_file)
    expected_filename <- glue::glue("{EXPECTED_FINALIST_FILES_DIR}/{run_name}/{base_filename}.rds")
    readRDS(expected_filename)
  })
}

.update_expected_finalist_df <- function(vcf_files, run_name, expected_finalist_df) {
  purrr::map_dfr(vcf_files, function(vcf_file) {
    base_filename <- fs::path_file(vcf_file)
    expected_filename <- glue::glue("{EXPECTED_FINALIST_FILES_DIR}/{run_name}/{base_filename}.rds")
    fs::dir_create(fs::path_dir(expected_filename), recurse = TRUE)
    expected_finalist_df_sample <- expected_finalist_df %>%
      dplyr::filter(stringr::str_detect(fs::path_file(vcf_file), as.character(Name)))
    message(glue::glue("Updated {expected_filename} (sample {unique(expected_finalist_df_sample$Name)})"))
    saveRDS(expected_finalist_df_sample, expected_filename)
  })

  invisible(NULL)
}

.prepare_finalist_df <- function(finalist_df) {
  finalist_df %>%
    # those are context-dependent (number of converted samples)
    dplyr::select(-c(`Variant occurred`, `Variant occurred [%]`)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) %>%
    dplyr::arrange(SAMPLE_variant_key)
}

assert_finalist_df_equality <- function(run_name, vc_tool = c("haplotypecaller", "mutect2")) {
  vc_tool <- rlang::arg_match(vc_tool)

  if (vc_tool == "haplotypecaller") {
    glob <- HAPLOTYPECALLER_VCF_GLOB
  } else if (vc_tool == "mutect2") {
    glob <- MUTECT2_VCF_GLOB
  }

  vcf_files <- .get_vcf_files(run_name, glob)
  finalist_df_orig <- convert_vcf_files_to_finalist(vcf_files, vc_tool = vc_tool, vcf_annotation_tool = "vep")
  finalist_df <- .prepare_finalist_df(finalist_df_orig)
  expected_finalist_df <- .load_expected_finalist_df(vcf_files, run_name) %>%
    .prepare_finalist_df()

  df_diff <- diffdf::diffdf(expected_finalist_df, finalist_df)
  df_diff_invalid <- diffdf::diffdf_has_issues(df_diff)

  if (df_diff_invalid) {
    warning("Found differences between expected and output Finalist dataframe")
    print(df_diff)
    # .update_expected_finalist_df(vcf_files, run_name, finalist_df_orig)
    expect_false(df_diff_invalid)
  }

  finalist_df <- finalist_df[, colnames(expected_finalist_df)]

  expect_true(all.equal(expected_finalist_df, finalist_df))
}

assert_csv_equal <- function(run_name, vc_tool, glob) {
  cli::cli_alert_info("{run_name} / {vc_tool} / {glob}")
  withr::with_tempfile("out_file", fileext = ".csv", code = {
    run_convert_vcf_cli(c("--vc_tool", vc_tool, "--output_file", out_file, .get_vcf_files(run_name, glob)))
    expect_snapshot_file(out_file, name = glue("{run_name}_{vc_tool}.csv"))

    parsed_vcf_csv_gz_file <- get_parsed_vcf_file_name(out_file)
    parsed_vcf_csv_file <- fs::path_ext_remove(parsed_vcf_csv_gz_file)
    .decompress_gzip(parsed_vcf_csv_gz_file, parsed_vcf_csv_file)
    expect_snapshot_file(parsed_vcf_csv_file, name = glue("{run_name}_{vc_tool}_original.csv"))
  })
}

.decompress_gzip <- function(infile, outfile) {
  con_in  <- gzfile(infile, "rb")
  con_out <- file(outfile, "wb")

  repeat {
    buf <- readBin(con_in, what = raw(), n = 65536)
    if (length(buf) == 0) break
    writeBin(buf, con_out)
  }

  close(con_in)
  close(con_out)

  invisible(NULL)
}

test_that("Conversion to FinalistDX format works (old Sarek 3.1.* VCF/VEP format, haplotypecaller)", {
  assert_finalist_df_equality("2021_09_MR-Mikro4", vc_tool = "haplotypecaller")
})

test_that("Conversion to FinalistDX format works (new Sarek 3.4.2 VCF/VEP format, haplotypecaller)", {
  assert_finalist_df_equality("2024_05_24_Panel_run_2", vc_tool = "haplotypecaller")
})

test_that("Conversion to FinalistDX format works using CLI (old Sarek 3.1.* VCF/VEP format, all VC tools except mutect2, CSV output)", {
  run_name <- "2021_09_MR-Mikro4"
  assert_csv_equal(run_name, "haplotypecaller", HAPLOTYPECALLER_VCF_GLOB)
  assert_csv_equal(run_name, "deepvariant", DEEPVARIANT_VCF_GLOB)
  assert_csv_equal(run_name, "strelka", STRELKA_VCF_GLOB)
})

test_that("Conversion to FinalistDX format works using CLI (new Sarek 3.4.2 VCF/VEP format, haplotypecaller, CSV output)", {
  run_name <- "2024_05_24_Panel_run_2"
  assert_csv_equal(run_name, "haplotypecaller", HAPLOTYPECALLER_VCF_GLOB)
})

test_that("Conversion to FinalistDX format works (new Sarek 3.4.2 VCF/VEP format, mutect2)", {
  assert_finalist_df_equality("2025_04_03_Panel_run_23_BRCA_255_tkane", vc_tool = "mutect2")
})

test_that("Conversion to FinalistDX format works using CLI (old Sarek 3.1.* VCF/VEP format, mutect2, CSV output)", {
  run_name <- "2025_04_03_Panel_run_23_BRCA_255_tkane"
  assert_csv_equal(run_name, "mutect2", MUTECT2_VCF_GLOB)
})
