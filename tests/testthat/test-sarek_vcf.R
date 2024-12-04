TEST_VCF_FILES_DIR <- fs::path("test_data/vcf")

HAPLOTYPECALLER_VCF_GLOB <- "*.haplotypecaller.filtered_VEP.ann.vcf.gz"
DEEPVARIANT_VCF_GLOB <- "*.deepvariant_VEP.ann.vcf.gz"
STRELKA_VCF_GLOB <- "*.strelka.variants_VEP.ann.vcf.gz"

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

.prepare_finalist_df <- function(finalist_df) {
  finalist_df %>%
    # those are context-dependent (number of converted samples)
    dplyr::select(-c(`Variant occurred`, `Variant occurred [%]`)) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) %>%
    dplyr::arrange(SAMPLE_variant_key)
}

assert_haplotypecaller_finalist_df_equality <- function(run_name) {
  vcf_files <- .get_vcf_files(run_name, HAPLOTYPECALLER_VCF_GLOB)
  finalist_df <- convert_vcf_files_to_finalist(vcf_files, vc_tool = "haplotypecaller", vcf_annotation_tool = "vep") %>%
    .prepare_finalist_df()

  expected_finalist_df <- .load_expected_finalist_df(vcf_files, run_name) %>%
    .prepare_finalist_df()

  df_diff <- diffdf::diffdf(expected_finalist_df, finalist_df)
  df_diff_invalid <- diffdf::diffdf_has_issues(df_diff)

  if (df_diff_invalid) {
    warning("Found differences between expected and output Finalist dataframe")
    print(df_diff)
    expect_false(df_diff_invalid)
  }

  finalist_df <- finalist_df[, colnames(expected_finalist_df)]

  expect_true(all.equal(expected_finalist_df, finalist_df))
}

assert_finalist_csv_equal <- function(run_name, vc_tool, glob) {
  cli::cli_alert_info("{run_name} / {vc_tool} / {glob}")
  withr::with_tempfile("out_file", fileext = ".csv", code = {
    run_convert_vcf_cli(c("--vc_tool", vc_tool, "--output_file", out_file, .get_vcf_files(run_name, glob)))
    expect_snapshot_file(out_file, name = glue("{run_name}_{vc_tool}.csv"))
  })
}

test_that("Conversion to FinalistDX format works (old Sarek 3.1.* VCF/VEP format, haplotypecaller)", {
  assert_haplotypecaller_finalist_df_equality("2021_09_MR-Mikro4")
})

test_that("Conversion to FinalistDX format works (new Sarek 3.4.2 VCF/VEP format, haplotypecaller)", {
  assert_haplotypecaller_finalist_df_equality("2024_05_24_Panel_run_2")
})

test_that("Conversion to FinalistDX format works using CLI (old Sarek 3.1.* VCF/VEP format, all VC tools, CSV output)", {
  run_name <- "2021_09_MR-Mikro4"
  assert_finalist_csv_equal(run_name, "haplotypecaller", HAPLOTYPECALLER_VCF_GLOB)
  assert_finalist_csv_equal(run_name, "deepvariant", DEEPVARIANT_VCF_GLOB)
  assert_finalist_csv_equal(run_name, "strelka", STRELKA_VCF_GLOB)
})

test_that("Conversion to FinalistDX format works using CLI (new Sarek 3.4.2 VCF/VEP format, haplotypecaller, CSV output)", {
  run_name <- "2024_05_24_Panel_run_2"
  assert_finalist_csv_equal(run_name, "haplotypecaller", HAPLOTYPECALLER_VCF_GLOB)
})
