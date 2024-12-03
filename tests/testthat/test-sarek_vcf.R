TEST_VCF_FILES_DIR <- fs::path("test_data/vcf")
EXPECTED_FINALIST_FILES_DIR <- fs::path("expected_data/finalist")

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

assert_finalist_df_equality <- function(run_name) {
  vcf_files <- fs::dir_ls(TEST_VCF_FILES_DIR / run_name, glob = "*.haplotypecaller.filtered_VEP.ann.vcf.gz")
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

test_that("Conversion to FinalistDX format works (old Sarek 3.1.* VCF/VEP format)", {
  assert_finalist_df_equality("2021_09_MR-Mikro4")
})

test_that("Conversion to FinalistDX format works (new Sarek 3.4.2 VCF/VEP format)", {
  assert_finalist_df_equality("2024_05_24_Panel_run_2")
})
