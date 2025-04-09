.validate_vcf_files <- function(vcf_files) {
  vcf_files_missing <- vcf_files[!fs::file_exists(vcf_files)]
  if (!rlang::is_empty(vcf_files_missing)) {
    cli::cli_abort("Cannot find the following files: {.file {vcf_files_missing}}")
  }

  duplicated_files <- vcf_files[duplicated(vcf_files)]
  if (length(duplicated_files) > 0) {
    cli::cli_abort("There are duplicated input files: {.file {duplicated_files}}")
  }
}

#' @export
extract_sample_names <- function(files, vc_tool) {
  base_names <- fs::path_file(files)
  sample_names <- stringr::str_extract(base_names, glue("(.+)\\.{vc_tool}(_|\\.).+\\.vcf"), group = 1)

  files_with_na_sample_names <- files[is.na(sample_names)]
  if (length(files_with_na_sample_names) > 0) {
    cli::cli_abort("Unable to extract sample names from files (vc_tool = '{.field {vc_tool}}'): {.file {files_with_na_sample_names}}")
  }

  sample_names
}

.check_duplicated_samples <- function(vcf_files, sample_names) {
  duplicates <- suppressMessages(
    tibble::tibble(vcf_file = vcf_files, sample_name = sample_names) %>%
      janitor::get_dupes()
  )

  if (nrow(duplicates) > 0) {
    cli::cli_alert_warning("Found duplicated samples (the same sample name but different files): {.file {unique(duplicates$vcf_file)}}")
  }
}

#' @title Read VCF file into tidy tibble
#'
#' @description Internally, this function uses a [modified version](https://github.com/gorgitko/vcfR) of the vcfR package
#' that is able to parse annotated VCF (by VEP or snpEff tools) into a tidy tibble.
#' It is necessary to select the proper annotation tool of origin (see the `vcf_annotation_tool` parameter)
#' in order to correctly parse the VCF.
#'
#' Extra steps in this function:
#' 1. Filtering of variants (see [filter_vcf()]) if `filter_variants` is `TRUE`
#' 2. Creation of `variant_key` and `SAMPLE_variant_key` column to uniquely identify a variant and also among multiple samples,
#'    respectively (see [create_variant_keys()])
#' 3. Creation of variants stats (`gt_GT`, `gt_AD` etc., see [calculate_variant_stats()])
#'
#' @param file A character scalar: path to file VCF file (can be also gzipped)
#' @param sample_name A character scalar: sample name to fill in `SAMPLE` column of the output tibble
#'
#' @inheritParams filter_variants_param
#' @inheritParams vc_tool_param
#' @inheritParams vcf_annotation_tool_param
#'
#' @return A tibble, with `meta` attribute storing some additional VCF metadata
#'
#' @export
read_vcf <- function(
    file,
    sample_name,
    filter_variants = TRUE,
    vc_tool = c("haplotypecaller", "deepvariant", "strelka"),
    vcf_annotation_tool = c("vep", "common", "snpeff")) {
  vcf_annotation_tool <- rlang::arg_match(vcf_annotation_tool)

  vcf_tidy <- vcfR::read.vcfR(file, verbose = FALSE) %>%
    vcfR::vcfR2tidy(single_frame = TRUE, verbose = FALSE)

  if (vcf_annotation_tool == "snpeff") {
    vcf_df <- vcfR::separate_ann_snpeff(vcf_tidy$dat, vcf_tidy$meta)
  } else if (vcf_annotation_tool == "vep") {
    vcf_df <- vcfR::separate_ann_vep(vcf_tidy$dat, vcf_tidy$meta)
  } else {
    vcf_df <- vcf_tidy$dat
  }

  attr(vcf_df, "meta") <- vcf_tidy$meta

  vcf_df <- vcf_df %>%
    dplyr::mutate(SAMPLE = !!sample_name, input_file = !!file) %>%
    dplyr::relocate(SAMPLE, input_file)

  if (filter_variants) {
    cli::cli_alert_info("Filtering variants")
    vcf_df <- filter_vcf(vcf_df, vc_tool = vc_tool)
    cli::cli_alert_success("Done")
  }

  vcf_df %>%
    create_variant_keys() %>%
    calculate_variant_stats(vc_tool = vc_tool)
}

#' @export
filter_vcf <- function(vcf_df, vc_tool = c("haplotypecaller", "deepvariant", "strelka")) {
  vc_tool <- rlang::arg_match(vc_tool)

  if (vc_tool == "haplotypecaller") {
    ## currently not filtered
    vcf_df
  } else if (vc_tool == "deepvariant") {
    vcf_df %>%
      dplyr::filter(FILTER != "RefCall", gt_GT != "0/0")
  } else if (vc_tool == "strelka") {
    vcf_df %>%
      dplyr::filter(FILTER == "PASS", stringr::str_detect(gt_GT, stringr::fixed("/")))
  }
}

#' @export
create_variant_keys <- function(vcf_df) {
  vcf_df %>%
    dplyr::mutate(
      variant_key = stringr::str_c(CHROM, POS, REF, ALT, sep = "_"),
      SAMPLE_variant_key = stringr::str_c(SAMPLE, variant_key, sep = "_")
    )
}

#' @export
calculate_variant_stats <- function(vcf_df, vc_tool = c("haplotypecaller", "deepvariant", "strelka")) {
  vc_tool <- rlang::arg_match(vc_tool)

  if (vc_tool %in% c("haplotypecaller", "strelka")) {
    vcf_df %>%
      dplyr::mutate(
        .get_alt_allele_stats_haplotypecaller_strelka(vcf_df$gt_AD, vcf_df$gt_GT)
      )
  } else if (vc_tool == "deepvariant") {
    vcf_df %>%
      dplyr::mutate(
        .get_alt_allele_stats_deepvariant(vcf_df$gt_AD, vcf_df$gt_GT, vcf_df$gt_VAF)
      )
  }
}

.get_alt_allele_stats_haplotypecaller_strelka <- function(gt_AD, gt_GT) {
  assertthat::assert_that(length(gt_AD) == length(gt_GT))

  gt_AD_split <- .split_by_char(gt_AD, as_matrix = FALSE)
  gt_AD_split_m <- .split_by_char(gt_AD, as_matrix = TRUE)
  gt_GT_split <- .split_by_char(gt_GT, char = "/", as_matrix = FALSE)

  tibble::tibble(i = seq_len(length(gt_AD)), gt_GT = gt_GT, gt_GT_split = gt_GT_split, gt_AD_split = gt_AD_split) %>%
    dplyr::group_by(gt_GT) %>%
    dplyr::group_map(function(data, key) {
      gt_GT_split <- data$gt_GT_split[[1]]
      if (0 %in% gt_GT_split) {
        alt_allele_col_indices <- gt_GT_split[2] + 1
        alt_depths <- purrr::map_int(data$gt_AD_split, alt_allele_col_indices)
      } else {
        alt_allele_col_indices <- which(gt_GT_split > 0 & !duplicated(gt_GT_split)) + 1
        alt_depths <- purrr::map_int(data$gt_GT_split, ~ sum(.[alt_allele_col_indices]))
      }

      gt_AD_split_m_subset <- gt_AD_split_m[data$i, , drop = FALSE]
      gt_AD_alt_sums <- rowSums(gt_AD_split_m_subset[, alt_allele_col_indices, drop = FALSE])
      vaf <- gt_AD_alt_sums / rowSums(gt_AD_split_m_subset)

      data %>%
        dplyr::mutate(
          gt_VAF = vaf,
          gt_AD_alt = gt_AD_alt_sums
        )
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(gt_AD_ref = purrr::map_int(gt_AD_split, 1)) %>%
    dplyr::arrange(i)
}

.get_alt_allele_stats_deepvariant <- function(gt_AD, gt_GT, gt_VAF) {
  assertthat::assert_that(length(gt_AD) == length(gt_GT))
  assertthat::assert_that(length(gt_AD) == length(gt_VAF))

  gt_AD_split <- .split_by_char(gt_AD, as_matrix = FALSE)
  gt_split <- .split_by_char(gt_GT, char = "/", as_matrix = FALSE)
  vaf_split <- .split_by_char(gt_VAF, as_matrix = FALSE, convert_to = "double")
  df <- tibble::tibble(
    i = seq_len(length(gt_AD)),
    gt_GT = gt_GT,
    gt_split = gt_split,
    vaf_split = vaf_split,
    gt_AD_split = gt_AD_split
  )

  df %>%
    dplyr::group_by(gt_GT) %>%
    dplyr::group_map(function(data, key) {
      gt_split <- data$gt_split[[1]]
      if (0 %in% gt_split) {
        alt_allele_col_indices <- gt_split[2] + 1
        alt_depth <- purrr::map_int(data$gt_AD_split, alt_allele_col_indices)
        vaf <- purrr::map_dbl(data$vaf_split, alt_allele_col_indices - 1)
      } else {
        alt_allele_col_indices <- which(data$gt_split[[1]] > 0 & !duplicated(data$gt_split[[1]])) + 1
        alt_depth <- purrr::map_int(data$gt_AD_split, ~ sum(.[alt_allele_col_indices]))
        vaf <- purrr::map_dbl(data$vaf_split, ~ sum(.[alt_allele_col_indices - 1]))
      }

      data %>%
        dplyr::mutate(
          gt_VAF = vaf,
          gt_AD_alt = alt_depth
        )
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(gt_AD_ref = purrr::map_int(gt_AD_split, 1)) %>%
    dplyr::arrange(i)
}

#' @title Convert a VCF tibble to the FinalistDX format
#'
#' @param vcf_df A tibble: obtained with [read_vcf()]
#'
#' @return A tibble
#'
#' @export
convert_vcf_df_to_finalist <- function(vcf_df) {
  n_samples <- unique(vcf_df$SAMPLE) %>% length()

  ## count variant occurrences within all samples
  n_ids <- vcf_df %>%
    ## distinct IDs within each sample as there are multiple annotations per a single ID
    dplyr::distinct(SAMPLE, CHROM, POS, ID, REF, ALT) %>%
    dplyr::group_by(CHROM, POS, ID, REF, ALT) %>%
    dplyr::tally() %>%
    dplyr::mutate(
      `Variant occurred` = glue("{n}/{n_samples}") %>% as.character(),
      `Variant occurred [%]` = scales::percent(n / n_samples, suffix = "", accuracy = 0.01)
    )

  vcf_df %>%
    dplyr::left_join(n_ids, by = c("CHROM", "POS", "ID", "REF", "ALT")) %>%
    tidyr::extract(sift, c("sift_class", "sift_score"), regex = "^(.+)\\((.*)\\)$", remove = FALSE, convert = TRUE) %>%
    tidyr::extract(poly_phen, c("poly_phen_class", "poly_phen_score"), regex = "^(.+)\\((.*)\\)$", remove = FALSE, convert = TRUE) %>%
    dplyr::mutate(
      Genotype = dplyr::case_when(
        gt_GT == "0/0" ~ "HOM_ref",
        gt_GT == "1/1" ~ "HOM",
        stringr::str_detect(gt_GT, "0/[1-9]") ~ "HET",
        stringr::str_detect(gt_GT, "[1-9]/[2-9]") ~ "HET_alt",
        TRUE ~ gt_GT
      ),
      Canonical = dplyr::if_else(canonical == "YES", "YES", "NO"),
      `ClinVar Allele ID` = NA_character_,
      `ClinVar disease name` = NA_character_,
      `ClinVar review status` = NA_character_,
      `LRT pred` = NA_real_,
      `CADD phred` = NA_real_,
      `DANN score` = NA_real_,
      `FATHMM pred` = NA_real_,
      `HGVS offset` = NA,
      codonpos = NA,
      refcodon = NA,
      `Depth of reference forward` = NA_integer_,
      `Depth of reference reverse` = NA_integer_,
      `Depth of alternate forward` = NA_integer_,
      `Depth of alternate reverse` = NA_integer_
    ) %>%
    dplyr::select(
      Name = SAMPLE,
      `Variant occurred`,
      `Variant occurred [%]`,
      `Read Depth` = gt_DP,
      `Variant Frequency` = gt_VAF,
      Symbol = symbol,
      Chr = CHROM,
      Genotype,
      Variant_class = variant_class,
      Coordinate = POS,
      Feature = feature,
      `MANE Select` = mane_select,
      `MANE Plus Clinical` = mane_plus_clinical,
      Canonical,
      Exon = exon,
      Intron = intron,
      HGVSc = hgv_sc,
      HGVSp = hgv_sp,
      `VEP dbSNP ID` = ID,
      `dbSNP ID` = ID,
      `ClinVar Allele ID`,
      Impact = impact,
      Consequence = consequence,
      `Clinical significance` = clin_sig,
      `ClinVar disease name`,
      `ClinVar review status`,
      SIFT = sift_class,
      `SIFT score` = sift_score,
      `PolyPhen` = poly_phen_class,
      `PolyPhen score` = poly_phen_score,
      `LRT pred`,
      `CADD phred`,
      `DANN score`,
      `FATHMM pred`,
      `gnomAD AF` = dplyr::any_of(c("gnom_ad_af", "gnom_a_de_af")),
      `gnomAD NFE AF` = dplyr::any_of(c("gnom_ad_nfe_af", "gnom_a_de_nfe_af")),
      `MAX AF` = max_af,
      `MAX AF POPS` = max_af_pops,
      `cDNA position` = c_dna_position,
      `CDS position` = cds_position,
      `Protein position` = protein_position,
      `Amino acids` = amino_acids,
      Codons = codons,
      Distance = distance,
      Strand = strand,
      `HGVS offset`,
      Somatic = somatic,
      Pheno = pheno,
      codonpos,
      refcodon,
      Reference = REF,
      Alternate = ALT,
      `Reference depth` = gt_AD_ref,
      `Alternate depth` = gt_AD_alt,
      `Depth of reference forward`,
      `Depth of reference reverse`,
      `Depth of alternate forward`,
      `Depth of alternate reverse`,
      `gnomAD AF_` = dplyr::any_of(c("gnom_ad_af", "gnom_a_de_af")),
      `gnomAD AFR AF` = dplyr::any_of(c("gnom_ad_afr_af", "gnom_a_de_afr_af")),
      `gnomAD AMR AF` = dplyr::any_of(c("gnom_ad_amr_af", "gnom_a_de_amr_af")),
      `gnomAD ASJ AF` = dplyr::any_of(c("gnom_ad_asj_af", "gnom_a_de_asj_af")),
      `gnomAD EAS AF` = dplyr::any_of(c("gnom_ad_eas_af", "gnom_a_de_eas_af")),
      `gnomAD FIN AF` = dplyr::any_of(c("gnom_ad_fin_af", "gnom_a_de_fin_af")),
      `gnomAD OTH AF` = dplyr::any_of(c("gnom_ad_oth_af", "gnom_a_de_oth_af")),
      `gnomAD SAS AF` = dplyr::any_of(c("gnom_ad_sas_af", "gnom_a_de_sas_af")),
      `Ensembl transcript ID` = feature,
      `Pubmed ID` = pubmed,
      SAMPLE_variant_key
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          `Name`, `Symbol`, `Chr`, `Genotype`, `Variant_class`, `Feature`, `Canonical`, `Impact`, `Consequence`,
          `Clinical significance`, `SIFT`, `PolyPhen`, `Strand`, `Somatic`, `Pheno`, `Ensembl transcript ID`
        ),
        factor
      ),
      dplyr::across(
        c(
          `Variant occurred [%]`, `Variant Frequency`, `SIFT score`, `PolyPhen score`, `gnomAD AF`,
          `gnomAD NFE AF`, `MAX AF`, `gnomAD AF_`, `gnomAD AFR AF`, `gnomAD AMR AF`, `gnomAD ASJ AF`, `gnomAD EAS AF`,
          `gnomAD FIN AF`, `gnomAD OTH AF`, `gnomAD SAS AF`
        ),
        as.numeric
      ),
      dplyr::across(
        c(
          `Read Depth`, `Coordinate`, `Distance`, `Reference depth`, `Alternate depth`
        ),
        as.integer
      )
    ) %>%
    dplyr::mutate(`Variant Frequency` = `Variant Frequency` * 100) %>%
    dplyr::arrange(Name, Chr, Coordinate)
}

#' @export
get_per_sample_variant_stats <- function(vcf_df) {
  vcf_df %>%
    dplyr::group_by(SAMPLE, input_file) %>%
    dplyr::group_map(function(data, key) {
      data %>%
        get_variant_stats() %>%
        dplyr::mutate(SAMPLE = key$SAMPLE, input_file = key$input_file) %>%
        dplyr::relocate(SAMPLE, input_file)
    }) %>%
    dplyr::bind_rows()
}

#' @export
get_variant_stats <- function(vcf_df) {
  vcf_df <- dplyr::distinct(vcf_df, SAMPLE_variant_key, .keep_all = TRUE)
  tibble::tibble(
    n_variants = nrow(vcf_df),
    n_canonical_variants = sum(vcf_df$canonical == "YES", na.rm = TRUE)
  )
}

#' @title Convert VCF files into FinalistDX format and export to Excel
#'
#' @description This function provides the following steps:
#' 1. Extraction of sample names from VCF files (see [extract_sample_names()])
#' 2. Conversion of VCF files to tidy tibble (see [read_vcf()])
#' 3. Conversion to FinalistDX format (see [convert_vcf_df_to_finalist()])
#'
#' @param vcf_files A character scalar/vector: path(s) to VCF files
#' @inheritParams filter_variants_param
#' @inheritParams vc_tool_param
#' @inheritParams vcf_annotation_tool_param
#'
#' @return A tibble
#'
#' @export
convert_vcf_files_to_finalist <- function(
    vcf_files,
    filter_variants = TRUE,
    vc_tool = c("haplotypecaller", "deepvariant", "strelka"),
    vcf_annotation_tool = c("vep", "common", "snpeff")) {
  vc_tool <- rlang::arg_match(vc_tool)

  .validate_vcf_files(vcf_files)
  sample_names <- extract_sample_names(vcf_files, vc_tool = vc_tool)
  .check_duplicated_samples(vcf_files, sample_names)

  cli::cli_alert_info("Starting the conversion")
  agg_df <- purrr::map2_dfr(vcf_files, sample_names, function(vcf_file, sample_name) {
    cli::cli_alert_info("Processing sample {.field {sample_name}} ({.file {vcf_file}})")
    read_vcf(vcf_file, sample_name, filter_variants = filter_variants, vc_tool = vc_tool, vcf_annotation_tool = vcf_annotation_tool)
  })
  cli::cli_alert_success("Done")

  cli::cli_alert_info("Per-sample variant numbers:")
  get_per_sample_variant_stats(agg_df) %>%
    print(n = Inf)

  cli::cli_alert_info("Converting to FinalistDX format")
  agg_df <- convert_vcf_df_to_finalist(agg_df)
  cli::cli_alert_success("Done")

  agg_df
}
