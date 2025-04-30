test_that(".parse_genotype parses basic genotypes correctly", {
  expect_equal(.parse_genotype("0/0"), "HOM_ref")
  expect_equal(.parse_genotype("0|0"), "HOM_ref")
  expect_equal(.parse_genotype("1/1"), "HOM")
  expect_equal(.parse_genotype("0/1"), "HET")
  expect_equal(.parse_genotype("0|1"), "HET")
  expect_equal(.parse_genotype("1/2"), "HET_alt")
  expect_equal(.parse_genotype("2|3"), "HET_alt")
})

test_that(".parse_genotype parses multi-allelic genotypes correctly", {
  expect_equal(.parse_genotype("0/1/2"), "MULTI_HET")
  expect_equal(.parse_genotype("0|2|3"), "MULTI_HET")
  expect_equal(.parse_genotype("1/2/3"), "MULTI_ALT")
  expect_equal(.parse_genotype("2|3|4"), "MULTI_ALT")
})

test_that(".parse_genotype handles missing and unknown cases", {
  expect_equal(.parse_genotype("."), NA_character_)
  expect_equal(.parse_genotype(NA_character_), NA_character_)
  expect_equal(.parse_genotype("0/5"), "HET")
  expect_equal(.parse_genotype("5/5/5"), "MULTI_ALT")
})
