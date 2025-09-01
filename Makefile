check:
	Rscript -e 'rcmdcheck::rcmdcheck(error_on = "error", args = c("--no-manual", "--no-tests"), build_args = "--no-resave-data")'

test:
	Rscript -e 'devtools::test()'
