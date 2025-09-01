check:
	Rscript -e 'rcmdcheck::rcmdcheck(error_on = "error", args = "--no-manual", build_args = "--no-resave-data")'
