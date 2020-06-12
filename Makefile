check: build
	R CMD check --as-cran brm_1.1.tar.gz

build:
	R CMD build .
