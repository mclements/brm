check: build
	R-devel CMD check --as-cran brm_1.1.tar.gz

build:
	R-devel CMD build .
