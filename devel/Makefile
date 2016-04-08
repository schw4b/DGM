doc:
	cd mdm; R -e 'devtools::document()'; cd ..
	echo "importFrom("utils", "combn")" >> ./mdm/NAMESPACE

build:
	R CMD build mdm

check:
	R CMD check --as-cran  mdmwarwick_1.0.tar.gz 
