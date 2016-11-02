doc:
	rm -rf ./mdm/NAMESPACE mdmwarwick*tar.gz mdmwarwick.Rcheck
	cd mdm;	R -e 'devtools::document()'; cd ..
	echo "importFrom(utils, combn, read.table)"   >> ./mdm/NAMESPACE
	echo "importFrom(Rcpp, evalCpp)"      >> ./mdm/NAMESPACE
	echo "importFrom(igraph, graph.adjacency, plot.igraph)" >> ./mdm/NAMESPACE
	echo "importFrom(graphics, image, mtext, par)" >> ./mdm/NAMESPACE
	echo "importFrom(grDevices, heat.colors)" >> ./mdm/NAMESPACE
	echo "importFrom(fields, image.plot)" >> ./mdm/NAMESPACE
	echo "importFrom("graphics", "grid")" >> ./mdm/NAMESPACE
	echo "importFrom("stats", "binom.test", "p.adjust")" >> ./mdm/NAMESPACE
	echo "importFrom("utils", "glob2rx")" >> ./mdm/NAMESPACE
	echo "importFrom("stats", "cor")" >> ./mdm/NAMESPACE
	echo "importFrom("stats", "var")" >> ./mdm/NAMESPACE

	#echo "exportPattern("^[[:alpha:]]+")" >> ./mdm/NAMESPACE
	cd mdm; sed -i '1iuseDynLib(mdmwarwick)' NAMESPACE; cd ..

build:
	R CMD build mdm

check:
	R CMD check --as-cran  mdmwarwick_1.3.tar.gz 
