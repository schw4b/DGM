doc:
	rm -rf ./mdm/NAMESPACE mdmwarwick*tar.gz mdmwarwick.Rcheck
	cd mdm;	R -e 'devtools::document()'; cd ..
	echo "importFrom(utils, combn)"   >> ./mdm/NAMESPACE
	echo "importFrom(Rcpp, evalCpp)"      >> ./mdm/NAMESPACE
	echo "importFrom(stats, binom.test, "p.adjust")" >> ./mdm/NAMESPACE
	echo "importFrom(utils, glob2rx)" >> ./mdm/NAMESPACE
	echo "importFrom(stats, cor)" >> ./mdm/NAMESPACE
	echo "importFrom(stats, var)" >> ./mdm/NAMESPACE
	echo "importFrom(stats, quantile)" >> ./mdm/NAMESPACE
	echo "importFrom(data.table, fread)" >> ./mdm/NAMESPACE
	echo "importFrom(reshape2, melt)" >> ./mdm/NAMESPACE
	echo "importFrom(ggplot2, ggplot, aes_string, geom_tile, scale_fill_gradient2, theme, element_blank, element_text, scale_y_reverse, ggtitle, guides)" >> ./mdm/NAMESPACE
	cd mdm; sed -i '1iuseDynLib(mdmwarwick)' NAMESPACE; cd ..

build:
	R CMD build mdm

check:
	R CMD check --as-cran $(file)
