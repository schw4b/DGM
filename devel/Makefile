doc:
	rm -rf ./multdyn/NAMESPACE multdyn*tar.gz multdyn.Rcheck
	cd multdyn;	R -e 'devtools::document()'; cd ..
	echo "importFrom(utils, combn)"   >> ./multdyn/NAMESPACE
	echo "importFrom(Rcpp, evalCpp)"      >> ./multdyn/NAMESPACE
	echo "importFrom(stats, binom.test, "p.adjust")" >> ./multdyn/NAMESPACE
	echo "importFrom(utils, glob2rx)" >> ./multdyn/NAMESPACE
	echo "importFrom(stats, cor)" >> ./multdyn/NAMESPACE
	echo "importFrom(stats, var)" >> ./multdyn/NAMESPACE
	echo "importFrom(stats, quantile)" >> ./multdyn/NAMESPACE
	echo "importFrom(data.table, fread)" >> ./multdyn/NAMESPACE
	echo "importFrom(reshape2, melt)" >> ./multdyn/NAMESPACE
	echo "importFrom(ggplot2, ggplot, aes_string, geom_tile, scale_fill_gradient2, theme, element_blank, element_text, scale_y_reverse, ggtitle, guides)" >> ./multdyn/NAMESPACE
	cd multdyn; sed -i '1iuseDynLib(multdyn)' NAMESPACE; cd ..

build:
	R CMD build multdyn

check:
	R CMD check --as-cran $(file)
