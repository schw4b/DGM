clean:
	rm -rf ./DGM/NAMESPACE DGM*tar.gz DGM.Rcheck
	echo "# Generated by roxygen2: do not edit by hand" > ./DGM/NAMESPACE



namespace:
	echo "importFrom(utils, combn)"   >> ./DGM/NAMESPACE
	echo "importFrom(stats, binom.test, "p.adjust")" >> ./DGM/NAMESPACE
	echo "importFrom(utils, glob2rx)" >> ./DGM/NAMESPACE
	echo "importFrom(stats, cor)" >> ./DGM/NAMESPACE
	echo "importFrom(stats, var)" >> ./DGM/NAMESPACE
	echo "importFrom(stats, quantile)" >> ./DGM/NAMESPACE
	echo "importFrom(data.table, fread)" >> ./DGM/NAMESPACE
	echo "importFrom(reshape2, melt)" >> ./DGM/NAMESPACE
	echo "importFrom(ggplot2, ggplot, aes_string, geom_tile, scale_fill_gradient2, theme, element_blank, element_text, scale_y_reverse, scale_x_continuous, ggtitle, guides, waiver, guide_colorbar)" >> ./DGM/NAMESPACE
	echo "importFrom(stats, pnorm, t.test)" >> ./DGM/NAMESPACE
	echo "importFrom(coin, oneway_test, approximate, statistic, pvalue)" >> ./DGM/NAMESPACE
	sed -i '1s/^/importFrom(Rcpp, evalCpp)\n/' ./DGM/NAMESPACE	
	sed -i '1s/^/exportPattern("^[[:alpha:]]+")\n/' ./DGM/NAMESPACE
	sed -i '1s/^/useDynLib(DGM, .registration = TRUE)\n/' ./DGM/NAMESPACE

build:
	R CMD build DGM

check:
	R CMD check --as-cran $(file)
