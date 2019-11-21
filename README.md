# RecursiveConsensusClustering
RCC, an unsurpersived clustering algorithm for novel sub type discovery from bulk and single-cell transcriptomic data sets

For sample input and ouput files please download data from: https://www.msctr.org/2019/05/30/recursive-consensus-clustering/

To download the latest version of the package please use the following commands:

install.packages(c("devtools", "clv", "fields", "matrixStats", "data.table", "cluster", "clue", "circlize", "gdata"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("scran", "ComplexHeatmap"))

devtools::install_github("hoxo-m/pforeach")

devtools::install_github("MSCTR/RecursiveConsensusClustering")
