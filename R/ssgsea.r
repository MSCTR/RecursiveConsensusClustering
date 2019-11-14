ssgseaFC = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)

    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos

      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

      step_cdf_diff = step_cdf_pos - step_cdf_neg

      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes

      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })

  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))

  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

ssgsea = function(exprsFile, infoFile, geneSet, attribute, clus, outputDir){
  library(matrixStats)
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)

  data = as.data.frame(fread(exprsFile))
  rownames(data) = data[,1]
  data = data[,-1]
  data = as.matrix(data)

  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]

  data = t(data)
  data = data[order(rownames(data)), ]
  info = info[order(rownames(info)),]

  gene_set = read.csv(geneSet)
  gene_sets = as.list(as.data.frame(gene_set))

  res = ssgseaFC(t(data), gene_sets, scale = TRUE, norm = FALSE)
  setwd(outputDir)
  types = unique(info[, attribute])
  pdf("ssgsea.pdf", height = 20, width = 15)
  for(j in 1:length(types)){
    subInfo = subset(info, info[,attribute] == as.character(types[j]))
    mat = t(res)
    mat = subset(mat, rownames(mat) %in% rownames(subInfo))
    mat = cbind(subInfo[, clus], mat)
    colnames(mat)[1] = "Clusters"
    mat = mat[order(mat[,1]),]
    df = as.data.frame(mat[,1])
    mat = mat[,-1]
    mat = t(mat)
    zscore_mat = (mat- rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
    df[,1] = as.character(df[,1])
    ha = HeatmapAnnotation(df = df)
    ht = Heatmap(zscore_mat, cluster_columns = F, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")), top_annotation = ha, column_title = as.character(types[j]))
    draw(ht)
  }
  dev.off()
}
