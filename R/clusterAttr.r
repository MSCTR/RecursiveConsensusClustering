clusterAttr = function(infoFile, attribute, minSample, category, algorithm, outputDir){
  data = read.csv(infoFile)
  types = unique(data[, attribute])
  pvalMat = matrix(NA, nrow = length(types), ncol = 2)
  for(subtypes in 1:length(types)){
    pvalMat[subtypes, 1] = as.character(types[subtypes])
    info = subset(data, data[, attribute] == as.character(types[subtypes]))
    freq = as.data.frame(table(info[, algorithm]))
    freq = subset(freq, freq[,2] > 10)
    if(nrow(freq) < 2){ next }
    info = subset(info, info[, algorithm] %in% freq[,1])
    subClus = unique(info[, algorithm])
    if(length(unique(subClus)) < 2){ next }

    subClass = unique(info[, category])
    if(length(unique(subClass)) < 2){ next }
    mat = matrix(0, nrow = length(subClass), ncol = length(subClus))
    for(i in 1:nrow(mat)){
      cat = subClass[i]
      for(j in 1:ncol(mat)){
        clus = subClus[j]
        n = subset(info, (info[, algorithm] == clus) & (info[, category] == as.character(cat)))
        mat[i,j] = nrow(n)
      }
    }
    colnames(mat) = paste0("K_", subClus)
    rownames(mat) = subClass
    if(ncol(as.data.frame(mat[rowSums(mat) > 0, ])) < 2){ next }
    ft = fisher.test(mat, simulate.p.value=T, alternative = "two.sided")
    pvalMat[subtypes, 1] = as.character(types[subtypes])
    pvalMat[subtypes, 2] = ft$p.value
  }
  setwd(outputDir)
  colnames(pvalMat) = c(attribute, "pval")
  write.csv(pvalMat, file = "FE.csv", row.names = F)
}
