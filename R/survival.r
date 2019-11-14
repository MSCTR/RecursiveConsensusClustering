survival = function(infoFile, attribute, clus, minSample, outputDir){

  library(survminer)
  library(survival)

  data = read.csv(infoFile)
  types = unique(data[, attribute])
  setwd(outputDir)
  pdf("survival.pdf", height = 10, width = 15)
  for(j in 1:length(types)){
    info = subset(data, data[, attribute] == as.character(types[j]))
    freq = as.data.frame(table(info[,clus]))
    freq = subset(freq, freq[,2] > minSample)
    if(nrow(freq) < 2){ next }
    info = subset(info, info[,clus] %in% freq[,1])
    Cluster = paste0(info[,attribute], "_", info[,clus])
    newSub = cbind(info, Cluster)
    fit = survfit(Surv(Delay, Event) ~ RCC, data = newSub)
    ha = ggsurvplot(fit, data = newSub, risk.table = TRUE, pval = T, title = as.character(types[j]))
    print(ha)
  }
dev.off()
}
