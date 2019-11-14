DEG = function(fileName, infoFile, genesFile){
  data = as.data.frame(fread(fileName))
  rownames(data) = data[,1]
  data = data[,-1]
  #data = t(data)

  info = read.csv(infoFile)
  rownames(info) = info[,1]
  info = info[,-1]
  genes = read.csv(genesFile)
  genes = genes[,1]

  mat = subset(data, rownames(data) %in% genes)
  mat = t(mat)
  ord_mat = mat[order(rownames(mat)), ]
  ord_info = info[order(rownames(info)),]

  identical(rownames(ord_mat), rownames(ord_info))
  ord_mat = t(ord_mat)
  Levels = which(grepl("Level_*", colnames(ord_info)))
  previousLevel = 0
  clusArray = rep(1, nrow(ord_info))

  deg_info = ord_info

  level_markers = NULL

  for(li in 1:length(Levels)){
    ord_info = cbind(ord_info, clusArray)
    deg_info = cbind(deg_info, clusArray)
    colNums = Levels[li]

    if(li == 1){
      a = findMarkers(as.matrix(ord_mat), ord_info[,colNums], direction = "up", lfc = 1)
      mat = NULL
      for(i in 1:length(a)){

        clusterInfo = a[[i]]
        clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
        lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
        lfc = apply(lfc_mat, 1, min)
        clusterInfo = cbind(clusterInfo, lfc)
        clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)

        if(nrow(clusterInfo) == 0){
          mat = qpcR:::cbind.na(mat, rep("0", 5))
        } else {
          mat = qpcR:::cbind.na(mat, rownames(clusterInfo))
          level_markers = append(level_markers, rownames(clusterInfo))
        }
        for(rows in 1:nrow(deg_info)){
          if(deg_info[rows,colNums] == i){
            deg_info[rows,colNums] = nrow(clusterInfo)
          }
        }
      }

      mat = mat[,-1]
      colnames(mat) = names(a)
      write.csv(mat, file = paste0(colnames(ord_info)[colNums], ".csv"))
    }else{
      clusters = unique(clusArray)

      for(newClus in 1:length(clusters)){

        newInfo = subset(ord_info, ord_info[,ncol(ord_info)] == clusters[newClus])
        if((length(unique(newInfo[,colNums])) == 1) && (unique(newInfo[,colNums]) == 0)){
          next
        }
        newMat = t(ord_mat)
        newMat = subset(newMat, rownames(newMat) %in% rownames(newInfo))
        newMat = t(newMat)
        a = findMarkers(as.matrix(newMat), newInfo[,colNums], direction = "up", lfc = 1)
        mat = NULL
        for(i in 1:length(a)){
          clusterInfo = a[[i]]
          clusterInfo = subset(clusterInfo, clusterInfo$FDR < 0.01)
          lfc_mat = as.data.frame(clusterInfo[,4:ncol(clusterInfo)])
          lfc = apply(lfc_mat, 1, min)
          clusterInfo = cbind(clusterInfo, lfc)
          clusterInfo = subset(clusterInfo, clusterInfo$lfc > 1)

          if(nrow(clusterInfo) == 0){
            mat = qpcR:::cbind.na(mat, rep("0", 5))
          } else {
            mat = qpcR:::cbind.na(mat, rownames(clusterInfo))
            level_markers = append(level_markers, rownames(clusterInfo))
          }
          write.csv(a[[i]], file = paste0(colnames(ord_info)[colNums], "_k", clusters[newClus],"k_",i, ".csv"))

          for(rows in 1:nrow(deg_info)){
            if((deg_info[rows, ncol(deg_info)] == clusters[newClus]) && (deg_info[rows,colNums] == i)){
              deg_info[rows,colNums] = nrow(clusterInfo)
            }
            if((deg_info[rows,colNums] == 0) && (ord_info[rows, colNums] == 0)){
              deg_info[rows,colNums] = -1
            }
          }
        }
        mat = mat[,-1]
        colnames(mat) = names(a)
        write.csv(mat, file = paste0(colnames(ord_info)[colNums], "_k", clusters[newClus],".csv"))
      }
    }
    previousLevel = previousLevel + 1
    clusArray = paste0(clusArray, ord_info[, colNums])
    ord_info = ord_info[,-(ncol(ord_info))]
    deg_info = deg_info[,-(ncol(deg_info))]
  }

  colnames(deg_info) = gsub("Level_", "marker_", colnames(deg_info))
  Clusters = deg_info$Clusters
  deg_info = deg_info[,-(ncol(deg_info))]
  deg_info = cbind(deg_info, ord_info[, Levels])
  deg_info = cbind(deg_info, Clusters)

  write.csv(deg_info, file = "markerInfo.csv")
  write.csv(level_markers, file = "levelMarkers.csv")
}

tracking = function(fileName, attribute){
  library(gdata)
  library(circlize)
  library(randomcoloR)
  library(assertr)
  library(ComplexHeatmap)
  library(stringr)

  clusterTrackingPlot = function(m, lenCol){
    plot(NULL,xlim=c(-0.1,1),ylim=c(0,1),axes=FALSE,ylab="",xlab = "")
    for(i in 1:nrow(m)){

      rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/(nrow(m) + 1),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(i-1)/(nrow(m) + 1),ncol(m)), col=m[i,],border=NA)
      rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/(nrow(m) + 1),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-i/(nrow(m) + 1),ncol(m)), col="black",border=T, lwd= 10)
    }
    i = lenCol
    rect(xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/(nrow(m) + 1),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(nrow(m)-1)/(nrow(m)),ncol(m)), col=m[i,],border=NA)
    xl = seq(0,1-1/ncol(m),by=1/ncol(m))
    ypos = seq(1,0,by=-1/(nrow(m)+1)) -1/(2*(nrow(m)))
    ypos=ypos[-length(ypos)]
    text(x=-0.1,y=ypos[-length(ypos)],labels=rownames(m), cex=3)
  }

  info = read.csv(fileName)
  rownames(info) = info[,1]
  info = info[,-1]

  colNums = which(grepl("Level*", colnames(info)))
  markerColNums = which(grepl("marker_*", colnames(info)))

  newInfo = info[, c(markerColNums,colNums,ncol(info))]
  colNum = which(attribute == colnames(info))
  newInfo = cbind(newInfo, info[,colNum])

  for(i in 1:(ncol(newInfo) - 1)){
    newInfo[,i] = str_pad(newInfo[,i], 2, pad = "0")
  }
  subMat = newInfo[,-c((ncol(newInfo) - 1))]
  concMat = subMat[,-(1:length(markerColNums))]
  a = col_concat(concMat, sep = "")
  newInfo = cbind(newInfo, a)

  mat = newInfo[order(newInfo[, ncol(newInfo)]),]
  write.csv(mat, file = "sampleOrder.csv")

  mat = mat[,-ncol(mat)]
  df = as.data.frame(mat[, ncol(mat)])
  ha  = HeatmapAnnotation(df = df)
  mat = mat[,-ncol(mat)]
  newMarkers = which(grepl("marker_*", colnames(mat)))
  markerDF = mat[, newMarkers]
  mat = mat[,-(newMarkers)]

  markerDF = data.frame(sapply(markerDF, function(x) as.numeric(as.character(x))))
  markerDF = log2(markerDF + 3)
  markerDF = round(markerDF,0)
  numMark = unmatrix(markerDF, byrow = T)
  uniqueMark = append(c(1,2),as.numeric(sort(unique(numMark))))
  uniqueMark = unique(uniqueMark)
  minMark = uniqueMark[2]
  maxMark = uniqueMark[length(uniqueMark)]
  colfunc <- colorRampPalette(c("white", "red"))
  red_palette = colfunc((maxMark - minMark) + 2)
  if(uniqueMark[1] == 1){red_palette[1] = "#A8A8A8"}

  for(rows in 1:nrow(markerDF)){
    for(cols in 1:ncol(markerDF)){
      cCol = which(as.numeric(markerDF[rows,cols]) == uniqueMark)
      markerDF[rows,cols] = red_palette[cCol]
    }
  }

  types = unique(df[,1])
  dfNew = NULL

  for(i in 1:length(types)){
    attr_id = as.character(types[i])
    a = ifelse(df[,1] == attr_id, "#000000", "#ffffff")
    dfNew = cbind(dfNew, a)
  }
  cols = as.vector(types)
  colnames(dfNew) = cols


  color_mat = mat
  n = length(unique(mat$Clusters))
  nAttr = length(unique(df[,1]))
  palette <- distinctColorPalette(n)

  clus = unique(mat$Clusters)
  k = 1

  color_mat$Clusters = as.character(color_mat$Clusters)

  i = ncol(color_mat)

  for(j in 1:nrow(mat)){
    color_mat[j,i] = palette[as.numeric(mat[j,i])]
  }

  for(i in (ncol(mat) - 1):1){
    previousClus = mat[1,i]
    color_mat[1,i] = color_mat[1, i + 1]
    for(j in 2:nrow(mat)){

      if((mat[j,i] == previousClus) & (mat[j,i] != "00")){
        color_mat[j,i] = color_mat[j-1,i]
      } else {
        color_mat[j,i] = color_mat[j, i + 1]
      }
      previousClus = mat[j,i]
    }
  }

  palette2 <- distinctColorPalette(nAttr)
  df2 = rep(0, nrow(df))
  label = unique(df[,1])
  for(i in 1:length(label)){
    type = as.character(label)[i]
    for(j in 1:nrow(df)){
      if(df[j,1] == type){ df2[j] = palette2[i] }
    }
  }

  color_mat =  color_mat[,-ncol(color_mat)]
  m = t(color_mat)
  m = rbind(t(markerDF),m)
  m = rbind(m, df2)

  rownames(m)[nrow(m)] = "Attribute"
  m[length(colNums)*2,] = paste0(m[length(colNums)*2,], "66")

  mN = rbind(m, t(dfNew))

  pdf(paste0(attribute, "tracking.pdf"), height = 50, width = 100)
  clusterTrackingPlot(mN, length(colNums)*2)
  dev.off()

}

trackingPlot = function(exprsFile, clusterFile, geneFile, attribute, outputDir){
  setwd(outputDir)
  library(scran)
  library(data.table)

  DEG(exprsFile, clusterFile, geneFile)
  tracking("markerInfo.csv", attribute)
}
