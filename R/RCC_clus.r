ConsensusClusterPlus <- function( d=NULL,
                                  maxK = 3,
                                  reps=10,
                                  pItem=0.8,
                                  pFeature=1,
                                  clusterAlg="hc",
                                  title="untitled_consensus_cluster",
                                  innerLinkage="average",
                                  finalLinkage="average",
                                  distance="pearson",
                                  ml=NULL,
                                  tmyPal=NULL,
                                  seed=NULL,
                                  plot=NULL,
                                  writeTable=FALSE,
                                  weightsItem=NULL,
                                  weightsFeature=NULL,
                                  verbose=F,
				  corUse="everything" ) {
  ##description: runs consensus subsamples
  #if(is.null(seed)==TRUE){
  #  seed=timeSeed = as.numeric(Sys.time())
  #}
  #distance=ifelse( inherits(d,"dist"), attr( d, "method" ), "pearson" )


  if(is.null(ml)==TRUE){

    if ( ! class( d ) %in% c( "dist", "matrix", "ExpressionSet" ) ) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }

    if ( inherits( d, "dist" ) ) {
      ## if d is a distance matrix, fix a few things so that they don't cause problems with the analysis
      ##  Note, assumption is that if d is a distance matrix, the user doesn't want to sample over the row features
      if ( is.null( attr( d, "method" ) ) ) {
        attr( d, "method" ) <- distance <- "unknown - user-specified"
      }
      if ( is.null( distance ) || ( distance != attr( d, "method" ) ) ) {
        distance <- attr( d, "method" )
      }

      if ( ( ! is.null( pFeature ) ) && ( pFeature < 1 ) ) {
        message( "Cannot use the pFeatures parameter when specifying a distance matrix as the data object\n" )
        pFeature <- 1
      }
      if ( ! is.null( weightsFeature ) ) {
        message( "Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n" )
        weightsFeature <- NULL
      }
      if ( clusterAlg == "km" ) {
        message( "Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix")
        ##d <- as.matrix( d )  #this is now done w/in ccRun
      }
    } else {
      if ( is.null( distance ) ) {
        ## we should never get here, but just in case
        distance <- "pearson"
      }
    }

    if ( ( clusterAlg == "km" ) && inherits( distance, "character" ) && ( distance != "euclidean" ) ) {
      message( "Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the 'kmdist' option, or use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
      distance <- 'euclidean'
    }


    if ( inherits( d,"ExpressionSet" ) ) {
      d <- exprs(d)
    }

    ml <- ccRun( d=d,
                 maxK=maxK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
		 corUse=corUse)
  }
  res=list();

  ##make results directory
  if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
    dir.create(paste(title,sep=""))
  }

  ##write log file
  log <- matrix( ncol=2,
                 byrow=T,
                 c("title",title,
                   "maxK",maxK,
                   "input matrix rows",ifelse ( inherits( d, "matrix" ), nrow(d), "dist-mat" ),
                   "input matrix columns",ifelse ( inherits( d, "matrix" ), ncol(d), ncol( as.matrix(d) ) ),
                   "number of bootstraps",reps,
                   "item subsampling proportion",pItem,
                   "feature subsampling proportion",ifelse( is.null(pFeature), 1, pFeature ),
                   "cluster algorithm",clusterAlg,
                   "inner linkage type",innerLinkage,
                   "final linkage type",finalLinkage,
                   "correlation method",distance,
                   "plot",if(is.null(plot)) NA else plot,
                   "seed",if(is.null(seed)) NA else seed))
  colnames(log) = c("argument","value")
  if(writeTable){
    write.csv(file=paste(title,"/",title,".log.csv",sep=""), log,row.names=F)
  }
  if(is.null(plot)){
    ##nothing
  }else if(plot=="pngBMP"){
    bitmap(paste(title,"/","consensus%03d.png",sep=""))
  }else if(plot=="png"){
    png(paste(title,"/","consensus%03d.png",sep=""))

  }else if (plot=="pdf"){
    pdf(onefile=TRUE, paste(title,"/","consensus.pdf",sep=""))
  }else if (plot=="ps"){
    postscript(onefile=TRUE, paste(title,"/","consensus.ps",sep=""))
  }

  colorList=list()
  colorM = rbind() #matrix of colors.

  #18 colors for marking different clusters
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea", #magenta
               "#2ef4ca", #aqua
               "#f4cced", #pink,
               "#f4cc03", #lightorange
               "#05188a", #navy,
               "#e5a25a", #light brown
               "#06f106", #bright green
               "#85848f", #med gray
               "#000000", #black
               "#076f25", #dark green
               "#93cd7f",#lime green
               "#4d0776", #dark purple
               "#ffffff" #white
               )

  ##plot scale
  colBreaks=NA
  if(is.null(tmyPal)==TRUE){
    colBreaks=10
    tmyPal = myPal(colBreaks)
  }else{
    colBreaks=length(tmyPal)
  }
  sc = cbind(seq(0,1,by=1/( colBreaks) )); rownames(sc) = sc[,1]
  sc = cbind(sc,sc)
  heatmap(sc, Colv=NA, Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=rownames(sc),labCol=F,main="consensus matrix legend")

  for (tk in 2:maxK){
    if(verbose){
      message(paste("consensus ",tk))
    }
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    message("clustered")
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(class(d)=="dist"){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm

    colorList = setClusterColors(res[[tk-1]][[3]],ct,thisPal,colorList)
    pc = c
    pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.



    if(!is.null(plot) && plot=="pngBMP"){
         pc = pc[,hc$order ] #mod for no tree
        pc = rbind(pc,0)
	#no dendrogram if pngBMP
	oc = colorList[[1]][hc$order] #mod for no tree
	heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", col = tmyPal, na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=",
                tk, sep = ""), ColSideCol = oc)
    }else{
        pc = rbind(pc,0)
    #former with tree:
    	heatmap(pc, Colv=as.dendrogram(hc), Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=F,labCol=F,mar=c(5,5),main=paste("consensus matrix k=",tk,sep="") , ColSideCol=colorList[[1]])
    }

    legend("topright",legend=unique(ct),fill=unique(colorList[[1]]),horiz=FALSE )

    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]],clrs=colorList)
    colorM = rbind(colorM,colorList[[1]])
  }
  CDF(ml)
  clusterTrackingPlot(colorM[,res[[length(res)]]$consensusTree$order])
  if(is.null(plot)==FALSE){
    dev.off();
  }
  res[[1]] = colorM
  if(writeTable){
    for(i in 2:length(res)){
      write.csv(file=paste(title,"/",title,".k=",i,".consensusMatrix.csv",sep=""), res[[i]]$consensusMatrix)
      write.table(file=paste(title,"/",title,".k=",i,".consensusClass.csv",sep=""), res[[i]]$consensusClass,col.names = F,sep=",")
    }
  }
  return(res)
}


calcICL = function(res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE){
  #calculates and plots cluster consensus and item consensus
  cc=rbind()
  cci = rbind()
  sumRes=list()
  colorsArr=c()

  #make results directory
  if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
	dir.create(paste(title,sep=""))
  }
  if(is.null(plot)){
    #to screen
  }else if(plot=="pdf"){
    pdf(onefile=TRUE, paste(title,"/","icl.pdf",sep=""))
  }else if(plot=="ps"){
    postscript(onefile=TRUE, paste(title,"/","icl.ps",sep=""))
  }else if (plot=="png"){
    png(paste(title,"/","icl%03d.png",sep=""))
  }else if (plot=="pngBMP"){
    bitmap(paste(title,"/","icl%03d.png",sep=""))
  }

  par(mfrow=c(3,1),mar=c(4,3,2,0))

  for (k in 2:length(res)){ #each k
    eiCols = c();
    o = res[[k]]
    m = o$consensusMatrix
    m = triangle(m,mode=2)
    for (ci in sort(unique(o$consensusClass))){ #each cluster in k
	items = which(o$consensusClass==ci)
	nk = length(items)
	mk = sum( m[items,items], na.rm=T)/((nk*(nk-1))/2)
	cc=rbind(cc,c(k,ci,mk)) #cluster-consensus

      for (ei in rev(res[[2]]$consensusTree$order) ){
		denom = if (ei %in% items) { nk - 1} else { nk }
        	mei = sum( c(m[ei,items],m[items,ei]), na.rm=T)/denom  # mean item consensus to a cluster.
		cci = rbind(cci,c(k,ci,ei,mei)) #cluster, cluster index, item index, item-consensus
      }
      eiCols = c(eiCols, rep(ci,length(o$consensusClass)) )
    }

	  cck = cci[which(cci[,1]==k),] #only plot the new k data.

	  #group by item, order by cluster i
	  w=lapply(split(cck,cck[,3]), function(x) { y=matrix(unlist(x),ncol=4); y[order(y[,2]),4] })
	  q = matrix(as.numeric(unlist(w)),ncol=length(w),byrow=F)
	  q = q[,res[[2]]$consensusTree$order] #order by leave order of k=2
 	  #q is a matrix of k rows and sample columns, values are item consensus of sample to the cluster.

	  thisColors = unique(cbind(res[[k]]$consensusClass,res[[k]]$clrs[[1]]))
	  thisColors=thisColors[order(as.numeric(thisColors[,1])),2]
	  colorsArr=c(colorsArr,thisColors)
	  sumRes[[k]] = rankedBarPlot(q,thisColors,cc=res[[k]]$consensusClass[res[[2]]$consensusTree$order],paste("k=",k,sep="") )
  }

  ys=cs=lab=c()
  lastk=cc[1,1]
  for(i in 1:length(colorsArr)){
    if(lastk != cc[i,1]){
      ys=c(ys,0,0)
      cs=c(cs,NA,NA)
      lastk=cc[i,1]
      lab=c(lab,NA,NA)
    }
    ys=c(ys,cc[i,3])
    cs=c(cs,colorsArr[i])
    lab=c(lab,cc[i,1])
  }
  names(ys) = lab
  par(mfrow=c(3,1),mar=c(4,3,2,0))
  barplot(ys,col=cs,border=cs,main="cluster-consensus",ylim=c(0,1),las=1)
  if(is.null(plot)==FALSE){
	  dev.off()
  }
  colnames(cc) = c("k","cluster","clusterConsensus")
  colnames(cci) = c("k","cluster","item","itemConsensus")
  cci[,"item"] = names(res[[2]]$consensusClass)[ cci[,"item"] ]
  #type cci
  cci = data.frame( k=as.numeric(cci[,"k"]), cluster=as.numeric(cci[,"cluster"]), item=cci[,"item"], itemConsensus=as.numeric(cci[,"itemConsensus"]))

  #write to file.
  if(writeTable){
	write.csv(file=paste(title,"/",title,".summary.cluster.consensus.csv",sep=""),row.names=F, cc)
	write.csv(file=paste(title,"/",title,".summary.item.consensus.csv",sep=""), row.names=F, cc)
  }
  return(list(clusterConsensus=cc,itemConsensus=cci))
}


ccRun <- function( d=d,
                   maxK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
		   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);

  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL

  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                      "pearson", "spearman" )

  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
         ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
         ( is.null( pFeature ) ||
           ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) stop("unsupported distance.")

	if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
	}else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }


  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    seedN = (runif(1, min = 1, max = 100000000))/10000
	set.seed(seedN)
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )

    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as well), or
      ##   2) the submat (if using km)

      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  ) stop("unsupported distance.")
	if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
	}else{
	  if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
	  }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
	  }
	}
        attr( this_dist, "method" ) <- distance
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
            ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }

          this_dist <- sample_x$submat
        }
      }
    }

    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] )

    ##use samples for each k
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)

      }else if(clusterAlg=="kmdist"){
	this_assignment = kmeans(this_dist, k, iter.max = 10^9, nstart = 1, algorithm = c("Hartigan-Wong") )$cluster

      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        iterMax = 10^9
        nStart = 10
        #print(paste0("Seed: ", seedN))
        #print(paste0("max iteration = ", iterMax))
        #print(paste0("number of starts = ", nStart))
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = iterMax,
                                   nstart = nStart,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance,
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ##add to tally
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
    }
  }


  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  message("end fraction")
  return(res)
}


connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId

  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}



sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
##  to reduce memory overhead


  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )

  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

CDF=function(ml,breaks=100){
  #plot CDF distribution
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF",main="consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)

    #empirical CDF distribution. default number of breaks is 100
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)

    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
       thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
       bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.8,0.5,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)

  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
}


myPal = function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors = function(past_ct,ct,colorU,colorList){
	#description: sets common color of clusters between different K
	newColors = c()
	if(length(colorList)==0){
		#k==2
		newColors = colorU[ct]
		colori=2
	}else{
		newColors = rep(NULL,length(ct))
		colori = colorList[[2]]
		mo=table(past_ct,ct)
		m=mo/apply(mo,1,sum)
			for(tci in 1:ncol(m)){ # for each cluster
				maxC = max(m[,tci])
				pci = which(m[,tci] == maxC)
				if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
				#if new column maximum is unique, same cell is row maximum and is also unique
				##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
					newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
				}else{ #add new color.
					colori=colori+1
					newColors[which(ct==tci)] = colorU[colori]
				}
			}
	}
	return(list(newColors,colori,unique(newColors) ))
}

clusterTrackingPlot = function(m){
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and values are cluster assignments.
  plot(NULL,xlim=c(-0.1,1),ylim=c(0,1),axes=FALSE,xlab="samples",ylab="k",main="tracking plot")
  for(i in 1:nrow(m)){
    rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/nrow(m),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(i-1)/nrow(m),ncol(m)), col=m[i,],border=NA)
  }
  #hatch lines to indicate samples
  xl = seq(0,1-1/ncol(m),by=1/ncol(m))
  segments(  xl, rep(-0.1,ncol(m)) , xl, rep(0,ncol(m)), col="black")    #** alt white and black color?
  ypos = seq(1,0,by=-1/nrow(m))-1/(2*nrow(m))
  text(x=-0.1,y=ypos[-length(ypos)],labels=seq(2,nrow(m)+1,by=1))
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m


  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half

  fm = t(nm)+nm
  diag(fm) = diag(m)

  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]

  if(mode==1){
    return(vm) #vector
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }

}


rankedBarPlot=function(d,myc,cc,title){
	colors = rbind() #each row is a barplot series
	byRank = cbind()

	spaceh = 0.1 #space between bars
	for(i in 1:ncol(d)){
	  byRank = cbind(byRank,sort(d[,i],na.last=F))
	  colors = rbind(colors,order(d[,i],na.last=F))
	}
	maxH = max(c(1.5,apply(byRank,2,sum)),na.rm=T) #maximum height of graph

	#barplot largest to smallest so that smallest is in front.
	barp = barplot( apply(byRank,2,sum) ,  col=myc[colors[,1]] ,space=spaceh,ylim=c(0,maxH),main=paste("item-consensus", title),border=NA,las=1  )
	for(i in 2:nrow(byRank)){
	  barplot( apply(matrix(byRank[i:nrow(byRank),],ncol=ncol(byRank))  ,2,sum), space=spaceh,col=myc[colors[,i]],ylim=c(0,maxH), add=T,border=NA,las=1  )
	}
	xr=seq(spaceh,ncol(d)+ncol(d)*spaceh,(ncol(d)+ncol(d)*spaceh)/ncol(d)  )
	#class labels as asterisks
	text("*",x=xr+0.5,y=maxH,col=myc[cc],cex=1.4) #rect(xr,1.4,xr+1,1.5,col=myc[cc] )
}


########################################################################################################################################################################################################


variability = function(zscore_mat,info){
	library(clv)
	library(fields)

	mat = t(zscore_mat)
	mat = mat[order(rownames(mat)),]
	info = info[order(rownames(info)),]

	clusterStat = cls.scatt.data(mat, info[,2], dist = "euclidean")
	inter = clusterStat$intercls.complete


	medianInter = apply(inter,2,summary)
	medianInter = medianInter[2,]
	intra = clusterStat$intracls.complete

	minCount = floor(min(clusterStat$intracls.complete))
	maxCount = ceiling(max(max(clusterStat$intracls.complete), max(clusterStat$intercls.complete)))
	#problem = subset(medianInter, medianInter <= intra)
	pdf(paste0(as.numeric(Sys.time()),"clusterVariability.pdf"), width = 20, height = 15)
		boxplot(clusterStat$intercls.complete, ylim = c(minCount, maxCount))
		par(new = T)
		boxplot(clusterStat$intracls.complete, border = "red", ylim = c(minCount, maxCount))
	dev.off()
	rejectPro = ifelse(medianInter <= intra , 1, 0)
	if(max(rejectPro) == 1){
		reject = 1
		return(reject)
	} else{
		reject = 0
		return(reject)
	}
}


########################################################################################################################################################################################################


row.oneway.anova = function(Y,grplbl){
	ugrps<-unique(grplbl)
	ngrps<-length(ugrps) #number of groups
	ngenes<-dim(Y)[1]
	GrandM<-rowMeans(Y)   # overall mean

	SST<-rowSums((Y-GrandM)^2) # total sum of squares for each gene

	grp.mean<-matrix(NA,ngenes,ngrps)  # group mean matrix, rows for genes, each column for a different group
	grp.SSW<-matrix(NA,ngenes,ngrps)  # within-group sums of squares for each gene
	n<-rep(NA,ngrps)  # vector with group-specific sample sizes
	for (i in 1:ngrps)
	{
	grp.mtch<-(grplbl==ugrps[i])
	n[i]<-sum(grp.mtch)
	grp.mean[,i]<-rowMeans(Y[,grp.mtch])
	grp.SSW[,i]<-rowSums((Y[,grp.mtch]-grp.mean[,i])^2)
	}

	df1<-(ngrps-1)
	df2<-sum(n)-df1-1

	SSW<-rowSums(grp.SSW)
	SSB<-SST-SSW
	MSE<-SSW/df2
	MSB<-SSB/df1

	F.stat<-MSB/MSE
	pval<-1-pf(F.stat,df1,df2)
	res<-cbind.data.frame(stat=F.stat,pval=pval)
	res$FDR = p.adjust(res$pval, "BH", nrow(res))
	return(res)
}


########################################################################################################################################################################################################


sigGenes = function(selected_mat, info){
	aggr_FUN  <- mean
	combi_FUN <- function(x,y) "-"(x,y)
	pasteC <- function(x,y) paste(x,y,sep=" - ")

	fold = function(x, f, aggr_FUN = colMeans, combi_FUN = '-'){
		f = as.factor(f)
		i = split(1:nrow(x), f)
		x = sapply(i, function(i){ aggr_FUN(x[i,])})
		x = t(x)
		j = combn(levels(f), 2)
		ret = combi_FUN(x[j[1,],], x[j[2,],])
		rownames(ret) = paste(j[1,], j[2,], sep = '-')
		ret
	}

	mat = t(selected_mat)
	mat = subset(mat, rownames(mat) %in% rownames(info))
	mat = mat[order(rownames(mat)),]
	info = info[order(rownames(info)),]
	mat = t(mat)

	randClus = rep(0, nrow(info))
	info = cbind(info, randClus)
	clusFreq = as.data.frame(table(info[,2]))
	rIndex = sample(1:nrow(info),nrow(info), replace= FALSE)

	count = 1
	for(jh in 1:nrow(clusFreq)){
		jCount = (count + clusFreq[jh,2]) - 1
		freq = rIndex[count:jCount]
		count = (jCount+ 1)
		for(num in 1:length(freq)){
			cNum = freq[num]
			info[cNum,4] = jh
		}
	}
	info[,5] = paste0("grp", info[,4])
	rStatMat = row.oneway.anova(mat, factor(info[,5]))
	if(length(unique(info[,4])) < 3){
		matS = cbind(info[,4], t(mat))
		meanMat = aggregate(matS[,2:ncol(matS)], list(matS[,1]), mean)
		meanMat = t(meanMat)
		meanMat = meanMat[-1,]
		maxFC = meanMat[,1] - meanMat[,2]
		rStatMat = cbind(rStatMat, maxFC)
		rMinGenes = nrow(subset(rStatMat, (rStatMat[,3] < 0.01) & (rStatMat[,4] > 1)))
	} else {
		res = fold(t(mat), info[,4])
		res = t(res)
		maxFC = apply(res,1, max)
		rStatMat = cbind(rStatMat, maxFC)
		rMinGenes = nrow(subset(rStatMat, (rStatMat[,3] < 0.01) & (rStatMat[,4] > 1)))
	}
	statMat = row.oneway.anova(mat, factor(info[,3]))
	if(length(unique(info[,2])) < 3){
		matS = cbind(info[,2], t(mat))
		meanMat = aggregate(matS[,2:ncol(matS)], list(matS[,1]), mean)
		meanMat = t(meanMat)
		meanMat = meanMat[-1,]
		maxFC = meanMat[,1] - meanMat[,2]
		statMat = cbind(statMat, maxFC)
		minGenes = nrow(subset(statMat, (statMat[,3] < 0.01) & (statMat[,4] > 1)))
		return(minGenes)
	}
	res = fold(t(mat), info[,3])
	res = t(res)
	maxFC = apply(res,1, max)
	statMat = cbind(statMat, maxFC)
	minGenes = nrow(subset(statMat, (statMat[,3] < 0.01) & (statMat[,4] > 1)))
	return(minGenes)
}


########################################################################################################################################################################################################


stability = function(cluster_file, selected_cluster, output_dir){
	maxClus = max(cluster_file[,2])
	selected_mat_file = paste(output_dir,"/",output_dir,".k=",selected_cluster,".consensusMatrix.csv",sep="")
	mat_file = read.csv(selected_mat_file)
	mat_file = mat_file[,-1]
	colnames(mat_file) = rownames(cluster_file)
	rownames(mat_file) = rownames(cluster_file)
	stabMat = matrix(NA, ncol = 2, nrow = maxClus)

	for(i in 1:maxClus){
		subClus = subset(cluster_file, cluster_file[,2] != i)
		subClus2 = subset(cluster_file, cluster_file[,2] == i)
		subMat = subset(mat_file, (rownames(mat_file) %in% rownames(subClus2)))
		subMat = t(subMat)
		subMat = subset(subMat, rownames(subMat) %in% rownames(subClus))

		subMat2 = subset(mat_file, (rownames(mat_file) %in% rownames(subClus2)))
		subMat2 = t(subMat2)
		subMat2 = subset(subMat2, rownames(subMat2) %in% rownames(subClus2))

		stabMat[i,1] = mean(subMat2)
		stabMat[i,2] = mean(subMat)
	}
	colnames(stabMat) = c("sameClus", "otherClus")
	rownames(stabMat) = c(1:maxClus)
	return(stabMat)
}


########################################################################################################################################################################################################


CDF_RCC = function(zscore_mat,output_dir,maxK,times,rIN, selected_mat){
	zscore_mat = as.matrix(zscore_mat)
	maxK = maxK
	color = c("blue", "red", "yellow", "green", "hotpink", "orange", "brown", "purple", "cyan")

	output_dir = paste0(output_dir, "_",times)
	diffSlope = 0.0349066	#difference between slopes of 2 different Ks

	if(times == 1){
		pItem = 0.6
		pFeature = 0.8}

	if(times == 2){
		pItem = 0.6
		pFeature = 1}

	if(times == 3){
		pItem = 0.7
		pFeature = 0.8}

	if(times == 4){
		pItem = 0.7
		pFeature = 1}

	if(times == 5){
		pItem = 0.8
		pFeature = 0.8}

	if(times == 6){
		pItem = 0.8
		pFeature = 1}

	if(times == 7){
		pItem = 0.9
		pFeature = 0.8}

	if(times == 8){
		pItem = 0.9
		pFeature = 1}

	interval = 5
	sClusVal = 0.8
	oClusVal = 0.2
	tempClus = NULL
	allSlopes = NULL
	allLines = NULL
	sVal = NULL
	yVal = NULL
	fdrGenes = NULL
	minGenesR = round(((nrow(selected_mat) * minGenesP)/100),0)

	randSeed = (runif(1, min = 1, max = 10^(rIN[times])))/10000
	seeds = randSeed
	results = ConsensusClusterPlus(zscore_mat,maxK=maxK,
	reps=repCount,pItem=pItem,
	pFeature=pFeature,title= output_dir,clusterAlg=alg,
	innerLinkage= innerLinkage, finalLinkage=innerLinkage,
	distance=dis,plot="png",writeTable=T, seed = randSeed, verbose = F)
	clusSlopes = NULL

	files1 = list.files(output_dir, pattern = "*Class.csv", full.names=T)
	files2 = list.files(output_dir, pattern = "*atrix.csv", full.names=T)
	png(paste0(as.numeric(Sys.time()),times,".png"))
	for(i in 1:length(files1)){
		print(i)
		data = as.data.frame(fread(files2[i]))
		rownames(data) = data[,1]
		data = data[,-1]
		m = upperTriangle(data, diag = F, byrow = T)
		mat = as.vector(m)

		info = read.csv(files1[i], header = F)
		rownames(info) = info[,1]
		a = as.data.frame(table(info[,2]))

		s = sum((a[,2]*(a[,2] - 1)/2))
		s = 1 - (s/length(mat))
		s1 = s - 0.05
		s2 = s + 0.05

		cdf = 0
		cdfP = NULL
		grp_counts = NULL

		for(c in 0:100){
			cutoff  = c/100
			grp_counts = append(grp_counts,cutoff)
			cdfC = (length(subset(mat, mat <= cutoff)))/(length(mat))
			cdf = cdf + cdfC
			cdfP = append(cdfP, cdfC)
		}
		grp_mids = cdfP
		par(new = T)
		plot(cdfP, type = "l", col = color[i], xlim = c(1,100), ylim = c(0,1))
		abline(h = s1, col = color[i])
		abline(h = s2, col = color[i])

		grps = cbind(grp_counts, grp_mids)
		rownames(grps) = 1:101

		hGrps = subset(grps, (grps[,2] >= s1) & (grps[,2] <= s2))
		print(paste0("nrows:", nrow(hGrps)))
		Saxis = floor(s*100)
		if(nrow(hGrps) < 20){ next }
		print(dim(hGrps))
		start_point_I = as.numeric(rownames(hGrps)[1])
		end_point_I = as.numeric(rownames(hGrps)[nrow(hGrps)])

		counts = hGrps
		lm_r = lm(counts[,2] ~ counts[,1])
		slope = lm_r[[1]][[2]]

		choosen = 1
		start_point = start_point_I
		end_point = end_point_I
		if(slope > threshold){
			choosen = 0
			for(sp in start_point_I:(end_point_I - 20)){
				ep = sp + 20
				counts = grps[c(ep:sp),]
				lm_r = lm(counts[,2] ~ counts[,1])
				slope = lm_r[[1]][[2]]
				if(slope <= threshold){
					start_point = sp
					end_point = ep
					choosen = 1
					break
				}
			}
		}
		if(choosen == 0){ next }

		counts = grps[c(end_point:start_point),]
		lm_r = lm(counts[,2] ~ counts[,1])
		slope = lm_r[[1]][[2]]
		initial_slope = slope
		if(slope > threshold){
			next
		} else {
			for(k in start_point_I:end_point_I){
				end_point = end_point + interval
				if(end_point > end_point_I){
					end_point = end_point - interval
					k = end_point_I + 1
					break
				}
				counts = grps[c(start_point:end_point),]
				lm_r = lm(counts[,2] ~ counts[,1])
				new_slope = lm_r[[1]][[2]]
				if(abs(initial_slope - new_slope) > diffSlope){
					end_point = end_point - interval
					k = end_point_I + 1
					break
				} else {
						initial_slope = new_slope
				}
			}
			finalEP = end_point

			for(k in start_point_I:end_point_I){
				start_point = start_point - interval
				if(start_point < start_point_I){
					start_point = start_point + interval
					k = end_point_I + 1
					break
				}
				counts = grps[c(start_point:end_point),]
				lm_r = lm(counts[,2] ~ counts[,1])
				new_slope = lm_r[[1]][[2]]
				if(abs(initial_slope - new_slope) > diffSlope){
					start_point = start_point + interval
					k = end_point_I + 1
					break
				} else {
					initial_slope = new_slope
				}
			}
			finalSP = start_point
			if(finalSP > 30){ next }
			counts = grps[c(finalSP:finalEP),]
			lm_r = lm(counts[,2] ~ counts[,1])
			FS = lm_r[[1]][[2]]
			line_length = finalSP - finalEP

			info[,3] = paste0("grp", info[,2])
			freqClusInfo = as.data.frame(table(info[,2]))
			freqClusInfo = subset(freqClusInfo, freqClusInfo[,2] >= 3)
			clusters = freqClusInfo[,1]
			if(length(clusters) < 2) {
				minGenes = 0
			} else {
				sink(paste0(as.numeric(Sys.time()),"_",output_dir,"_cluster",i,".txt"))
				info = subset(info, info[,2] %in% clusters)
				minGenes = sigGenes(selected_mat, info)
				sink()
			}
			a = strsplit(files1[i], "=")[[1]][2]
			fdrGenes = append(fdrGenes, minGenes)
			tempClus = append(tempClus, as.numeric(strsplit(a, "[.]")[[1]][1]))
			allSlopes = append(allSlopes, FS)
			sVal = append(sVal,s)
			yVal = append(yVal, grps[as.character(Saxis),2])
			allLines = append(allLines, line_length)
		}
	}
	dev.off()

	if(is.null(tempClus)){
		Fclus = 0
		return(Fclus)
	}

	matStat = NULL
	allLines = abs(allLines)
	matStat = cbind(tempClus, allSlopes)
	matStat = cbind(matStat, allLines)
	matStat = cbind(matStat, sVal)
	matStat = cbind(matStat, yVal)


	weightAll = NULL

	for(rows in 1:nrow(matStat)){
		weight = 0
		if(matStat[rows,2] <= 0.0872665){ weight = weight + 1 }
		if(matStat[rows,3] >= 40){ weight = weight + 1 }
		weightAll = append(weightAll, weight)
	}
	matStat = cbind(matStat, weightAll)

	sClusAll = NULL
	oClusAll = NULL
	for(mClus in 1:nrow(matStat)){
		selected_cluster = matStat[mClus, 1]
		selected_cluster_file = paste(output_dir,"/",output_dir,".k=",selected_cluster,".consensusClass.csv",sep="")
		cluster_file = read.csv(selected_cluster_file,header = FALSE)
		matFC = stability(cluster_file, selected_cluster, output_dir)
		sClus = summary(matFC[,1])[2]
		sClusAll = append(sClusAll, sClus)
		oClus = summary(matFC[,2])[2]
		oClusAll = append(oClusAll, oClus)
	}
	matStat = cbind(matStat, sClusAll)
	matStat = cbind(matStat, oClusAll)
	matStat = cbind(matStat, fdrGenes)
	colnames(matStat) = c("tempClus","allSlopes", "allLines", "sVal", "yVal", "weight", "sClus", "oClus", "fdrGenes")
	write.csv(matStat, paste0(as.numeric(Sys.time()),times,"CDF.csv"))

	matStat_check = subset(matStat, (matStat[,2] < threshold) & (matStat[,3] >= min_line_len) & (matStat[,7] >= sClusVal) & (matStat[,8] <= oClusVal) & (matStat[,9] >= minGenesR))
	if(nrow(matStat_check) < 1){
		Fclus = 0
		return(Fclus)
	} else {
		matStat = matStat_check
		clus = which(matStat[,6] == max(matStat[,6]))
		print(paste0("clus: ", clus))
		if(length(clus) == 1){
			Fclus = matStat[clus,1]
			return(Fclus)
		}
		matStat = matStat[clus,]
		Fclus = matStat[,1]
		return(Fclus)
	}
}


########################################################################################################################################################################################################


RCC_clus = function(config_file){

	library(matrixStats)
	library(pforeach)
	library(data.table)
	library(cluster)
	library(clue)
	library(ComplexHeatmap)
	library(circlize)
	library(clv)
	library(fields)
	library(gdata)

  conf_file = read.csv(config_file,header = FALSE)
  expr_data = as.data.frame(fread(as.character(conf_file[1,2])))
  rownames(expr_data) = expr_data[,1]
  expr_data = expr_data[,-1]
  ssgsea_data = expr_data
  expr_data = t(expr_data)
  print("matrix file read")

  sample_ids = rownames(expr_data)
  cluster_num = rep(1,nrow(expr_data))
  sample_ids = cbind(sample_ids,cluster_num)
  rownames(sample_ids) = sample_ids[,1]
  sample_ids = as.data.frame(sample_ids)
  print("SampleInfo file read")

  sampleInfo = NULL
  m = NULL

  sampleInfo = as.data.frame(fread(as.character(conf_file[2,2])))
  rownames(sampleInfo) = sampleInfo[,1]
  sampleInfo = sampleInfo[,-1]
  init_cols = colnames(sampleInfo)
  init_colNum = ncol(sampleInfo)

  common_samples = intersect(rownames(sample_ids),rownames(sampleInfo))
  sampleInfo = subset(sampleInfo,rownames(sampleInfo) %in% common_samples)
  sampleInfo_col_num = ncol(sampleInfo)


  tcol_num = ncol(sampleInfo) + 1

  alg <<- "km"		#algorithm used for concensus clustering
  dis <<- "euclidean"		#distance matrix to be used
  output_dir = "RCCs"
  repCount <<- 100
  innerLinkage <<- "ward.D2"
  corUse <<- "everything"
  threshold <<- as.numeric(as.vector(conf_file[3,2]))
  min_samples = as.numeric(as.vector(conf_file[4,2]))
  min_line_len <<- as.numeric(as.vector(conf_file[5,2]))
  var_genes_percent = as.numeric(as.vector(conf_file[6,2]))
  minGenesP <<- as.numeric(as.vector(conf_file[7,2]))
  dataset = as.character(conf_file[8,2])
  finalOut <<- as.character(conf_file[9,2])
  zscore_mat = NULL
  colnames(sampleInfo) = init_cols

	setwd(finalOut)
	recursive_consensus = function(sample_ids,col_num,threshold,output_dir){
		system("rm -r RCC*")
		maxK = min(10,ceiling(nrow(sample_ids)/10))

		if((nrow(sample_ids) < min_samples)){
			selected_cluster = 0
		} else {
			if(maxK < 3){
				maxK = 3
			}
			expr_subset = NULL
			expr_subset = subset(expr_data,rownames(expr_data) %in% rownames(sample_ids))
			expr_subset = as.matrix(t(expr_subset))
			var_col = apply(expr_subset, 1, var)
			var_mat = NULL
			var_mat = cbind(expr_subset,var_col)
			colnames(var_mat)[ncol(var_mat)] = "Variance"

			ord_var_mat = var_mat[order(var_mat[,ncol(var_mat)],decreasing = TRUE),]
			genes_selected = (nrow(ord_var_mat) * var_genes_percent)/100
			genes_selected = round(genes_selected,0)
			if(dataset == "bulk"){
				if(genes_selected < 500){
					genes_selected = min(nrow(expr_subset),500)
				}
			}
			selected_mat <<- ord_var_mat[c(1:genes_selected),-c(ncol(ord_var_mat))]

			zscore_mat <<- (selected_mat - rowMeans(selected_mat))/(rowSds(as.matrix(selected_mat)))[row(selected_mat)]
			rIN = round(runif(8, min =4, max = 12),0)
			iterations = 8
			cluster_parallel <<- pforeach(times = 1:iterations, .combine = append, .parallel = T) ({
				cluster_slopes1 = CDF_RCC(zscore_mat,output_dir,maxK, times, rIN, selected_mat)
			})

			if(length(cluster_parallel) == 0){
				selected_cluster = 0
			} else if(is.null(cluster_parallel)){
				selected_cluster = 0
			} else if((length(cluster_parallel) == 1)){
			 if(unique(cluster_parallel) == 0){
			   selected_cluster = 0
			 }
			} else {
				clusFreq = as.data.frame(table(cluster_parallel))
				maxFreq = which(clusFreq[,2] == max(clusFreq[,2]))
				if(length(maxFreq) > 1){
					selected_cluster = clusFreq[maxFreq,1]
					selected_cluster = as.numeric(as.vector(selected_cluster))
					if(max(clusFreq[,2]) > 1){
						selected_cluster = max(selected_cluster)
					} else {
						selected_cluster = min(selected_cluster)
					}
				} else { selected_cluster = clusFreq[maxFreq,1] }
			}
		}
		selected_cluster = as.numeric(as.character(selected_cluster))
		timeSelected = which(cluster_parallel == selected_cluster)
		timeSelected = timeSelected[1]


		if(selected_cluster > 0){
			output_dir = paste0(output_dir,"_",timeSelected)

			zscore_file_name = paste("2018",runif(1, min = 0, max = 500),"colVar",selected_cluster,".csv",sep="")
			write.csv(rownames(zscore_mat),file = zscore_file_name,row.names = F)
			selected_cluster_file = paste(output_dir,"/",output_dir,".k=",selected_cluster,".consensusClass.csv",sep="")
			cluster_file <<- read.csv(selected_cluster_file,header = FALSE)
			rownames(cluster_file) = cluster_file[,1]
			cluster_file[,3] = paste0("grp", cluster_file[,2])

			for(h in 1:nrow(cluster_file)){
				sample_name = as.character(rownames(cluster_file)[h])
				sampleInfo[sample_name,col_num] <<- cluster_file[sample_name,2]
			}

			sampleLevel = abs(sampleInfo_col_num - col_num)

			fileName = paste0("Level", sampleLevel)

			if(sampleLevel > 1){
				for(g in 1:(sampleLevel - 1)){
					sampleSel = as.character(rownames(cluster_file)[1])
					a = unique(sampleInfo[sampleSel, sampleInfo_col_num + g])
					fileName = paste(fileName,"_k",a, sep = "")
				}
			}

			if(fileName == "Level1"){
				matFC = stability(cluster_file, selected_cluster, output_dir)
				write.csv((matFC), file  = paste("stability_", fileName, ".csv", sep = ""))
			} else {
				matFC = stability(cluster_file, selected_cluster, output_dir)
				write.csv((matFC), file  = paste("stability_", fileName, ".csv", sep = ""))
			}
			m = col_num + 1

			for(b in 1:selected_cluster){
				sample_ids = subset(cluster_file,cluster_file[,2] == b)
				recursive_consensus(sample_ids,m,threshold,output_dir)
			}
		}
	}

	recursive_consensus(sample_ids,tcol_num,threshold,output_dir)

	number_of_levels = abs(sampleInfo_col_num - ncol(sampleInfo))
	sampleInfo[,c((ncol(sampleInfo)-number_of_levels):ncol(sampleInfo))][is.na(sampleInfo[,c((ncol(sampleInfo)-number_of_levels):ncol(sampleInfo))])] = 0
	colnames(sampleInfo)[(sampleInfo_col_num + 1): ncol(sampleInfo)] = paste0("Level_",1:number_of_levels)
	colNum = c(1:number_of_levels)
	cols = paste("Level_", colNum, sep = "")

	if(length(cols) == 1){
		Clusters = sampleInfo[,ncol(sampleInfo)]
		sampleInfo = cbind(sampleInfo,Clusters)
		colnames(sampleInfo)[ncol(sampleInfo)] = "Clusters"
	} else {
		colnames(sampleInfo)[c(((ncol(sampleInfo) - number_of_levels)+1):ncol(sampleInfo))] = cols
		sampleInfo[,(ncol(sampleInfo) + 1)] = do.call(paste0, sampleInfo[c(cols)])
		colnames(sampleInfo)[ncol(sampleInfo)] = "Concatenated_clusters"

		len = length(unique(sampleInfo[,ncol(sampleInfo)]))
		cluster_found = unique(sampleInfo[,ncol(sampleInfo)])

		for(k in 1:len){
			for(j in 1:nrow(sampleInfo)){
				if(sampleInfo[j,ncol(sampleInfo)] == cluster_found[k]){
					sampleInfo[j,ncol(sampleInfo)] = k
				}
			}
		}

		colnames(sampleInfo)[ncol(sampleInfo)] = "Clusters"
		sampleInfo <<- sampleInfo
	}
	write.csv(sampleInfo,file= "OutputRCC.csv")

	system("rm -rf RCC*")
	system("cat *Var*csv > genesUsed.csv")
	system("rm -rf lm* *cluster*txt *png *CDF.csv stability* *Var* *mds*")
}
