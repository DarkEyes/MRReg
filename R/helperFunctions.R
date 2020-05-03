remappingClusterInx<-function(clsLayer)
{
  nL<- dim(clsLayer)[2]
  nIdv<-dim(clsLayer)[1]
  nclsLayer<-matrix(1, nIdv, nL)
  clsNameMappingTable<-list()
  for(i in seq(1,nL))
  {
    currLayer<-clsLayer[,i]
    currLayerList<-unique(currLayer)
    #print(sprintf("rmmfunc L%d",i) )
    #print(currLayerList)
    for( j in seq(1,length(currLayerList)))
    {
      targetInxvec <-currLayer== currLayerList[j]
      nclsLayer[targetInxvec,i] <-j
    }
    clsNameMappingTable[[i]]<-currLayerList
  }

  return(list("nclsLayer"=nclsLayer,"clsNameMappingTable"=clsNameMappingTable))
}

# MDL
getRealStoreBits<-function(realVec, insigThs=1)
{
  normalizedRealVec<-realVec
  normalizedRealVec[abs(normalizedRealVec)==insigThs]<-2 # for val = 1 need 2 bits
  normalizedRealVec[abs(normalizedRealVec)<insigThs]<-1
  realVecBits <- sum( ceiling(log2(abs(normalizedRealVec)) )+1)
  return(list("realVecBits"=realVecBits))
}

getModelInfoRecRatio<-function(h1ResidualVec,h2ResidualVec,h1coeff,h2coeff)
{
  h1ResBits<-getRealStoreBits(h1ResidualVec)$realVecBits
  h2ResBits<-getRealStoreBits(h2ResidualVec)$realVecBits
  h1modelBits <- getRealStoreBits(h1coeff)$realVecBits
  h2modelBits <- getRealStoreBits(h2coeff)$realVecBits
  h1TotalBits <- (h1ResBits+h1modelBits )
  h2TotalBits <- (h2ResBits+h2modelBits )
  modelInfoRecBitsRatio<-(  h1TotalBits -  h2TotalBits )/h1TotalBits
  return(list("modelInfoRecRatio"=modelInfoRecBitsRatio,"h1ResBits"=h1ResBits,"h2ResBits"=h2ResBits,"h1modelBits"=h1modelBits,"h2modelBits"=h2modelBits))
}

# CV functions
crossValEstFunc<-function(X,Y,clsInxVec)
{
  # X is a vector of independent vars, Y is a vector of a dependent var.
  # clsInxVec is a vector of clustering indices of rows in X and Y.
  # clsInxVec is always the local index vector!

  # return Coefficient of Determination R^2
  clsVec<-unique(clsInxVec)
  predY<-numeric(length(Y))

  #print(length(clsVec))
  for(cls in clsVec)
  {
    # cls is the validate cluster and the rest are training part
    currClsInx<-clsInxVec == cls
    currReverseClsInx <- !currClsInx # training

    #print(sprintf("Cls%d, Length%d",cls,sum(currClsInx)))

    x<-X[currReverseClsInx,]
    y<-Y[currReverseClsInx]

    df = data.frame(y, x)


    model1<- lm(y ~ ., data = df)

    if(sum(currClsInx)!=1)
    {
      predYtmp<-predict(model1, data.frame(X[currClsInx,]) )
    }
    else
      predYtmp<-0

    predY[currClsInx]<-predYtmp
  }
  return(list("r2"=cor(Y,predY)^2))
}

crossVal10FoldEstFunc<-function(X,Y)
{
  # X is a vector of independent vars, Y is a vector of a dependent var.
  # clsInxVec is a vector of clustering indices of rows in X and Y.
  # clsInxVec is always the local index vector!

  # return Coefficient of Determination R^2

  train.control <- trainControl(method = "cv", number = 10)

  x<-X
  y<-Y

  df = data.frame(y, x)


  model1 <- train(y ~., data = df, method = "lm",
                  trControl = train.control)

  return(list("r2"=model1$results$Rsquared))
}
