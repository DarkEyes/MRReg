# mixture model support functions
getRMSEFromMixtureModel<-function(x,y,k,polyDegree = 1, expFlag = FALSE)
{
  df<- data.frame(y,x)
  N<-length(y)
  df<- data.frame(y,x)

  if(polyDegree == 1 && expFlag == FALSE)
    mixmodel<-flexmix(y~.,data = df, k=k, control = list(minprior=0.02) )
  else
    mixmodel<-flexmix(formula = makePolyFormula(df,degree = polyDegree, expFlag= expFlag),data = df, k=k, control = list(minprior=0.02) )

  mixmodel<-flexmix(y~.,data = df, k=k, control = list(minprior=0.02) )
  yh<-predict(mixmodel, data.frame( (x) ) )
  clsVec<-mixmodel@cluster
  residuals<-list()
  for(i in seq(1,N))
  {
    cls<-clsVec[i]
    currResidual<- y[i] - yh[[cls]][i]
    residuals<-append(residuals,currResidual)
  }
  RMSE<-sqrt(mean(as.numeric(residuals)^2))
  return(list("clsVec"=clsVec, "RMSE"=RMSE) )
}

#mixOut<-getRMSEFromMixtureModel(x,y,k=4)
#TrueClsVec<-DataT$cls
#PredClsVec<-mixOut$clsVec
getMixturePartitionFscore<-function(TrueClsVec,PredClsVec)
{
  # TrueCopt[k,j]

  TrueclsMembers<-unique(TrueClsVec)
  M<-length(TrueclsMembers)
  Fscore<-0
  TP<-0
  FP<-0
  FN<-0
  for(trueCls in TrueclsMembers)
  {
    inxVec<-TrueClsVec == trueCls
    predCls<-median(PredClsVec[inxVec])
    predInxVec<- PredClsVec == predCls
    TP<- TP+ sum( predInxVec & inxVec)
    FP<- FP+sum( predInxVec & ! inxVec)
    FN<- FN+ sum( !predInxVec & inxVec)
  }
  prcVal<-TP/(TP+FP)
  recall<-TP/(TP+FN)
  Fscore<-2*(prcVal*recall)/(recall+prcVal)
  return(list("prcVal"=prcVal,"recall"=recall,"Fscore"=Fscore) )
}

makePolyFormula<-function(df,degree=2, expFlag= FALSE)
{
  varnames <- paste("X", 1:(length(df)-1), sep="")
  for(i in seq(length(varnames)))
  {
    if(expFlag == FALSE)
      varnames[i]<-sprintf("%s^%d",varnames[i],degree)
    else
      varnames[i]<-sprintf("exp(%s)",varnames[i],degree)
  }

  polyFormula <- as.formula(paste("y ~ ", paste(varnames, collapse= "+")))
  return(polyFormula)
}
# support function
getResidualFromCopt<-function(Copt,models)
{
  M<-dim(Copt)[1]
  residuals<-list()
  for(i in seq(1,M))
  {
    residuals<-append(residuals,models[[Copt[i,1]]][[Copt[i,2]]]$residuals)
  }
  return(list("residuals"=as.numeric(residuals) ) )
}
# ======== Evaluation parts for Copt
getPartitionFscore<-function(TrueCopt,Copt,clsLayer)
{
  # TrueCopt[k,j]

  TrueCopt<-matrix(TrueCopt,ncol=2)
  #Copt<-matrix(Copt,ncol=4)
  N1<-dim(TrueCopt)[1]
  N2<-dim(Copt)[1]
  TrueFlag<-logical(N1)
  PredFlag<-logical(N2)
  Fscore<-0
  TP<-0
  FP<-0
  FN<-0

  for(i1 in seq(1,N1))
  {
    currTrCls<-TrueCopt[i1,1:2]
    for(i2 in seq(1,N2))
    {
      currPdCls<-Copt[i2,1:2]
      if( sum(currTrCls ==currPdCls) ==2 ) # true cls is in pred cls set
      {
        TrueFlag[i1] = TRUE
        PredFlag[i2] = TRUE
        break
      }
    }
  }
  # Find true positive and false negative
  for(i1 in seq(1,N1))
  {
    currTrCls<-TrueCopt[i1,1:2]
    if(TrueFlag[i1]==TRUE)
    {
      TP<-TP+ sum(clsLayer[,currTrCls[1]] == currTrCls[2])
    }
    else
    {
      FN<-FN+ sum(clsLayer[,currTrCls[1]] == currTrCls[2])
    }
  }
  # Find  false positive
  for(i2 in seq(1,N2))
  {
    currPdCls<-Copt[i2,1:2]
    # print(currPdCls)
    if(PredFlag[i2]==FALSE)
    {
      FP<-FP+ sum(clsLayer[,currPdCls[1]] == currPdCls[2])
    }
  }

  prcVal<-TP/(TP+FP)
  recall<-TP/(TP+FN)
  Fscore<-2*(prcVal*recall)/(recall+prcVal)
  return(list("prcVal"=prcVal,"recall"=recall,"Fscore"=Fscore))
}

# ====== Based line: Greedy algorithm
greedyAlgo<-function(DataT,out)
{

  N<-length(DataT$Y)
  flagVec<-logical(N)
  nL<-dim(DataT$clsLayer)[2]
  Vc <-cbind(numeric(N),numeric(N))
  mSquareErr<-list()
  sortVec<-list()
  sortCls<-list()
  l<-1
  for(k in seq(1,nL) )
  {
    currLayer<-unique(DataT$clsLayer[,k])
    layerMSquareErr<-list()
    for(j in seq(1,length(currLayer) ) )
    {
      sortVec[[l]]<-mean(out$models[[k]][[j]]$residuals^2)
      sortCls[[l]]<-c(j,k)
      l=l+1
    }

  }
  sortVec<-as.numeric(sortVec)
  orderVec<-order(sortVec)

  for(i in orderVec)
  {
    currCls<-sortCls[[i]]
    mark<-unique(DataT$clsLayer[,currCls[2]])
    inxVec<-DataT$clsLayer[,currCls[2]] == mark[currCls[1] ]
    if(sum(flagVec[inxVec])==0)
    {
      flagVec[inxVec]<-TRUE
      k<-currCls[2] # Layer
      j<-currCls[1] # Cls in Layer
      Vc[inxVec,1]<-k
      Vc[inxVec,2]<-j

      print(sprintf("Greedy: Layer%d, cls%d",k,j) )
      if(is.null(out$models[[k]][[j]]$R2cv))
      {
        inxFilterVec<-DataT$clsLayer[,k] == j # selecting only members of jth cluster in ith layer
        out$models[[k]][[j]]$R2cv<- crossVal10FoldEstFunc(DataT$X[inxFilterVec,],DataT$Y[inxFilterVec])$r2
      }
    }
  }
  CoptGreedy<-unique(Vc)
  return(list("Copt"=CoptGreedy,"out"=out))
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
