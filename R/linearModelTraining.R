#'@title  linearModelTraining
#'
#'@description
#' linearModelTraining is a support function for training linear models for partitions in all layers.
#'
#'@param DataT contains a multiresolution dataset s.t.
#' \code{DataT$X[i,d]} is a value of feature \code{d} of individual \code{i},
#' \code{DataT$Y[i]} is value of target variable of individual \code{i} that we want to fit \code{DataT$Y ~ DataT$X} in linear model, and
#' \code{clsLayer[i,j]} is a cluster ID of individual \code{i} at layer \code{j}; \code{clsLayer[i,1]} is the first layer that everyone typically belongs to a single cluster.
#'@param insigThs is a threshold to determine whether a magnitude of a feature coefficient is enough so that the feature is designated as a selected feature.
#'@param alpha is a significance level to determine whether a magnitude of a feature coefficient is enough so that the feature is designated as a selected feature.
#'@param messageFlag is a flag. If it is true, the function shows the text regarding the progress of computing.
#'@param polyDegree is a degree of polynomial function that is used to fit the data.
#'If it is greater than 1, the polynomial formula is used in \code{lm()} instead of \code{"y=."}.
#'@param expFlag is an exponential flag to control the formula for data fitting.
#'If it is true, then the exp() formula is used in \code{lm()} instead of \code{"y=."}.
#'
#'@return This function returns \code{models} and \code{DataT}.
#'
#' \item{ \code{models[[j]][[k]]} }{ is a linear model of a cluster ID \code{k} at the layer \code{j}.
#'  The \code{models[[j]][[k]]$selFeatureSet} represents a set of selected-feature indices of the model where the feature index 1 is the intercept,
#'   and the feature index \code{d} is the (d-1)th variable \code{DataT$X[,d-1]}. }
#' \item{ DataT }{ is a \code{DataT} with \code{DataT$nNodes}, which is a number of total models from all layers. }
#'
#'@examples
#'# Running linearModelTraining using simulation data
#' DataT<-SimpleSimulation(100,type=1)
#' obj<-linearModelTraining(DataT)
#'
#'@importFrom caret trainControl train
#'@importFrom stats cor lm predict rnorm runif as.formula
#'@export
#'
linearModelTraining<-function(DataT,insigThs=1e-8,alpha=0.05,messageFlag=FALSE,polyDegree = 1, expFlag= FALSE)
{

  out1<-remappingClusterInx(DataT$clsLayer)
  DataT$clsLayer<-out1$nclsLayer
  DataT$clsNameMappingTable<-out1$clsNameMappingTable
  out<-dim(DataT$clsLayer)
  N<-out[[1]]
  nL<-out[[2]]
  models<-list()
  residualMat<-matrix(0,N,nL)
  nNodes <-0
  IDk<-1
  for(inx in seq(1,nL))
  {
    currLayer<- DataT$clsLayer[,inx]
    currLayerList<-unique(currLayer)
    nCls<- length( currLayerList )
    nNodes<-nNodes + nCls
    submodels<-list()
    for(inx2 in seq(1,nCls))
    {
      inxFilterVec<- currLayer == currLayerList[inx2]
      x<-DataT$X[inxFilterVec,]
      y<-DataT$Y[inxFilterVec,]
      df = data.frame(y, x)
      #===========================st Linear model
      if(polyDegree == 1 && expFlag == FALSE)
        submodels[[inx2]] <- lm(y ~ ., data = df)
      else
        submodels[[inx2]] <- lm(formula = makePolyFormula(df,degree = polyDegree, expFlag= expFlag), data = df)
      submodels[[inx2]] $ID<-IDk
      IDk<-IDk +1
      #===========================fn Linear model
      if(inx<nL)
        submodels[[inx2]]$ChildrenCls<- unique(DataT$clsLayer[inxFilterVec,inx+1])
      else
        submodels[[inx2]]$ChildrenCls<- NULL

      if(inx>1)
      {
        submodels[[inx2]]$ParentCls<- unique(DataT$clsLayer[inxFilterVec,inx-1])
      }
      else
        submodels[[inx2]]$ParentCls<- NULL

      submodels[[inx2]]$optFlag<-FALSE

      SelectedFeatures<- summary(submodels[[inx2]])$coefficients[,4] <= alpha
      sigFeature<- abs(summary(submodels[[inx2]])$coefficient[,1])>=insigThs

      SelectedFeatures<- SelectedFeatures & sigFeature
      selFeatureSet<-1:length(SelectedFeatures)
      selFeatureSet<- selFeatureSet[SelectedFeatures]
      submodels[[inx2]]$selFeatureSet<-selFeatureSet

      if(length(y)>1)
        residualMat[inxFilterVec,inx] <- submodels[[inx2]]$residuals
      else
        residualMat[inxFilterVec,inx] <-sum(submodels[[inx2]]$residuals)

      submodels[[inx2]]$clsName<-sprintf("C%g",DataT$clsNameMappingTable[[inx]][[inx2]])
      if(messageFlag == TRUE)
        message(sprintf("Training -> Layer%d,Cls:%s",inx,submodels[[inx2]]$clsName))
    }

    models[[inx]]<-submodels
  }
  DataT$nNodes<-nNodes
  return(list("models"=models,"DataT"=DataT))
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
