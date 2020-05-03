#'
#'
#'
#'@export
#'
linearModelTraining<-function(DataT,insigThs=1e-8,alpha=0.05)
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
      submodels[[inx2]] <- lm(y ~ ., data = df)
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

      SelectedFeatures<- summary(submodels[[inx2]])$coefficients[,4] < alpha
      sigFeature<- abs(summary(submodels[[inx2]])$coefficient[,1])>insigThs

      SelectedFeatures<- SelectedFeatures & sigFeature
      selFeatureSet<-1:length(SelectedFeatures)
      selFeatureSet<- selFeatureSet[SelectedFeatures]
      submodels[[inx2]]$selFeatureSet<-selFeatureSet

      if(length(y)>1)
        residualMat[inxFilterVec,inx] <- submodels[[inx2]]$residuals
      else
        residualMat[inxFilterVec,inx] <-sum(submodels[[inx2]]$residuals)

      submodels[[inx2]]$clsName<-sprintf("C%g",DataT$clsNameMappingTable[[inx]][[inx2]])
      print(sprintf("Layer%d,Cls:%s",inx,submodels[[inx2]]$clsName))
    }

    models[[inx]]<-submodels
  }
  DataT$nNodes<-nNodes
  return(list("models"=models,"DataT"=DataT))
}
