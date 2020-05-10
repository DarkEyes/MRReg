#'@title  FindMaxHomoOptimalPartitions
#'
#'@description
#' FindMaxHomoOptimalPartitions is a main function for inferring optimal homogeneous clusters from a multiresolution dataset \code{DataT}.
#'
#'@param DataT contains a multiresolution dataset s.t.
#' \code{DataT$X[i,d]} is a value of feature \code{d} of individual \code{i},
#' \code{DataT$Y[i]} is value of target variable of individual \code{i} that we want to fit \code{DataT$Y ~ DataT$X} in linear model, and
#' \code{clsLayer[i,j]} is a cluster ID of individual \code{i} at layer \code{j}; \code{clsLayer[i,1]} is the first layer that everyone typically belongs to a single cluster.
#'@param gamma is a threshold to ...
#'@param insigThs is a threshold to determine whether a magnitude of a feature coefficient is enough so that the feature is designated as a selected feature.
#'@param alpha is a significance level to determine whether a magnitude of a feature coefficient is enough so that the feature is designated as a selected feature.
#'@param minInvs is a minimum number of individuals for a cluster to be considered for inferring eta(C)cv, otherwise, eta(C)cv=0.
#'@param messageFlag is a flag. If it is true, the function shows the text regarding the progress of computing.
#'@param polyDegree is a degree of polynomial function that is used to fit the data.
#'If it is greater than 1, the polynomial formula is used in \code{lm()} instead of \code{"y=."}.
#'@param expFlag is an exponential flag to control the formula for data fitting.
#'If it is true, then the exp() formula is used in \code{lm()} instead of \code{"y=."}.
#'
#'@return This function returns \code{Copt}, \code{models}, \code{nNodes},  \code{invOptCls}, and \code{minR2cv}.
#'
#' \item{ \code{Copt[p,1]} }{ is equal to \code{k} implies a cluster that is a pth member of the maximal homogeneous partition is at kth layer and the cluster name in kth layer is \code{Copt[p,2]} }
#' \item{ \code{Copt[p,3]} }{ is "Model Information Reduction Ratio" \code{I({C},H0,Hlin)} of \code{p}th member of the maximal homogeneous partition: positive means the linear model is better than the null model.}
#' \item{ \code{Copt[p,4]} }{  is the squared correlation between predicted and real Y in CV step ( eta(C)cv ) of pth member of the maximal homogeneous partition. The greater \code{Copt[p,4]}, the higher homogeneous degree of this cluster.}
#' \item{ \code{ models[[k]][[j]]$clustInfoRecRatio} }{ is the "Cluster Information Reduction Ratio"  \code{I(Cj,Cjchildren,H)} between the \code{j}th cluster in \code{k}th layer
#' and its children clusters in \code{(k+1)}th layer: positive means current cluster is better than its children clusters.
#' Hence, we should keep this cluster at the member of maximal homogeneous partition instead of its children.}
#' \item{ \code{models[[j]][[k]]} }{ is a linear model of a cluster ID \code{k} at the layer \code{j}.
#'  The \code{models[[j]][[k]]$selFeatureSet} represents a set of selected-feature indices of the model where the feature index 1 is the intercept,
#'   and the feature index \code{d} is the (d-1)th variable \code{DataT$X[,d-1]}. }
#' \item{ \code{invOptCls[i,1]} }{ is the layer of optimal cluster of individual \code{i}. The optimal cluster of \code{i} is \code{invOptCls[i,2]}. }
#' \item{minR2cv}{ is the value of eta(C)cv from the cluster that has the lowest eta(C)cv.  }
#' \item{DataT}{is an updated \code{DataT} with the helper variables for plotting and printing results.}
#'
#'@examples
#'# Running FindMaxHomoOptimalPartitions using simulation data
#' DataT<-SimpleSimulation(100,type=1)
#' obj<-FindMaxHomoOptimalPartitions(DataT,gamma=0.05)
#'
#'@export
#'
FindMaxHomoOptimalPartitions<-function(DataT,gamma=0.05,insigThs=1e-8,alpha=0.05,minInvs=99,polyDegree = 1,expFlag=FALSE, messageFlag=FALSE)
{
  minR2cv<-Inf
  out<-dim(DataT$clsLayer)
  N<-out[[1]]
  nL<-out[[2]]
  Vc <-cbind(numeric(N),numeric(N)) # individual optimal clusters
  out<-linearModelTraining(DataT,insigThs,alpha,polyDegree=polyDegree,expFlag=expFlag,messageFlag=messageFlag) # training the models
  models<-out$models
  DataT<-out$DataT

  #options( warn = -1 )
  for(k in seq(1,nL)) # search each layer from top to buttom
  {

    nCls<-length(models[[k]]) # number of cluster in ith layer
    currLayerList<- unique(DataT$clsLayer[,k])
    #message("\014")
    for(j in seq(1,nCls)) # search each cluster in the layer k
    {

      inxFilterVec<-DataT$clsLayer[,k] == currLayerList[j] # selecting only members of jth cluster in ith layer
      H0coeff<-mean(DataT$Y[inxFilterVec],na.rm = TRUE)
      H0Residual<- DataT$Y[inxFilterVec]-H0coeff
      H1coeff<-models[[k]][[j]]$coefficients
      HlinResidual<-models[[k]][[j]]$residuals

      H0Residual[is.na(H0Residual)]<-0
      H1coeff[is.na(H1coeff)]<-0
      HlinResidual[is.na(HlinResidual)]<-0

      # positive modelInfoRecRatio  == H1 better than H0 or I({Cj,k },H0,Hlin)
      modelInfoRecRatio<- getModelInfoRecRatio(H0Residual,HlinResidual,H0coeff,H1coeff)$modelInfoRecRatio


      models[[k]][[j]]$modelInfoRecRatio<-modelInfoRecRatio # modelInfoRecRatio>0 mean better than H0
      if(min(Vc[inxFilterVec,1])==0 )
      {
        if(!is.null(models[[k]][[j]]$ChildrenCls))
        {
          HparentCoeff<-H1coeff #1 H1coeff
          HparentResidual<-HlinResidual #2 HlinResidual
          HparentCoeff[is.na(HparentCoeff)]<-0
          HparentResidual[is.na(HparentResidual)]<-0


          #======= cross-validation
          HchildrenCoeff<-list()
          HchildrenResidual<-list()

          clsInxVec<-DataT$clsLayer[inxFilterVec,k+1]
          #3 inxFilterVec #4 DataT, #5 models
          if(length(unique(clsInxVec))>1)
          {
            R2cv<-crossValEstFunc(DataT$X[inxFilterVec,],DataT$Y[inxFilterVec],clsInxVec)$r2 # eta_CV  the squared correlation between predicted and real Y in CV step.


            for(chCls in models[[k]][[j]]$ChildrenCls)
            {
              HchildrenCoeff<-append(HchildrenCoeff,t(as.list(models[[k+1]][[chCls]]$coefficients)) )
              HchildrenResidual<-append(HchildrenResidual,as.list(models[[k+1]][[chCls]]$residuals) )
            }
            HchildrenCoeff<-as.numeric(HchildrenCoeff)
            HchildrenCoeff[is.na(HchildrenCoeff)]<-0

            HchildrenResidual<-as.numeric(HchildrenResidual)
            HchildrenResidual[is.na(HchildrenResidual)]<-0
            #I(C′, {Cj,k },Hlin)
            clustInfoRecRatio<- getModelInfoRecRatio(HchildrenResidual,HparentResidual,HchildrenCoeff,HparentCoeff)$modelInfoRecRatio
            models[[k]][[j]]$clustInfoRecRatio<-clustInfoRecRatio

          }else
          {
            R2cv<-0
            clustInfoRecRatio<-0
          }
          models[[k]][[j]]$R2cv<-R2cv

          # end cross-validation

          #Append Cj,k to C∗ if I(C′, {Cj,k },Hlin) > 0 and eta(Cj,k )cv ≥ γ ;
          if(clustInfoRecRatio >0 && R2cv>=gamma )
          {
            minR2cv<-min(c(minR2cv,R2cv))
            Vc[inxFilterVec,1] <- k
            Vc[inxFilterVec,2] <- j
          }

        }else # No child
        {
          if(sum(inxFilterVec)>=minInvs)
          {
            R2cv<-crossVal10FoldEstFunc(DataT$X[inxFilterVec,],DataT$Y[inxFilterVec])$r2
          }
          else
          {
            R2cv<-0
          }
          minR2cv<-min(c(minR2cv,R2cv))
          models[[k]][[j]]$R2cv<-R2cv
          Vc[inxFilterVec,1] <- k
          Vc[inxFilterVec,2] <- j
        }
        if(messageFlag == TRUE)
          message(sprintf("Calculating Layer%d,Cls%d:modelInfoRecRatio %f, R2cv %f",k,j,modelInfoRecRatio,R2cv))
      }
    }
  }

  Copt<-unique(Vc) # a set of optimal clusters
  M<-dim(Copt)[1]
  Copt<-cbind(Copt,numeric(M),numeric(M))
  for(j in seq(1,M))
  {
    Copt[j,3] <- models[[Copt[j,1]]][[Copt[j,2]]]$modelInfoRecRatio
    Copt[j,4] <- models[[Copt[j,1]]][[Copt[j,2]]]$R2cv
  }


  return(list("Copt"=Copt,"models"=models,"minR2cv"=minR2cv,"invOptCls"=Vc,"DataT"=DataT) )
}

#'@title  PrintOptimalClustersResult
#'
#'@description
#' PrintOptimalClustersResult is a support function for printing the optimal clusters from FindMaxHomoOptimalPartitions function.
#'
#'@param resObj is an object list, which is the output of FindMaxHomoOptimalPartitions function
#'@param selFeature is a flag. If it is true, then the function shows the selected feature(s) of each optimal cluster.
#'
#'@return No return value, called for printing optimal clusters.
#'
#'@examples
#'# Running FindMaxHomoOptimalPartitions using simulation data
#' DataT<-SimpleSimulation(100,type=1)
#' obj<-FindMaxHomoOptimalPartitions(DataT,gamma=0.05)
#'# Printing the result
#' PrintOptimalClustersResult(obj)
#'
#'@export
#'
PrintOptimalClustersResult<-function(resObj, selFeature= FALSE)
{
  Copt<-resObj$Copt
  models<-resObj$models
  M<-dim(Copt)[1]
  #cat("\014")
  print("========== List of Optimal Clusters ==========")
  for(i in seq(1,M))
  {
    clsName<-models[[Copt[i,1]]][[Copt[i,2]]]$clsName
    clustInfoRecRatio<-models[[Copt[i,1]]][[Copt[i,2]]]$clustInfoRecRatio
    if(is.null(clustInfoRecRatio))
      clustInfoRecRatio<-NA
    print(sprintf("Layer%d,ClS-%s:clustInfoRecRatio=%.2f,modelInfoRecRatio=%.2f, eta(C)cv=%.2f",Copt[i,1],clsName,clustInfoRecRatio,Copt[i,3],Copt[i,4]) )
    if(selFeature==TRUE)
    {
      print("Selected features")
      print(models[[Copt[i,1]]][[Copt[i,2]]]$selFeatureSet)
    }
  }
  print(sprintf("min eta(C)cv:%f",resObj$minR2cv))
}
