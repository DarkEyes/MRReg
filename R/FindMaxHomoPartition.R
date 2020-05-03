#'
#'@export
#'
FindMaxHomoPartition<-function(DataT,gamma=0.05)
{
  minR2cv<-Inf
  out<-dim(DataT$clsLayer)
  N<-out[[1]]
  nL<-out[[2]]
  Vc <-cbind(numeric(N),numeric(N)) # individual optimal clusters
  out<-linearModelTraining(DataT)
  models<-out$models
  DataT<-out$DataT

  #options( warn = -1 )
  for(k in seq(1,nL))
  {

    nCls<-length(models[[k]]) # number of cluster in ith layer
    currLayerList<- unique(DataT$clsLayer[,k])
    cat("\014")
    for(j in seq(1,nCls))
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
            R2cv<-crossValEstFunc(DataT$X[inxFilterVec,],DataT$Y[inxFilterVec],clsInxVec)$r2

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
          #         print(sprintf("R2cv:%f",R2cv))
          models[[k]][[j]]$R2cv<-R2cv

          # end cross-validation

          #Append Cj,k to C∗ if I(C′, {Cj,k },Hlin) > 0 and η(Cj,k )cv ≥ γ ;
          if(clustInfoRecRatio >0 && R2cv>=gamma )
          {
            minR2cv<-min(c(minR2cv,R2cv))
            Vc[inxFilterVec,1] <- k
            Vc[inxFilterVec,2] <- j
          }

        }else # No child
        {
          if(sum(inxFilterVec)>100)
          {
            R2cv<-crossVal10FoldEstFunc(DataT$X[inxFilterVec,],DataT$Y[inxFilterVec])$r2
          }
          else
          {
            R2cv<-0
          }
          #         print(sprintf("R2cv:%f",R2cv))
          minR2cv<-min(c(minR2cv,R2cv))
          models[[k]][[j]]$R2cv<-R2cv
          Vc[inxFilterVec,1] <- k
          Vc[inxFilterVec,2] <- j
        }
        print(sprintf("Calculating Layer%d,Cls%d:modelInfoRecRatio %f, R2cv %f",k,j,modelInfoRecRatio,R2cv))
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


  return(list("Copt"=Copt,"models"=models,"minR2cv"=minR2cv,"DataT"=DataT,"Vc"=Vc) )
}

#'
#'@export
#'
PrintPartitionResult<-function(resObj)
{
  Copt<-resObj$Copt
  models<-resObj$models
  M<-dim(Copt)[1]
  cat("\014")
  print("========== List of Optimal Clusters ==========")
  for(i in seq(1,M))
  {
    clsName<-models[[Copt[i,1]]][[Copt[i,2]]]$clsName
    print(sprintf("Layer%d,ClS-%s:modelInfoRecRatio=%.2f, R2cv=%.2f",Copt[i,1],clsName,Copt[i,3],Copt[i,4]) )
  }
  print(sprintf("minR2cv:%f",resObj$minR2cv))
  #print(Copt)
  #options( warn = 1 )
}
