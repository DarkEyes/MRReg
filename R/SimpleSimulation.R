#'
#'@export
#'
SimpleSimulation<-function(indvN=10000,type=1)
{
  if(type ==1)
    DataT<-clusterSimpleGenT1Func(indvN)
  else if(type ==2)
    DataT<-clusterSimpleGenT2Func(indvN)
  else if(type ==3)
    DataT<-clusterSimpleGenT3Func(indvN)
  else
    DataT<-clusterSimpleGenT4Func(indvN) # Type of simulation datasets
  return(DataT)
}

LinearSimFunc <- function(clsList,stdList,dumDim) {
  totalDim <-length(clsList) + dumDim
  totalRow <- sum(clsList)
  mat<- matrix( rnorm(totalRow*totalDim,mean=0,sd=1), totalRow, totalDim)
  Y <-matrix(0, totalRow, 1)
  cls<-matrix(0, totalRow, 1)
  k<-1
  for(i in seq(1, length(clsList), by = 1) )
  {
    n<-clsList[i]
    const1<-runif(1, min=2, max=10)
    for( j in seq(k,k+n-1) )
    {
      cls[j] <-i
      Y[j]<- const1*mat[j,i]
    }
    k<-k+n
  }
  return(list("Y"=Y,"mat"=mat,"cls"=cls))
}


clusterSimpleGenT1Func <- function(indvN) {
  N<-indvN*8
  clsLayer<-matrix(1, N, 3)
  nL2<-indvN*2
  for (i in seq(1,4) ) # fill 2nd layer
  {
    st<-(i-1)*nL2+1
    fn<-i*nL2
    clsLayer[ st:fn, 2]<-i;

  }
  for (i in seq(1,8) ) # fill 3nd layer
  {
    st<-(i-1)*indvN+1
    fn<-i*indvN
    clsLayer[ st:fn, 3]<-i;

  }
  clsList<-cbind(N)
  stdList<-c(1)
  dumDim<-19
  out<- LinearSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"cls"=out$cls))
}

clusterSimpleGenT2Func <- function(indvN) {
  N<-indvN*8
  clsLayer<-matrix(1, N, 3)
  nL2<-indvN*2
  for (i in seq(1,4) ) # fill 2nd layer
  {
    st<-(i-1)*nL2+1
    fn<-i*nL2
    clsLayer[ st:fn, 2]<-i;

  }
  for (i in seq(1,8) ) # fill 3nd layer
  {
    st<-(i-1)*indvN+1
    fn<-i*indvN
    clsLayer[ st:fn, 3]<-i;

  }
  clsList<-cbind(nL2,nL2,nL2,nL2)
  stdList<-c(1,1,1,1)
  dumDim<-16
  out<- LinearSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"cls"=out$cls))
}

clusterSimpleGenT3Func <- function(indvN) {
  N<-indvN*8
  clsLayer<-matrix(1, N, 3)
  nL2<-indvN*2
  for (i in seq(1,4) ) # fill 2nd layer
  {
    st<-(i-1)*nL2+1
    fn<-i*nL2
    clsLayer[ st:fn, 2]<-i;

  }
  for (i in seq(1,8) ) # fill 3nd layer
  {
    st<-(i-1)*indvN+1
    fn<-i*indvN
    clsLayer[ st:fn, 3]<-i;

  }
  clsList<-cbind(indvN,indvN,indvN,indvN,indvN,indvN,indvN,indvN)
  stdList<-c(1,1,1,1,1,1,1,1)
  dumDim<-12
  out<- LinearSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"cls"=out$cls))
}

clusterSimpleGenT4Func <- function(indvN) {
  Data1<-clusterSimpleGenT1Func(indvN)
  Data2<-clusterSimpleGenT2Func(indvN)
  Data3<-clusterSimpleGenT3Func(indvN)

  clsLayer<-rbind(Data1$clsLayer,Data2$clsLayer+10,Data3$clsLayer+20)
  rN<-dim(clsLayer)[1]
  clsLayer<-cbind(numeric(rN)+1,clsLayer)
  Y<-rbind(Data1$Y,Data2$Y,Data3$Y)
  X<-rbind(Data1$X,Data2$X,Data3$X)
  cls<-rbind(Data1$cls,Data2$cls,Data3$cls)
  return(list("clsLayer"=clsLayer,"Y"=Y,"X"=X,"cls"=cls))
}
