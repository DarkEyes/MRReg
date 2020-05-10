#'@title  SimpleSimulation
#'
#'@description
#' SimpleSimulation is a support function for generating multiresolution datasets.
#'
#' All simulation types have three layers except the type 6 has four layers.
#'
#' The type-1 simulation has all individuals belong to the same homogeneous partition in the first layer.
#'
#' The type-2 simulation has four homogeneous partitions in a second layer. Each partition has its own models.
#'
#' The type-3 simulation has eight homogeneous partitions in a third layer. Each partition has its own models
#'
#' The type-4 simulation has one homogeneous partition in a second layer, four homogeneous partitions in a third layer,
#'  and eight homogeneous partitions in a fourth layer. Each partition has its own model.
#'
#' The type-5 simulation is similar to type-4 simulation but Y=h(X) is an exponential function.
#'
#' The type-6 simulation is similar to type-4 simulation but Y=h(X) is a polynomial function with \code{degree} parameter.
#'
#'@param indvN is a number of individuals per homogeneous partition.
#'@param type is a type of simulation dataset. There are four types.
#'@param degree is a degree parameter of a polynomial function for type-5 simulation
#'
#'@return The function returns a multiresolution dataset.
#'\item{\code{DataT$X[i,d]} }{ is a value of feature \code{d} of individual \code{i} }
#'\item{\code{DataT$Y[i]} }{ is value of target variable of individual \code{i} that
#'  we want to fit \code{DataT$Y ~ DataT$X} in linear model}
#'\item{\code{clsLayer[i,j]} }{ is a cluster ID of individual \code{i} at layer \code{j};
#' \code{clsLayer[i,1]} is the first layer that everyone typically belongs to a single cluster. }
#'\item{\code{DataT$TrueFeature[i]}}{ is equal to \code{d} if a true feature is \code{DataT$X[i,d-1]} that \code{DataT$Y[i]} is dependent with.
#' Note that \code{d = 1} is reserved for the intercept value in a linear model. }
#'
#'@examples
#'# Running SimpleSimulation to generate a dataset.
#' DataT<-SimpleSimulation(100,type=1)
#'
#'
#'@export
#'
SimpleSimulation<-function(indvN=10000,type=1,degree=2)
{
  if(type ==1)
    DataT<-clusterSimpleGenT1Func(indvN)
  else if(type ==2)
    DataT<-clusterSimpleGenT2Func(indvN)
  else if(type ==3)
    DataT<-clusterSimpleGenT3Func(indvN)
  else if(type ==4)
    DataT<-clusterSimpleGenT4Func(indvN)
  else if(type ==5)
    DataT<-clusterSimpleGenT5Func(indvN)
  else
    DataT<-clusterSimpleGenT6Func(indvN,degree=degree) # Type of simulation datasets

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

polySimFunc <- function(clsList,stdList,dumDim, degree =2) {
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
      Y[j]<- const1*(mat[j,i])^degree
    }
    k<-k+n
  }
  return(list("Y"=Y,"mat"=mat,"cls"=cls))
}

expSimFunc <- function(clsList,stdList,dumDim) {
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
      Y[j]<- const1*exp(mat[j,i])
    }
    k<-k+n
  }
  return(list("Y"=Y,"mat"=mat,"cls"=cls))
}

clusterSimpleGenT1Func <- function(indvN,nonLinearFlag= 0,degree=2) {
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
  if(nonLinearFlag == 0)
    out<- LinearSimFunc(clsList,stdList,dumDim)
  else if(nonLinearFlag ==1)
    out<- polySimFunc(clsList,stdList,dumDim,degree=degree)
  else
    out<- expSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"TrueFeature"=out$cls+1))
}

clusterSimpleGenT2Func <- function(indvN,nonLinearFlag= 0,degree=2) {
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
  if(nonLinearFlag == 0)
    out<- LinearSimFunc(clsList,stdList,dumDim)
  else if(nonLinearFlag ==1)
    out<- polySimFunc(clsList,stdList,dumDim,degree=degree)
  else
    out<- expSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"TrueFeature"=out$cls+1))
}

clusterSimpleGenT3Func <- function(indvN,nonLinearFlag= 0,degree=2) {
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
  if(nonLinearFlag == 0)
    out<- LinearSimFunc(clsList,stdList,dumDim)
  else if(nonLinearFlag ==1)
    out<- polySimFunc(clsList,stdList,dumDim,degree=degree)
  else
    out<- expSimFunc(clsList,stdList,dumDim)
  return(list("clsLayer"=clsLayer,"Y"=out$Y,"X"=out$mat,"TrueFeature"=out$cls+1))
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
  TrueFeature<-rbind(Data1$TrueFeature,Data2$TrueFeature,Data3$TrueFeature)
  return(list("clsLayer"=clsLayer,"Y"=Y,"X"=X,"TrueFeature"=TrueFeature))
}


clusterSimpleGenT5Func <- function(indvN) {
  Data1<-clusterSimpleGenT1Func(indvN,nonLinearFlag= 2)
  Data2<-clusterSimpleGenT2Func(indvN,nonLinearFlag= 2)
  Data3<-clusterSimpleGenT3Func(indvN,nonLinearFlag= 2)

  clsLayer<-rbind(Data1$clsLayer,Data2$clsLayer+10,Data3$clsLayer+20)
  rN<-dim(clsLayer)[1]
  clsLayer<-cbind(numeric(rN)+1,clsLayer)
  Y<-rbind(Data1$Y,Data2$Y,Data3$Y)
  X<-rbind(Data1$X,Data2$X,Data3$X)
  TrueFeature<-rbind(Data1$TrueFeature,Data2$TrueFeature,Data3$TrueFeature)
  return(list("clsLayer"=clsLayer,"Y"=Y,"X"=X,"TrueFeature"=TrueFeature))
}

clusterSimpleGenT6Func <- function(indvN,degree=2) {
  Data1<-clusterSimpleGenT1Func(indvN,nonLinearFlag= 1,degree=degree)
  Data2<-clusterSimpleGenT2Func(indvN,nonLinearFlag= 1,degree=degree)
  Data3<-clusterSimpleGenT3Func(indvN,nonLinearFlag= 1,degree=degree)

  clsLayer<-rbind(Data1$clsLayer,Data2$clsLayer+10,Data3$clsLayer+20)
  rN<-dim(clsLayer)[1]
  clsLayer<-cbind(numeric(rN)+1,clsLayer)
  Y<-rbind(Data1$Y,Data2$Y,Data3$Y)
  X<-rbind(Data1$X,Data2$X,Data3$X)
  TrueFeature<-rbind(Data1$TrueFeature,Data2$TrueFeature,Data3$TrueFeature)
  return(list("clsLayer"=clsLayer,"Y"=Y,"X"=X,"TrueFeature"=TrueFeature))
}
