#'@title  plotOptimalClustersTree
#'
#'@description
#' plotOptimalClustersTree is a support function for plotting the hierarchical tree of optimal clusters
#' from FindMaxHomoOptimalPartitions function.
#'
#' The red node(s) are the optimal homogeneous clusters while the gray nodes are non-optimal clusters.
#'
#'@param resObj is an object list, which is the output of FindMaxHomoOptimalPartitions function
#'
#'@return No return value, called for plotting the hierarchical tree of optimal clusters.
#'
#'@examples
#'# Running FindMaxHomoOptimalPartitions using simulation data
#' DataT<-SimpleSimulation(100,type=1)
#' obj<-FindMaxHomoOptimalPartitions(DataT,gamma=0.05)
#'# Plotting the result
#' plotOptimalClustersTree(obj)
#'
#'@importFrom graphics plot
#'@import igraph
#'@export
#'
plotOptimalClustersTree<-function(resObj)
{

  adjMat<- matrix(0,resObj$DataT$nNodes,resObj$DataT$nNodes)
  nameList<-list()
  graphOPTflag<-logical()


  mainCOPT<-resObj$Copt
  nL<-length(resObj$models)
  k<-1
  optCount<-0

  for(i in seq(1,nL)) # ith Layer
  {

    # =============
    flag1<-FALSE
    if( is.element(i,mainCOPT[,1] ) ) # this layer contains opt?
    {
      inxVec<-mainCOPT[,1] == i
      flag1<-TRUE

    }
    currLayer<- resObj$DataT$clsLayer[,i]
    nCls<-length(unique(currLayer) )

    for(j in seq(1,nCls)) # Nodes in ith layer
    {
      if(i>1) # build a graph
      {
        cID<-resObj$models[[i]][[j]]$ID
        parentCls<-resObj$models[[i]][[j]]$ParentCls
        pID<-resObj$models[[i-1]][[parentCls]]$ID
        adjMat[cID,pID]<- 1
        i
      }

      nameList[k]<-sprintf("L%dC%g",i,resObj$DataT$clsNameMappingTable[[i]][[j]])

      if( flag1 )
      {

        if(is.element(j,mainCOPT[inxVec,2])) # if found opt
        {
          graphOPTflag[k]<-TRUE

          optCount<-optCount+1
        }
        else
          graphOPTflag[k]<-FALSE
      }
      else
        graphOPTflag[k]<-FALSE
      k<-k+1
    }
  }


  g1 <- graph_from_adjacency_matrix( adjMat   ) %>%
    set_vertex_attr("label", value = nameList)


  V(g1)$color <- ifelse(graphOPTflag == TRUE, "red", "gray")

  #set_vertex_attr(g1,"label", value = nameList)
  V(g1)$label.cex = 0.8
  return(plot(g1, layout =  layout.auto,edge.arrow.size=0.25,vertex.label.color = "black"))
}
