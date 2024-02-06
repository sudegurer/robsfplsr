daprpr <- function(Data, type = c("classical", "robust"))
{
  
  type <- match.arg(type)
  
  if(type == "classical"){
    Data.center <- apply(Data, 2, mean)
  }else if(type == "robust"){
    if(dim(Data)[2] == 1){
      Data.center <- median(Data)
    }else{
      Data.center <- l1median(Data)
    }
    
  }
  
  Data.scaled <- (Data - matrix(Data.center, nrow=dim(Data)[1], ncol=dim(Data)[2], byrow=TRUE))
  
  if(type == "classical"){
    Data.scale <- apply(Data, 2, sd)
  }else if(type == "robust"){
    Data.scale <- apply(Data, 2, qn)
  }
  
  Data.scaled <- Data.scaled/matrix(Data.scale, nrow=dim(Data)[1], ncol=dim(Data)[2], byrow=TRUE)
  
  attr(Data.scaled,"Center") <- Data.center
  attr(Data.scaled,"Scale") <- Data.scale
  return(Data.scaled)
}
