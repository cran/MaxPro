
CandPoints<-function(N,p_cont,l_disnum=NULL,l_nom=NULL){
  
  p1<-p_cont
  p2<-length(l_disnum) 
  p3<-length(l_nom)
  
  ContD<-(apply(matrix(rep(1:N,p1),ncol=p1),2,sample)-.5)/N
  
  DisNumD<-matrix(0,nrow=N,ncol=p2)
  if(p2>0){
    for(ww in 1:p2){
      DisNumD[,ww]<-sample(rep(((1:l_disnum[ww])-1)/(l_disnum[ww]-1),ceiling(N/l_disnum[ww])),N)
    }
  }
  
  NomD<-matrix(0,nrow=N,ncol=p3)
  if(p3>0){
    for(ww in 1:p3){
      NomD[,ww]<-sample(rep(1:l_nom[ww],ceiling(N/l_nom[ww])),N)
    }
  }
  
  val<-cbind(ContD,DisNumD,NomD)
  return(val)
}
