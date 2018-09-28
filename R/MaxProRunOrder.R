 
MaxProRunOrder<-function(Design,p_nom=0,initial_row=1){
  
  ExistDesign<-as.matrix(Design)[initial_row,]
  CandDesign<-as.matrix(Design)[-initial_row,]
  
  k<-ncol(Design)
  nTotal<-nrow(Design)
  
  lambda<-1/apply(Design,2,function(x) length(unique(x)))
  lambda[lambda==1/nTotal]<-0
  
  nExist<-length(initial_row)
  nCand<-nTotal-nExist
  nNew<-nCand
  
  ExistDVec<-c(t(ExistDesign))
  CandDVec<-c(t(CandDesign))
  
  s<-2
  
  t00<-Sys.time()
  aaa<-.C("MaxProAugment", as.double(lambda), as.integer(p_nom), as.integer(k),
          as.integer(nExist),as.integer(nNew),as.integer(nCand),
          as.double(ExistDVec),as.double(CandDVec),design=double(nTotal*k),
          as.integer(s),measure=double(1), measurecheck=integer(1), PACKAGE="MaxPro")
  t01<-Sys.time()
  
  time_rec<-t01-t00 
  
  Design=cbind(1:nTotal,matrix(aaa$design,ncol=k,nrow=nTotal,byrow=TRUE))
  
  val<-list(Design=Design,measure=(aaa$measure)^2,time_rec=time_rec)
  
  return(val)
  
}

