

MaxProMeasure<-function(Design, p_nom=0){
  Design<-as.matrix(Design)
  p<-ncol(Design)
  n<-nrow(Design)
  lambda<-1/apply(Design,2,function(x) length(unique(x)))
  lambda[lambda==1/n]<-0
  DVec<-c(t(Design))
  s<-2
  aaa<-.C("MaxProMeasure", as.double(lambda), as.double(DVec), as.integer(p_nom), as.integer(p),
          as.integer(n), as.integer(s), measure=double(1), PACKAGE="MaxPro")
  val<-(aaa$measure)^2
  return(val)
} 
