

MaxProQQ<-
  function(InitialDesign, p_nom=0, temp0=0, nstarts=1, itermax=400, total_iter=1000000){
    
    DD<-as.matrix(InitialDesign)
    
    m<-nrow(DD)
    k<-ncol(DD)
    s=2
    
    lambda<-1/apply(DD,2,function(x) length(unique(x)))
    lambda[lambda==1/m]<-0
    
    localopm<-1
    
    predesign=c(t(DD))
    
    if(temp0==0){
      aaaa<-.C("MaxProMeasure", as.double(lambda), as.double(predesign), as.integer(p_nom), as.integer(k),
              as.integer(m), as.integer(s), measure=double(1), PACKAGE="MaxPro")
      delta<-aaaa$measure*0.01
      temp0=-delta/log(0.99)
    }
    
    t00<-Sys.time()
    aaa<-.C("MaxProQQ", as.double(lambda), as.integer(p_nom), as.integer(m),as.integer(k),as.integer(localopm),as.double(predesign),as.integer(nstarts),
            as.integer(itermax),as.integer(total_iter),design=double(m*k),measure=double(1), 
            as.double(temp0),ntotalI=integer(1),as.integer(s),PACKAGE="MaxPro")
    t01<-Sys.time()
    
    time_rec=t01-t00  
    Design=matrix(aaa$design,ncol=k,nrow=m,byrow=TRUE)
    
    val<-list(Design=Design,temp0=temp0,measure=(aaa$measure)^2,time_rec=time_rec,ntotal=aaa$ntotalI)
    
    return(val)
  }
