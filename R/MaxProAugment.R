 
MaxProAugment<-
  function(ExistDesign, CandDesign, nNew, p_disnum=0, l_disnum=NULL, p_nom=0, l_nom=NULL){
    
    CandDesign<-as.matrix(CandDesign)
    k<-ncol(CandDesign)
    
    if(is.data.frame(ExistDesign)) ExistDesign<-as.matrix(ExistDesign)
    if(!is.matrix(ExistDesign)){
      ExistDesign<-matrix(ExistDesign,ncol=k)
    } else{
      if(ncol(ExistDesign)!=ncol(CandDesign))
        print("Input error: Candidate design points and existing design points have different dimensions.")
    }
    
    p_cont<-k-p_disnum-p_nom
    if(p_cont<0) print("Input error: Summation of p_disnum and p_nom exceeds the total number of columns.")
    
    lambda<-rep(0,p_cont)
    
    if((p_disnum!=0)&(p_nom==0)){
      if(is.null(l_disnum)) 
        l_disnum<-apply(rbind(ExistDesign,CandDesign)[,(p_cont+1):k],2, function(x) length(unique(x)))
      if(length(l_disnum)!=p_disnum) print("Input error: Length of l_disnum does not match p_disnum.")
      lambda<-c(lambda,1/l_disnum)
    }
    
    if((p_disnum==0)&(p_nom!=0)){
      if(is.null(l_nom)) 
        l_nom<-apply(rbind(ExistDesign,CandDesign)[,(p_cont+1):k],2, function(x) length(unique(x)))
      if(length(l_nom)!=p_nom) print("Input error: Length of l_nom does not match p_nom.")
      lambda<-c(lambda,1/l_nom)
    }
    
    if((p_disnum!=0)&(p_nom!=0)){
      if(is.null(l_nom)) 
        l_nom<-apply(rbind(ExistDesign,CandDesign)[,(k-p_nom+1):k],2, function(x) length(unique(x)))
      if(is.null(l_disnum)) 
        l_disnum<-apply(rbind(ExistDesign,CandDesign)[,(p_cont+1):(p_cont+p_disnum)],2, function(x) length(unique(x)))
      
      if(length(l_nom)!=p_nom) print("Input error: Length of l_nom does not match p_nom.")
      if(length(l_disnum)!=p_disnum) print("Input error: Length of l_disnum does not match p_disnum.")
      
      lambda<-c(lambda,1/l_disnum,1/l_nom)
    }
    
    
    nExist<-nrow(ExistDesign)
    nCand<-nrow(CandDesign)
    
    if(nCand<nNew)
      print("Input error: The number of candidate design points is less than the number of new points to add.")
    nTotal<-nExist+nNew
    
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
    Design=matrix(aaa$design,ncol=k,nrow=nTotal,byrow=TRUE)
    
    if(aaa$measurecheck==0) 
      print("Note: Not enough candidate rows. For a continuous factor, any new added level must be distinct from the existing ones in the design. If repeated levels are needed, please specify the factor as discrete numeric.")
    
    val<-list(Design=Design,measure=(aaa$measure)^2,time_rec=time_rec)
    
    return(val)
  }

