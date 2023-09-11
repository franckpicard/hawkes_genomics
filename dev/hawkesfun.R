library(parallel)
library(bedr)
options(scipen=50,digits=5)

compute_hawkes_overlap = function(K,delta,kernel,processes,path_data,path_results,ofile_prefix){

  files.regions.in = paste(path_data,"/",processes,"_preprocessed.bed",sep="")  

  file.option     = paste(rep(" -f",length(files.regions.in)),files.regions.in,collapse="")
  file.out.hawkes = paste(path_results,"/",ofile_prefix,"_forward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep="")
  command         = paste(path_soft," ",file.option," -dump-intermediate-values -histogram ",K," ",delta," -kernel ",kernel," -lambda 1 > ",file.out.hawkes,sep="")  
  system(command)

  file.option     = paste(rep(" -b",length(files.regions.in)),files.regions.in,collapse="")
  file.out.hawkes = paste(path_results,"/",ofile_prefix,"_backward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep="")
  command         = paste(path_soft," ",file.option," -dump-intermediate-values -histogram ",K," ",delta," -kernel ",kernel," -lambda 1 > ",file.out.hawkes,sep="")  
  system(command)
  
}

getparam_overlap<-function(M,K,ofile, normalize=FALSE){  
  # les colonnes sont en m

  fileName = ofile 
  Res = scan(ofile, n = (M*M*K+M),what="character",skip=1)  

  Res = t(matrix(as.numeric(Res),nrow=M,ncol=M*K+1))
  nu  = Res[1,]
  Res = Res[-1,]
  
  ss  = seq(1,M*K,by=K)
  ee  = c((ss-1)[2:M],M*K)
  a   = array(data=NA,dim=c(M,M,K))
  for (m in 1:M){
    for (ell in 1:M){
      a[m,ell,] = Res[ss[ell]:ee[ell],m]  # proba d'observer des occurrences de m en plus après avoir observé des occurrences de ell
    }
  }
  if (normalize){
    for (m in 1:M){
      for (ell in 1:M){
        a[m,ell,] = a[m,ell,]/nu[m]
      }
    }
  }
  return(a)  
}


results_hawkes_overlap = function(K,delta,kernel,processes,path_data,path_results,ofile_prefix){  

  M  = length(processes)
  af = getparam_overlap(M,K,paste(path_results,"/",ofile_prefix,"_forward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep=""),normalize=FALSE)
  ab = getparam_overlap(M,K,paste(path_results,"/",ofile_prefix,"_backward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep=""),normalize=FALSE)
  a  = merge_fb(af,ab)

}



compute_hawkes_histogram = function(K,delta,kernel,processes,path_data,path_results,ofile_prefix){

  files.regions.in = paste(path_data,"/",processes,"_preprocessed.bed",sep="")  

  file.option     = paste(rep(" -f",length(files.regions.in)),files.regions.in,collapse="")
  file.out.hawkes = paste(path_results,"/",ofile_prefix,"_forward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep="")
  command         = paste(path_soft," ",file.option," -histogram ",K," ",delta," -kernel ",kernel," -lambda 1 > ",file.out.hawkes,sep="")  
  system(command)

  file.option     = paste(rep(" -b",length(files.regions.in)),files.regions.in,collapse="")
  file.out.hawkes = paste(path_results,"/",ofile_prefix,"_backward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep="")
  command         = paste(path_soft," ",file.option," -histogram ",K," ",delta," -kernel ",kernel," -lambda 1 > ",file.out.hawkes,sep="")  
  system(command)
  
}

results_hawkes_histogram = function(K,delta,kernel,processes,path_data,path_results,ofile_prefix){  

  af = getparam(paste(path_results,"/",ofile_prefix,"_forward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep=""),normalize=TRUE)
  ab = getparam(paste(path_results,"/",ofile_prefix,"_backward_K_",K,"_delta_",delta,"_kernel_",kernel,".txt",sep=""),normalize=TRUE)
  a  = merge_fb(af,ab)

}

df2bedr <- function(df,check.valid=FALSE){
  ll = unique(paste(df[,1],paste(df[,2],df[,3],sep="-"),sep=":"))
  bedr.sort.region(ll,check.valid=check.valid)
}

get_regions_from_bed <-function(dd,maxlag,check.valid=FALSE){

	dd.bed   = mclapply(dd,df2bedr,mc.cores=8)
	full.bed = unlist(dd.bed)
	full.bed = bedr.sort.region(full.bed,check.valid=check.valid)
	full     = convert2bed(full.bed,check.valid=check.valid)
	minocc   = 1
  chrlist  = unique(full$chr)

	RR = mclapply(chrlist,FUN=function(chr){
	  df      = full[full[,1]==chr,]
	  u       = maxlag
	  lag     = 1
	  chrsize = max(df[,3]) - min(df[,2]) + u 
	  ss = seq(min(df[,2]),chrsize,by=u/lag)
	  ee = seq(min(df[,2])+u-1,chrsize+u,by=u/lag)
	  hh = which(ss<(chrsize-u))
	  ss = ss[hh]
	  ee = ee[1:length(ss)]
	  ws = ss
	  we = ee
	  cov.bed      = paste(chr,paste(ws,we,sep="-"),sep=":")
	  fullchr.bed  = df2bedr(df)
	  inter        = bedr(input = list(a = cov.bed, b = fullchr.bed), method = "intersect", params = "-loj -header",check.valid = check.valid)
	  inter        = inter[inter$V4!=".",]
	  inter        = inter[!(inter$index %in% names(which(table(inter$index)<=minocc))),]	  
	  return(unique(inter$index))
	},mc.cores=8)

	RR = bedr.merge.region(unlist(RR),distance=1,check.valid=check.valid)
		
	return(RR)

}


getparam<-function(ofile, normalize=FALSE){  
  # les colonnes sont en m
  Res = read.table(ofile)
  nu  = unlist(Res[1,])
  Res = Res[-1,]
  M   = ncol(Res)  
  ss  = seq(1,M*K,by=K)
  ee  = c((ss-1)[2:M],M*K)
  a   = array(data=NA,dim=c(M,M,K))
  for (m in 1:M){
    for (ell in 1:M){
      a[m,ell,] = Res[ss[ell]:ee[ell],m]  # proba d'observer des occurrences de m en plus après avoir observé des occurrences de ell
    }
  }
  if (normalize){
    for (m in 1:M){
      for (ell in 1:M){
        a[m,ell,] = a[m,ell,]/nu[m]
      }
    }
  }
  return(a)  
}


merge_fb <-function(af,ab){  
  K    = dim(af)[3]
  M    = dim(af)[1]
  aa   = array(data=NA,dim=c(M,M,2*K))
  for (m in 1:M){
    for (ell in 1:M){
     aa[m,ell,] = c(ab[m,ell,K:1],af[m,ell,1:K])
    }
  }
  return(aa)
}

W <-function(x,eta){1*(abs(x)<=eta)/eta}

hw <- function(a,K,delta,eta_m,eta_ell){
  w_ell = sapply(1:(2*K*delta),FUN=function(x){W(x,eta_ell)})
  w_m   = sapply(1:(2*K*delta),FUN=function(x){W(x,eta_m)})
  hh    = rep(a,each=delta)
  convolve(w_ell,convolve(w_m,hh))
}

