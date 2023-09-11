rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)


mytheme<- function(){
  theme(
    legend.title      = element_blank(),
    legend.background = element_blank(),
    legend.text       = element_text(size=10),
    legend.key.size   = unit(0.5,"cm"),
    axis.text.x       = element_text(size = 8,angle = 0,hjust=0.5,vjust=0.5),    
    axis.title.x      = element_text(size=8),
    axis.title.y      = element_blank(),
    axis.text.y       = element_text(size=8,hjust=0.5),
    plot.margin       = unit(rep(0.2,4), "lines"),
    panel.background  = element_rect(fill=NA,color="black"),
    #strip.background  = element_rect(fill="white",color="black"),
    strip.background = element_blank(),
    panel.border      = element_rect(fill = NA, color = NA),
    strip.text        = element_text(size=8),
    plot.title        = element_text(hjust=0.5),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    
  )
}

path_code =""
source(paste0(code_path,"hawkesfun.R"))

path_data    = "/home/picard/Dropbox/Hawkes-genomics/cecile/data"
path_results = "/home/picard/Dropbox/Hawkes-genomics/cecile/results"
path_soft    = "/home/picard/projets/hawkes/hawkes/./hawkes"

eta    = c(1400,1,560)
K      = 10
delta  = 10000
kernel = "heterogeneous_interval"

x  = c(1:(K*delta))
x  = c(-x[length(x):1],x)

compute_hawkes_histogram(K,delta,kernel,
	processes    = c("enh","prom","sin3"),
	path_data    = path_data,
	path_results = path_results,
	ofile_prefix = "hawkes_histogram")

a =results_hawkes_histogram(K,delta,kernel,
	processes    = c("enh","prom","sin3"),
	path_data    = path_data,
	path_results = path_results,
	ofile_prefix = "hawkes_histogram")

pp = c("enh","prom","sin3")
ww = mclapply(1:3,FUN=function(i){
  mclapply(1:3,FUN=function(j){
    hh = hw(a[i,j,]/sqrt(delta),K,delta,eta[i],eta[j])
    data.frame(h = hh[seq(1,length(x),by=100)], x = x[seq(1,length(x),by=100)], pp1 = pp[i], pp2 = pp[j])
  })
})
ww = Reduce("rbind",Reduce("rbind",ww))


gg = list()
h  = 0 
for (Fi in  c("enh","prom","sin3")){
  for (Fj in  c("enh","prom","sin3")){
    h       = h+1    
    gg[[h]] = ggplot(ww[ (ww$pp1 ==Fi) & (ww$pp2 == Fj),],aes(x=x/1000,y=h,color=method)) + 
                    geom_vline(xintercept=0, color="gray", size=0.5) + geom_line(color="black") + 
                    xlab("position kb") + ylab("") + xlab("")+
                    mytheme()
  }
}
gg[[1]] = gg[[1]] + ggtitle(expression("enhancer" %->% "enhancer | prom, sin3"))
gg[[2]] = gg[[2]] + ggtitle(expression("prom"  %->% "enhancer | sin3"))
gg[[3]] = gg[[3]] + ggtitle(expression("sin3" %->% "enhancer | prom")) + geom_hline(yintercept=0, color="gray", size=0.5)
gg[[4]] = gg[[4]] + ggtitle(expression("enhancer" %->% "prom | sin3"))
gg[[5]] = gg[[5]] + ggtitle(expression("prom" %->%  "prom | enhancer,sin3"))
gg[[6]] = gg[[6]] + ggtitle(expression("sin3" %->% "prom | enhancer")) + geom_hline(yintercept=0, color="gray", size=0.5)
gg[[7]] = gg[[7]] + ggtitle(expression("enhancer" %->% "sin3 | prom")) + geom_hline(yintercept=0, color="gray", size=0.5)
gg[[8]] = gg[[8]] + ggtitle(expression("prom" %->%  "sin3 | enhancer")) + geom_hline(yintercept=0, color="gray", size=0.5)
gg[[9]] = gg[[9]] + ggtitle(expression("sin3" %->% "sin3 | enhancer,prom"))+ geom_hline(yintercept=0, color="gray", size=0.5)
gg = grid.arrange(gg[[1]],gg[[2]],gg[[3]],gg[[4]],gg[[5]],gg[[6]],gg[[7]],gg[[8]],gg[[9]])


