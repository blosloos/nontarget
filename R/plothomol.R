plothomol<-
function(
	homol,
	xlim=FALSE,
	ylim=FALSE,
	plotlegend=TRUE
){

    ############################################################################
    # check inputs #############################################################
    if(xlim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
    if(ylim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
    if(length(homol[[5]])<1){stop("no homologue series found!")}
    ############################################################################
    
    ############################################################################
    # plot #####################################################################
    sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
    plot.new();
    if(xlim[1]!=FALSE & ylim[1]==FALSE){plot.window(xlim=xlim,ylim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));}
    if(ylim[1]!=FALSE & xlim[1]==FALSE){
      if(plotlegend==TRUE){
        plot.window(xlim=c(min(homol[[1]][,3]),max(homol[[1]][,3])*1.2),ylim=ylim);
      }else{
        plot.window(xlim=c(min(homol[[1]][,3]),max(homol[[1]][,3])),ylim=ylim);
      }
    }
    if(xlim[1]==FALSE & ylim[1]==FALSE){
      if(plotlegend==TRUE){
        plot.window(xlim=c(min(homol[[1]][,3]),max(homol[[1]][,3])*1.2),ylim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));
      }else{
        plot.window(xlim=c(min(homol[[1]][,3]),max(homol[[1]][,3])),ylim=c(min(homol[[1]][,1]),max(homol[[1]][,1])));    
      }
    }
    if(xlim[1]!=FALSE & ylim[1]!=FALSE){plot.window(xlim=xlim,ylim=ylim);}
    box();axis(1);axis(2);
    title(xlab="Retention time",ylab="m/z");
    points(homol[[1]][,3],homol[[1]][,1],cex=0.3,pch=19,col="lightgrey");
    ############################################################################
    
    ############################################################################
    this<-round(homol[[3]][,3],digits=2);
    that<-levels(as.factor(this));
    colo<-rainbow(length(that))
    for(i in 1:length(homol[[5]])){
	
	
      for(j in 2:length(homol[[5]][[i]])){
       for(k in 1:length(homol[[5]][[i]][j-1])){
        for(m in 1:length(homol[[5]][[i]][j])){
          lines (
            c(homol[[1]][homol[[5]][[i]][[j-1]][k],3],homol[[1]][homol[[5]][[i]][[j]][m],3]),
            c(homol[[1]][homol[[5]][[i]][[j-1]][k],1],homol[[1]][homol[[5]][[i]][[j]][m],1]),
            col=colo[that==this[i]],lwd=1.8
                );
          points(
             homol[[1]][homol[[5]][[i]][[j-1]][k],3],
             homol[[1]][homol[[5]][[i]][[j-1]][k],1],
             col=colo[that==this[i]],pch=19,cex=0.5
          );
          points(
             homol[[1]][homol[[5]][[i]][[j]][m],3],
             homol[[1]][homol[[5]][[i]][[j]][m],1],
             col=colo[that==this[i]],pch=19,cex=0.5
          );
        } 
       }
      }
	  #i<-(i+1)
	  
	  
    }
    # add a legend
    if(plotlegend==TRUE){
      plot.window(xlim=c(0,1),ylim=c(min(as.numeric(that)),max(as.numeric(that))));
      lines(c(0.95,0.95),c(min(as.numeric(that)),max(as.numeric(that))),col="lightgrey",lwd=6)
      it<-2
      for(i in 1:length(that)){
          points(0.95,as.numeric(that[i]),pch=19,col=colo[i])
          text(0.95,as.numeric(that[i]),labels=that[i],col=colo[i],cex=0.65,pos=it)
          if(it==2){it<-4}else{it<-2}
      }

    }
    ############################################################################
    
}
