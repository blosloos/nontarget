homol.search <-
function(
	peaklist,
	isotopes,
	elements=c("C","H","O"),
	charge=c(1,2),
	use_C=TRUE,
	minmz=5,
	maxmz=60,
	minrt=0.5,
	maxrt=2,
	ppm=TRUE,
	mztol=3.5,
    rttol=0.5,
	minlength=4,
	mzfilter=FALSE,
	vec_size=1E6
){
  
    ##########################################################################
    # (0) Issue warnings: check arguments ####################################
    if(ppm==TRUE & mztol>10){cat("Too big mztol?")};
    if(minlength<3){stop("invalid minlen argument. set minlen >= 3")};
    if(mztol<=0){warning("mztol should be >0!")};
	if(length(elements)<1){stop("specify elements")}
    charge<-abs(charge);
	if(!is.numeric(mzfilter) & mzfilter[1]!="FALSE"){stop("mzfilter must either be a numeric vector or set to FALSE")}
	if(any(is.na(match(elements,isotopes[,1])))){ stop("unknown elements") }
	if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
	if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
	if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
	##########################################################################
    # (1) retrieve feasible mass differences & all combinations thereof ######
	# (1.1) upper & lower mass defect / mass bound ###########################
	# only requires monoisotopic mass ########################################
	# relation to mass makes charges irrelevant ##############################
	# - not used if mzfilter is set ##########################################
	##########################################################################
	if(ppm==TRUE){
		delmz<-(mztol*max(peaklist[,1])/1e6);
	}else{
		delmz<-mztol;
	}
	if(mzfilter[1]==FALSE){
		delmass<-c();
		c_ratio<-c();
		for(i in 1:length(elements)){
			isos<-isotopes[isotopes[,1]==elements[i],]
			delmass<-c(delmass,
				isos[isos[,3]==min(isos[,3]),3]
			);
			c_ratio<-c(c_ratio,unique(isos[,5]));
		}	
		if( (min(delmass)/max(charge)) > (minmz/max(charge)) ){
			minmz<-min(delmass)
		}
		c_ratio[c_ratio==0]<-Inf # well...
		if(use_C){
			delmass<-( c( delmass-round(delmass,digits=0) ) / (delmass+(1/c_ratio*12)) );		
		}else{
			delmass<-( c( delmass-round(delmass,digits=0) ) / delmass );
		}
		maxup<-max(delmass);
		maxdown<-abs(min(delmass));
		# (1.2) combinatorial restrictions: nominal mass windows #################
		lowbound<-c(0)
		highbound<-c(0)
		atnominal<-min(1/charge)
		while(max(lowbound)<maxmz){
			lowbound<-c(lowbound,
				atnominal-(maxdown*atnominal)
			)
			highbound<-c(highbound,
				atnominal+(maxup*atnominal)
			)
			atnominal<-(atnominal+min(1/charge))
		}
		lowbound<-lowbound[-length(lowbound)]
		lowbound<-lowbound-delmz
		highbound<-highbound[-length(highbound)]
		highbound<-highbound+delmz
		highbound<-highbound[lowbound>=minmz]
		lowbound<-lowbound[lowbound>=minmz]	
		# (1.3) correct for bound overlaps #######################################
		if(length(lowbound)>1){
			new_lowbound<-lowbound[1]
			new_highbound<-highbound[1]
			for(i in 2:length(highbound)){
				if( new_highbound[length(new_highbound)] >= lowbound[i]){
					new_highbound[length(new_highbound)]<-highbound[i]
				}else{
					new_lowbound<-c(new_lowbound,lowbound[i])
					new_highbound<-c(new_highbound,highbound[i])
				}
			}
			lowbound<-new_lowbound
			highbound<-new_highbound
			# cbind(lowbound,highbound)
		}
	}else{	
		mzfilter<-mzfilter[order(mzfilter,decreasing=FALSE)]
		lowbound<-(mzfilter-delmz)	
		highbound<-(mzfilter+delmz)
	}	
    ##########################################################################
	# run calculation ########################################################
	##########################################################################
	ord<-order(peaklist[,1],decreasing=FALSE) # supply input by increas. m/z #
	peaklist2<-peaklist[ord,]
	if(ppm=="TRUE"){
		ppm2<-1
	}else{
		ppm2<-2
	}
	cat("Detecting:\n")
	pBar <- txtProgressBar( min = 0, max = length(highbound), style = 3 );
	inter<-as.numeric(interactive())
	relat <- .Call("homol_1",
		as.integer(vec_size),          	# storage size
		as.numeric(peaklist2[,1]), 		# mass
		as.numeric(peaklist2[,3]), 		# RT 
		as.integer(ord),		  		# order
		as.numeric(lowbound),	  		# lower bound mass window	
		as.numeric(highbound),	  		# higher bound mass window	 	
		as.numeric(minrt),		  		# minimum RT change
		as.numeric(maxrt),        		# maximum RT change
		as.numeric(mztol),		  		# mass tolerance
		as.integer(ppm2),	   	  		# ppm==1 -> TRUE
		as.numeric(rttol),		  		# RT tolerance		
		as.integer(minlength),    		# minimum length of HS
		as.integer(inter),
		pBar,
		PACKAGE="nontarget"
	)
	close(pBar)
    if(length(relat) == 0){
		stop("no series detected");
    }
	cat("\n done.\n")
	relat<-relat[order(relat[,5],peaklist[relat[,1],1]),] # sort by increasing mass
	names(relat)<-c("from_ID","to_ID","dmz","dRT","HS_ID")
    ##########################################################################
	
	##########################################################################
    # Generate outputs #######################################################
    ##########################################################################
    # generate peaklist with links & #####################################
    # generate component list with relevant m/z & RT increments ##########
    group1<-c();  # group ID
    group2<-c();  # peak IDs
    group3<-c();  # mzincr
    group4<-c();  # retincr
    group5<-c();  # retmin
    group6<-c();  # retmax
	len<-length(peaklist[,1])
    getit1<-rep("0",len);      	# (1) group ID
    getit2<-rep("0",len);      	# (2) level
    getit3<-rep("0",len);   	# (3) to which peak
    getit4<-rep("none",len);   	# (4) mean dmass
    getit5<-rep("none",len);    # (5) mean RT
	peakID<-list(0)
	if(length(relat[,1])>0){
		HS<-unique(relat[,5])
		for(i in 1:length(HS)){
			peakID[[i]]<-numeric(0)
			atgroup<-relat[relat[,5]==HS[i],];
			# entry to homol[[1]] = extended peaklist ########################	
			# ...those pointing
			for(j in 1:length(atgroup[,1])){
				getit1[atgroup[j,1]]<-paste(getit1[atgroup[j,1]],"/",as.character(HS[i]),sep="")
				getit2[atgroup[j,1]]<-paste(getit2[atgroup[j,1]],"/",as.character(j),sep="")
				getit3[atgroup[j,1]]<-paste(getit3[atgroup[j,1]],"/",as.character(atgroup[j,2]),sep="")
				getit4[atgroup[j,1]]<-paste(getit4[atgroup[j,1]],"/",as.character(round(mean(atgroup[,3]),digits=4)),sep="")
				getit5[atgroup[j,1]]<-paste(getit5[atgroup[j,1]],"/",as.character(round(mean(atgroup[,4]),digits=4)),sep="")
			}
			# ... those only pointed at
			for(p in 1:length(atgroup[,2])){
				if( !any(atgroup[,1]==atgroup[p,2]) ){
					getit1[atgroup[p,2]]<-paste(getit1[atgroup[p,2]],"/",as.character(HS[i]),sep="")
					getit2[atgroup[p,2]]<-paste(getit2[atgroup[p,2]],"/",as.character(j+p),sep="")
					getit4[atgroup[p,2]]<-paste(getit4[atgroup[p,2]],"/",as.character(round(mean(atgroup[,3]),digits=4)),sep="")
					getit5[atgroup[p,2]]<-paste(getit5[atgroup[p,2]],"/",as.character(round(mean(atgroup[,4]),digits=4)),sep="")
				}
				
			}			
			# entry to homol[[3]] = homologue series #########################
			# entry to homol[[5]] = peakIDs in list per level ################
			listit<-list(0)
			group1<-c(group1,HS[i]);
			atit<-unique(c(atgroup[,1],atgroup[,2]));
			atthat<-as.character(atit[1]);
			listit[[1]]<-atit[1];
			for( j in 2:length(atit)){
				atthat<-paste(atthat,",",as.character(atit[j]),sep="");
				listit[[j]]<-atit[j];
			}
			group2<-c(group2,atthat);
			group3<-c(group3,mean(atgroup[,3]));
			group4<-c(group4,mean(atgroup[,4]));
			group5<-c(group5,min(atgroup[,4]));
			group6<-c(group6,max(atgroup[,4]));
			peakID[[HS[i]]]<-listit
		}
		for(i in 1:length(getit1)){
			if(getit1[i]!="0"){
			  getit1[i]<-substr(getit1[i],3,nchar(getit1[i]))
			  getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))
			  getit4[i]<-substr(getit4[i],6,nchar(getit4[i]))
			  getit5[i]<-substr(getit5[i],6,nchar(getit5[i]))
			}
			if(getit3[i]!="0"){
			  getit3[i]<-substr(getit3[i],3,nchar(getit3[i]))
			}
		}
		grouped_samples<-data.frame(peaklist,seq(1,len,1),getit1,getit2,getit3,getit4,getit5);
		names(grouped_samples)<-c("mz","intensity","RT","peak ID", "group ID","series level","to ID","m/z increment","RT increment")
		grouping<-data.frame(group1,group2,group3,group4,group5,group6,group6-group5);
		names(grouping)<-c("group ID","peak IDs","m/z increment","RT increment","min. RT in series","max. RT in series","max.-min. RT");
		# store parameters used ##################################################
		parameters<-list(0)
		parameters[[1]]<-elements;
		parameters[[2]]<-charge
		parameters[[3]]<-use_C
		parameters[[4]]<-minmz
		parameters[[5]]<-maxmz
		parameters[[6]]<-minrt
		parameters[[7]]<-maxrt
		parameters[[8]]<-ppm	
		parameters[[9]]<-mztol
		parameters[[10]]<-rttol
		parameters[[11]]<-minlength
		names(parameters)<-c("elements","charge","use_C","minmz","maxmz","minrt","maxrt","ppm","mztol","rttol","minlength")
		##########################################################################
		homol<-list(grouped_samples,parameters,grouping,mzfilter,peakID,relat);
		names(homol)<-c("Homologue Series","Parameters","Peaks in homologue series","m/z-Filter used","Peaks per level","Single relations")
		##########################################################################
		return(homol); ###########################################################
		##########################################################################
  	}else{
		stop("no series detected");
	}
}















