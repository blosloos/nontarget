combine <-
function(
	pattern,
	adduct,
	homol = FALSE,
	rules = c(TRUE, TRUE, FALSE, FALSE),
	dont = FALSE
){

    ############################################################################
	# rules[1]: incorporate interfering peaks into a group? ####################
    # rules[2]: deal with interfering peaks also seperately ? ##################
    # rules[3]: remove single-peaked components? ###############################
    # rules[4]: remove components not part of a homologue series ? #############
    # warning(1): only one adduct found per main isotope group? ################
    # warning(2): adduct ranks ok? #############################################
    # warning(3): any interferring peaks? ######################################
    ############################################################################
    # (0.1) checks #############################################################
	if(!all(is.logical(rules))) stop("rules must be a vector of three logical values!")
	if(length(rules) != 4) stop("rules must be a vector of length 4!")
	if(!rules[1] & !rules[2]) warning("Beware: rules[1] and rules[2] both set to FALSE: peaks marked as interfering may not be combined into components!")
	if(!is.logical(dont)) stop("dont must be logical!")	
    if(length(pattern) > 1 & length(adduct) > 1) if(!identical(pattern[[1]][,1], adduct[[1]][,1])) stop("Different data sets pattern<->adduct used for combining!")
    if(length(homol) > 1 & length(adduct) > 1) if(!identical(homol[[1]][,1], adduct[[1]][,1])) stop("Different data sets homol<->adduct used for combining!")
    if(length(homol) > 1 & length(pattern) > 1) if(!identical(homol[[1]][,1], pattern[[1]][,1])) stop("Different data sets homol<->pattern used for combining!") 
    if(rules[4] & length(homol) < 1) stop("rule3 = TRUE not applicable if no correct homologue series are provided!")
    if(dont != FALSE){
		if(length(dont) > 4) stop("Invalid dont argument")
		if(any(dont > 4)) stop("Invalid dont argument")
		if(any(dont < 1)) stop("Invalid dont argument")
    }
    if(length(pattern) > 1 & length(adduct) > 1){
		if(pattern[[11]][1] != "FALSE"){
			if((adduct[[2]][4] == "negative" & any(as.numeric(pattern[[11]])>0)) || (adduct[[2]][4]=="positive" & any(as.numeric(pattern[[11]])<0))){
				warning("Are charges of adduct vs. pattern groups consistent?")
			}
		}
	}
    ###########################################################################
    # (0.2) local function definitions ########################################
    rankit <- function(a, b, inttol){
		# on input vector a ###################################################
		a_low <- a * (1 - (inttol * 2))
		a_up <- a * (1 + (inttol * 2))
		orda <- order(a, decreasing = TRUE)
		ord_a <- rep(0, length(a))
		dor <- length(a)
		ord_a[orda[1]] <- length(a)
		for(j in 2:length(a)){
           if(a_up[orda[j]] < a_low[orda[(j-1)]]) dor <- c(dor - 1)
           ord_a[orda[j]]<-dor;
		}
		# on input vector b ###################################################
		b_low <- b * (1 - (inttol * 2))
		b_up <- b * (1 + (inttol * 2))
		ordb <- order(b, decreasing = TRUE)
		ord_b <- rep(0, length(b))
		dor <- c(length(b))
		ord_b[ordb[1]] <- c(length(b))
		for(j in 2:length(b)){ 
			if(b_up[ordb[j]] < b_low[ordb[(j-1)]]) dor <- c(dor - 1)
			ord_b[ordb[j]] <- dor;
		}
		#######################################################################
		if(sum(ord_a) >= sum(ord_b))  return(all(ord_a >= ord_b)) else return(all(ord_a <= ord_b))
    }
    cat("done.");
    ##########################################################################

    ##########################################################################
    # combine # 1 ############################################################
    cat("\n (1) Assemble lists...");
    comp1a <- c();    # = ID of pattern group 
    comp1b <- c();    # = IDs of peaks in main pattern group
    comp2a <- c();    # = ID of adduct group(s)
    comp2b <- c();    # = IDs of peaks in adduct group
    comp2c <- c();    # = main adduct, i.e. that of comp1b
    comp2d <- c();    # = other adducts, i.e. those of comp2b
    comp3 <- c();     # = ID of homologue serie(s) 
    comp4 <- c();     # = IDs of interfering peaks
    comp5 <- c();     # = IDs of interfering pattern groups
    comp6 <- c();     # = IDs of interfering adduct groups
    comp7 <- c();     # = consistent? -> ok/(1,2,3,4,5,6)
    comp8 <- c();     # = isotope relations
    comp9 <- c();     # = charge levels
    cat("\n (2) Combine...");
    if(length(pattern[[1]]) > 1){
		no3<-rep(TRUE,length(pattern[[1]][,1]));
		intord<-order(pattern[[1]][,2],decreasing=TRUE);
    }else{
		no3<-rep(TRUE,length(adduct[[1]][,1]));
		intord<-order(adduct[[1]][,2],decreasing=TRUE);
    }
    if(length(homol[[1]])>1){ # save homologue IDs per peak ID ###############
		getit<-rep(0,length(no3));
		for(i in 1:length(homol[[3]][,1])){
			compit<-as.numeric(strsplit(as.character(homol[[3]][i,2]),",")[[1]]);
			for(m in 1:length(compit)){
				getit[compit[m]]<-paste(getit[compit[m]],"/",i,sep="");
			}
		}
    }
	#options(warn = 2) # enable this option to make R stop() on any warnings
    ##########################################################################
	cat("\n")
	pBar <- txtProgressBar( min = 0, max = length(intord), style = 3 )
    for(i in 1:length(intord)){ # along decreasing intensity #################
#if(i == 4) stop()
		setTxtProgressBar(pBar, i, title = NULL, label = NULL)
		if(no3[intord[i]] == FALSE) next # if not yet included in some group ##
        ######################################################################
        comp7 <- c(comp7,"-"); # initialized as consistent 
        ######################################################################
        # on isotope pattern #################################################
        ######################################################################
        if(length(pattern[[1]])>1){ # if available ...
            if(pattern[[1]][intord[i], "group ID"] != "0"){
                allpeaks <- intord[i]
                get1 <- as.numeric(strsplit(as.character(pattern[[1]][allpeaks, "group ID"]), "/")[[1]]); # group ID(s) of main peak
				# get1: if several groups, the first one should be the nesting group at highest charge	
				where <- unique(unlist(
					lapply(get1, function(x) which(grepl(paste("/", x, "/", sep = ""), as.character(pattern[[3]][,1]), fixed = TRUE)))
				))	
				# should be length(where)==1, unless group(s) is nested into several unmergeable groups of higher charge
				# if so, combine them and issue a warning = 5
				if(length(where) > 1) comp7[length(comp7)] <- paste0(comp7[length(comp7)], ",5");
				# get all peak IDs from these (possibly nested) group(s) ######
				get2 <- unique(unlist(strsplit(as.character(pattern[[3]][where, "peak IDs"]), ",")))
				# ID of pattern group(s) of this most intense peak ############
				comp1a <- c(comp1a, as.character(pattern[[1]][allpeaks, "group ID"])); 
				# IDs of all related peaks #################################### 
				compit <- paste(get2, collapse = ",")
                comp1b <- c(comp1b, compit);                                       # = IDs of peaks in pattern group
                allpeaks <- c(allpeaks, as.numeric(get2));
                allpeaks <- unique(allpeaks);				
                no3[allpeaks] <- FALSE;
				# get isotope relations among peaks ###########################
                get6 <- c()
                for(k in 1:length(allpeaks)){
                    if(pattern[[1]][allpeaks[k],7] != "0"){
						get3 <- as.numeric(strsplit(as.character(pattern[[1]][allpeaks[k], "to ID"]),"/")[[1]]);	# to ID
						get4 <- strsplit(as.character(pattern[[1]][allpeaks[k], "isotope(s)"]),"/")[[1]]; 			# isotope(s)
						get5 <- strsplit(as.character(pattern[[1]][allpeaks[k], "mass tolerance"]),"/")[[1]];		# mass tolerance
						for(y in 1:length(get3)){
							if(any(allpeaks == get3[y])){
								get6 <- c(get6, paste(get4[y], "(", get5[y], ")", sep = ""));
							}
						}
                    }
                }
				get6 <- paste(get6, collapse = "/")			
				comp8 <- c(comp8, get6);	# comp8 = isotope relations							
				get7 <- unique(unlist(strsplit(as.character(pattern[[3]][where, "charge level"]),"/")))
				get7 <- sort(get7, decreasing = TRUE)
				get7 <- paste(get7, collapse = "/")
				comp9 <- c(comp9, get7)	# comp9 = charge levels
				#if(length(comp8)!=length(comp9)){stop("\n debug_A")} 
				if(length(where)>1){
					comp7[length(comp7)]<-paste(comp7[length(comp7)],",4",sep="")
                }      
            }else{
				comp1a<-c(comp1a,"-");   # = ID of pattern group 
                comp1b<-c(comp1b,as.character(intord[i]));   # = IDs of peaks in pattern group
                comp8<-c(comp8,"-"); # = isotope relations found
                comp9<-c(comp9,"-"); # = isotope relations found
				#if(length(comp8)!=length(comp9)){stop("\n debug_B")}
                allpeaks<-c(intord[i]);
                no3[intord[i]]<-FALSE;
            }
        }else{
            comp1a <- c(comp1a, "-")   # = ID of pattern group 
            comp1b <- c(comp1b, as.character(intord[i]))   # = IDs of peaks in pattern group
            comp8 <- c(comp8, "-") # = isotope relations found
            comp9 <- c(comp9, "-") # = isotope relations found
			#if(length(comp8) != length(comp9)) stop("\n debug_C")
            allpeaks <- c(intord[i]);
            no3[intord[i]] <- FALSE;
        }
		######################################################################
        # on adduct pattern ##################################################
        ######################################################################
        if(length(adduct[[1]]) > 1){ # if available ...
            thesepeaks <- c();                                         
            addID <- c();
            addfrom <- c();
            addto <- c();
            for(m in 1:length(allpeaks)){
				get1 <- as.numeric(strsplit(as.character(adduct[[1]][allpeaks[m], 5]), "/")[[1]]);
				get3 <- strsplit(as.character(adduct[[1]][allpeaks[m], 7]), "//")[[1]];
				if(get1[1] != 0){
					for(n in 1:length(get1)){
						get2 <- as.numeric(strsplit(as.character(adduct[[3]][get1[n], 2]), ",")[[1]]);
						for(k in 1:length(get2)){
							if(any(allpeaks == get2[k]) == FALSE){
								thesepeaks <- c(thesepeaks, get2[k]);
							};
						};
						addID <- c(addID, get1[n]);
					};
					for(n in 1:length(get3)){
						get4 <- strsplit(get3[n], "<->")[[1]];
						addfrom <- c(addfrom, get4[1]);
						addto <- c(addto, get4[2]);
					};
				};
            };
            if(length(thesepeaks)>0){
                thesepeaks <- unique(thesepeaks);
                addID <- unique(addID);
                # store adduct group IDs #######################################
                compit <- paste(as.character(addID), collapse = "/");
                comp2a <- c(comp2a,compit);   # = ID of adduct group(s)
                # store adduct peak IDs ########################################
                compit <- paste(as.character(thesepeaks), collapse = ",");
                comp2b <- c(comp2b,compit);   # = IDs of peaks in adduct group
                # store adduct from ############################################
				compit <- paste(as.character(unique(addfrom)), collapse = ",");
                comp2c <- c(comp2c,compit);   # = main adduct, i.e. that of comp1b
                if(length(addfrom)>1) comp7[length(comp7)] <- paste0(comp7[length(comp7)], ",1")
                # store addut to ###############################################
				compit <- paste(as.character(unique(addto)), collapse = ",");
                comp2d <- c(comp2d,compit);   # = other adducts, i.e. those of comp2b
                ################################################################
                # check adduct plausbility #####################################
                # (1) same adducts? ... done ###################################
                if( length(addfrom) == 1 & length(allpeaks) > 1 ){ 
					# allpeaks>1 also ensures that pattern data set available! #
					# (0) derive raltion matrix ################################
					mat <- matrix(ncol = length(allpeaks), nrow = (length(addto) + 1), 0);
					colnames(mat) <- allpeaks
					rownames(mat) <- c(addfrom,addto)
					for(z in 1:length(allpeaks)){
						mat[1, z] <- adduct[[1]][allpeaks[z], 2];
						geta <- as.numeric(strsplit(as.character(adduct[[1]][allpeaks[z], 6]),"/")[[1]]);
						getb <- strsplit(as.character(adduct[[1]][allpeaks[z], 7]), "//")[[1]];
						if(geta[1] != 0){
							for(y in 1:length(geta)){
								mat[rownames(mat) == strsplit(getb[y], "<->")[[1]][2],z] <- adduct[[1]][geta[y], 2];
							};
						};
					};
					cutint <- pattern[[2]][[7]];
					inttol <- pattern[[2]][[6]];
					mat[(mat * (1 - 2 * inttol)) < cutint] <- 0
					# (2) adduct intensity ranks ok? ###########################
					for(z in 2:length(mat[,1])) if(rankit(mat[1,],mat[z,],inttol)==FALSE) comp7[length(comp7)]<-paste0(comp7[length(comp7)], ",2")
					# (3) relative adduct intensities ok? ######################
					# ... skipped ##############################################
                }
                ################################################################
                allpeaks <- c(allpeaks, thesepeaks);
                allpeaks <- unique(allpeaks);
                no3[allpeaks] <- FALSE; 
            }else{
                comp2a <- c(comp2a, "-");   # = ID of adduct group(s)
                comp2b <- c(comp2b, "-");   # = IDs of peaks in adduct group
                comp2c <- c(comp2c, "-");   # = main adduct, i.e. that of comp1b
                comp2d <- c(comp2d, "-");   # = other adducts, i.e. those of comp2b         
            }
        }else{
            comp2a <- c(comp2a, "-");   # = ID of adduct group(s)
            comp2b <- c(comp2b, "-");   # = IDs of peaks in adduct group
            comp2c <- c(comp2c, "-");   # = main adduct, i.e. that of comp1b
            comp2d <- c(comp2d, "-");   # = other adducts, i.e. those of comp2b         
        }
        ######################################################################
        # on interfering peaks ###############################################
        allpeaks <- unique(allpeaks)
        oldpeaks <- allpeaks;
        if((length(allpeaks) > 1) & rules[1]){
            IDpeak <- c("")
            IDpat <- c("")
            IDadd <- c("")             
            newpeaks1 <- allpeaks
            newpeaks2 <- c()		
            while(length(newpeaks1)){
				for(j in 1:length(newpeaks1)){
					##################################################################
					# find other pattern groups ######################################
					if(length(pattern)>1){
						that <- strsplit(as.character(pattern[[1]][newpeaks1[j], 5]), "/")[[1]];
						if(that[1] != "0"){
							for(n in 1:length(that)){
								if( !grepl(paste0("/", that[n], "/"), as.character(comp1a[i])) & grepl(paste0("/", that[n], "/"), IDpat)){
									this <- as.numeric(strsplit(as.character(pattern[[3]][grepl(paste0("/",that[n],"/"),pattern[[3]][,1]), 2]), ",")[[1]])
									for(m in 1:length(this)){
										doit <- FALSE;
										if(!any(oldpeaks == this[m])){
											IDpeak <- paste0(IDpeak, ",", as.character(this[m]));
											newpeaks2 <- c(newpeaks2, this[m]);
											doit <- TRUE;
										}
										if(doit) IDpat <- paste0(IDpat, paste0("/", as.character(that[n]), "/"));
									}
								}
							}
						}
					}
					# find other adduct groups #######################################
					if(length(adduct)>1){
						that <- strsplit(as.character(adduct[[1]][newpeaks1[j], 5]), "/")[[1]];
						if(that[1] != "0"){
							for(n in 1:length(that)){
								if(!grepl(paste0("/", that[n], "/"), as.character(comp2a[i])) & grepl(paste0("/", that[n], "/"), IDadd)){
									this <- as.numeric(strsplit(as.character(adduct[[3]][grepl(paste0("/", that[n], "/"), adduct[[3]][, 1]), 2]), ",")[[1]])
									for(m in 1:length(this)){
											doit <- FALSE;
											if(!any(oldpeaks == this[m])){
												IDpeak <- paste0(IDpeak, ",", as.character(this[m]));
												newpeaks2 <- c(newpeaks2, this[m]);
												doit <- TRUE;
											}
											if(doit) IDadd <- paste0(IDadd, paste0("/", as.character(that[n]), "/"));
											
									}
								}
							}
						}
					}
					##################################################################
				}		
				newpeaks2 <- unique(newpeaks2);
				oldpeaks <- c(oldpeaks, newpeaks2)
				newpeaks1 <- newpeaks2;
				newpeaks2 <- c();
            } # while newpeaks1
            if(nchar(IDpeak)){
				comp7[length(comp7)] <- paste0(comp7[length(comp7)], ",3")
				comp4 <- c(comp4, substr(IDpeak, 2, nchar(IDpeak)));
            }else{
				comp4 <- c(comp4, "-");
            }
			if(nchar(IDpat)){
				comp5 <- c(comp5, IDpat);
            }else{
				comp5 <- c(comp5, "-");
            }
            if(nchar(IDadd)){
				comp6 <- c(comp6, IDadd);
            }else{
				comp6 <- c(comp6, "-");
            }  
        }else{
			comp4 <- c(comp4,"-");
            comp5 <- c(comp5,"-");
            comp6 <- c(comp6,"-");
        }
        ######################################################################
        # rules[2]: deal with interfering peaks also seperately ? ############
        if(!rules[2]){
			allpeaks <- oldpeaks;
			no3[allpeaks] <- FALSE;
        }
        ######################################################################
        # on homologue series ################################################          
        if(length(homol[[1]])>1){ # if available ...
            compit<-c();
            for(m in 1:length(allpeaks)){
				if(getit[allpeaks[m]]!="0") compit <- paste0(compit, substr(getit[allpeaks[m]], 2, nchar(getit[allpeaks[m]])));
            }
            if(length(compit)) comp3 <- c(comp3, compit) else comp3<-c(comp3,"-");    # = ID of homologue serie(s)             
        }else{
			comp3<-c(comp3,"-");    # = ID of homologue serie(s) 
        }
        ######################################################################
#cat("*")
    }
	close(pBar)
	cat("\n")
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3)
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3,comp4,comp5,comp6,comp7)
    cat("done.");
    ############################################################################
        
    ############################################################################
    cat("\n (3) Apply rules 2&3...");
    # remove single-peaked components? #########################################
    if(rules[3]){
		getit <- c()
		for(i in 1:length(comp1b)){
			if(length(strsplit(comp1b[i], ",")[[1]]) == 1){
				if(comp2b[i] == "-" & comp3[i] == "-" ){
					getit <- c(getit, i);
				}
			}
		}
		comp1a <- comp1a[-getit]    # = ID of pattern group 
		comp1b <- comp1b[-getit]    # = IDs of peaks in pattern group
		comp2a <- comp2a[-getit]    # = ID of adduct group(s)
		comp2b <- comp2b[-getit]    # = IDs of peaks in adduct group
		comp2c <- comp2c[-getit]    # = main adduct, i.e. that of comp1b
		comp2d <- comp2d[-getit]    # = other adducts, i.e. those of comp2b
		comp3 <- comp3[-getit]      # = ID of homologue serie(s) 
		comp4 <- comp4[-getit]      # = IDs of interfering peaks
		comp5 <- comp5[-getit]      # = IDs of interfering pattern groups
		comp6 <- comp6[-getit]      # = IDs of interfering adduct groups
		comp7 <- comp7[-getit]      # = consistent? -> ok/(1,2,3,4,5,6)
		comp8 <- comp8[-getit]      # = isotope relations in main pattern group
		comp9 <- comp9[-getit]      # = charge level
    }
    ############################################################################
    # remove components not part of a homologue series ? #######################
    if(rules[4] & length(homol)>1){
		getit <- c()
		for(i in 1:length(comp3)) if(comp3[i] != "-") getit <- c(getit, i);
		comp1a <- comp1a[getit]    # = ID of pattern group 
		comp1b <- comp1b[getit]    # = IDs of peaks in pattern group
		comp2a <- comp2a[getit]    # = ID of adduct group(s)
		comp2b <- comp2b[getit]    # = IDs of peaks in adduct group
		comp2c <- comp2c[getit]    # = main adduct, i.e. that of comp1b
		comp2d <- comp2d[getit]    # = other adducts, i.e. those of comp2b      
		comp3 <- comp3[getit]      # = ID of homologue serie(s) 
		comp4 <- comp4[getit]      # = IDs of interfering peaks
		comp5 <- comp5[getit]      # = IDs of interfering pattern groups
		comp6 <- comp6[getit]      # = IDs of interfering adduct groups
		comp7 <- comp7[getit]      # = consistent? -> ok/(1,2,3,4,5,6)
		comp8 <- comp8[getit]      # = isotope relations in main pattern group
		comp9 <- comp9[getit]      # = charge level
    };
    # data.frame(comp1a,comp1b,comp2a,comp2b,comp3,comp4,comp5,comp6,comp7)
    ############################################################################
    if(!rules[3] & !rules[4]) cat("skipped.") else cat("done.")
    ############################################################################

    ############################################################################
    # summarize / return #######################################################
    cat("\n (4) Generate output...");
    # main list ################################################################
    comp7[comp7!="-"] <- substr(comp7[comp7 != "-"], 3, nchar(comp7[comp7 != "-"]))
    comps <- data.frame(seq(1, length(comp1a), 1), comp1a, comp1b, comp2a, comp2b, comp3, comp4, comp5, comp6, comp7, comp2c, comp2d, stringsAsFactors = FALSE)
    names(comps)<-c("Component ID |","ID pattern group |","ID pattern peaks |","ID adduct group(s) |","ID adduct peaks |",
    "ID homologue series |","ID interfering peaks |","ID interfering pattern group(s) |",
    "ID interfering adduct group(s) |","Warnings |","pattern group adduct|","adduct group adduct(s) |");
    # exclude components with warnings in comp 7 ###############################
    if(length(dont)){
		if(dont[1]){
			excl <- rep(TRUE, length(comp7));
			for(i in 1:length(comp7)){
				if(comp7[i] != "-"){  
					this <- as.numeric(strsplit(as.character(comp7[i]),",")[[1]])
					for(j in 1:length(this)) if(any(dont == this[j])) excl[i]<-FALSE
					
				
				}
			}
			comps <- comps[excl,]
			comp8 <- comp8[excl]
			comp9 <- comp9[excl]
		} # if done
    } # if done
    # sort by decreasing intensity of ALL peaks and ############################ 
    # calculate mean component size ############################################
    # index peaks in components ################################################
    int <- c();
    intID <- c();
    intmz <- c();
    intrt <- c();
    numb <- c();
    allin <- c();
    partin <- c();
    for(i in 1:length(comps[,1])){
		this <- c();
		if(comps[i,3] != "-"){
			this <- c(this, as.numeric(strsplit(as.character(comps[i,3]),",")[[1]]));
		}
		if(comps[i,5] != "-"){
			this <- c(this, as.numeric(strsplit(as.character(comps[i,5]),",")[[1]]));
		}
		if(comps[i,7] != "-"){
			this<-c(this, as.numeric(strsplit(as.character(comps[i,7]),",")[[1]]));
			partin<-c(partin,as.numeric(strsplit(as.character(comps[i,7]),",")[[1]]));
		}
		this<-unique(this);
		numb<-c(numb,length(this));
		if(length(pattern)>1){
			int<-c(int,max(pattern[[1]][this,2])[1]);
			intID<-c(intID,pattern[[1]][this,4][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);
			intmz<-c(intmz,pattern[[1]][this,1][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);
			intrt<-c(intrt,pattern[[1]][this,3][pattern[[1]][this,2]==max(pattern[[1]][this,2])[1]][1]);        
		}else{
			int<-c(int,max(adduct[[1]][this,2])[1]);
			intID<-c(intID,adduct[[1]][this,4][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);
			intmz<-c(intmz,adduct[[1]][this,1][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);
			intrt<-c(intrt,adduct[[1]][this,3][adduct[[1]][this,2]==max(adduct[[1]][this,2])[1]][1]);        
		}
		allin <- c(allin, this)
    }
    comps <- cbind(comps, intID, intmz, int, intrt, comp8, comp9)
    names(comps)[13:18] <- c("HI peak ID |", "HI m/z |", "Highest intensity (HI) |", "HI RT |", "Isotope(s) d(m/z) |", "z")
    comps <- comps[order(int, decreasing = TRUE),]
    comps[,1] <- seq(1, length(comps[,1]),1);
    numb <- mean(numb)
    if(length(pattern) > 1) no1 <- rep(FALSE,length(pattern[[1]][,1])) else no1 <- rep(FALSE,length(adduct[[1]][,1]));
    
    no1[allin] <- TRUE
    if(length(pattern) > 1) no2<-rep(FALSE,length(pattern[[1]][,1])) else no2 <- rep(FALSE,length(adduct[[1]][,1]))
    
    no2[partin] <- TRUE;
    out<-c(round(length(no1[no1])/length(no1),digits=3),round(length(no1[no1]),digits=0),round(length(no2[no2]),digits=0),round(numb,digits=2));
    names(out)<-c("Fraction of peaks in components","Number of peaks in components","Number of interfering peaks","Mean number of peaks in components");
    ############################################################################
    # get parameters ###########################################################
    if(length(pattern)>1){
		param <- pattern[[2]];
		param <- param[-1];
		if(length(adduct)>1){
			param[1] <- max(adduct[[2]][1],pattern[[2]][1],pattern[[2]][2]);
		}else{
			param[1] <- max(pattern[[2]][1],pattern[[2]][2]);
		}
    }else param <- adduct[[2]];
    ############################################################################
    comp<-list(comps,FALSE,FALSE,FALSE,no1,out,param);
    if(length(pattern)>1) comp[[2]] <- pattern[[1]] else comp[[2]] <- "No isotope pattern groups available."
    if(length(adduct)>1) comp[[3]] <- adduct[[1]] else comp[[3]] <- "No adduct groups available."
    if(length(homol)>1) comp[[4]] <- homol[[1]] else comp[[4]] <- "No homologue series available."
    names(comp) <- c("Components", "pattern peak list", "adduct peak list", "homologue list", "Peaks in components", "Summary", "Parameters");
    cat("done.\n\n");
    return(comp);
    ############################################################################
}
