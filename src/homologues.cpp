#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <vector>
#include <algorithm>

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RMATRIX2(m,i,j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RVECTOR(m,i) (REAL(m)[i])
#define RVECTOR2(m,i) (INTEGER(m)[i])
#define RRow(m) (INTEGER(GET_DIM(m))[0])
#define RCol(m) (INTEGER(GET_DIM(m))[1])

extern "C"{

/******************************************************************************/
/* filter homologue triplets **************************************************/
/******************************************************************************/

    SEXP homol_triplet(
        SEXP peaklist3,
		SEXP dist_ID,
		SEXP dist_dist,
		SEXP triplets,
		SEXP tupeldo,
		SEXP peaklist4,
        SEXP use,
        SEXP max_delmz,
        SEXP rttol,
        SEXP diagno
   ){

        PROTECT(peaklist3 = AS_NUMERIC(peaklist3));
        PROTECT(dist_ID = AS_NUMERIC(dist_ID));
        PROTECT(dist_dist = AS_NUMERIC(dist_dist));
        PROTECT(triplets = AS_NUMERIC(triplets));
        PROTECT(tupeldo = AS_NUMERIC(tupeldo));
        PROTECT(peaklist4 = AS_NUMERIC(peaklist4));
        PROTECT(use = AS_NUMERIC(use));
        PROTECT(max_delmz = AS_NUMERIC(max_delmz));
        PROTECT(rttol = AS_NUMERIC(rttol));

        int n,m,nrow;
        nrow = RRow(triplets);
        int leng_dist = LENGTH(dist_dist);
        double min_mass, min_mass_LB, min_mass_UB, max_mass, max_mass_LB, max_mass_UB, delmz_use;
        double rtdif1, rtdif1_LB, rtdif1_UB, rtdif2, rtdif2_LB, rtdif2_UB;
        delmz_use = RVECTOR(peaklist4, int(NUMERIC_VALUE(use) - 1));

        int a = 0, b = 0; // precheck for occurrence of positive AND negative m/z distances
        for(n = 0; n < leng_dist; n++){
            if(RVECTOR(dist_dist, n) < 0){
                a++;
            }else{
                b++;
            }
        }
        if((a == 0)||(b == 0)){
            UNPROTECT(9);
            return(R_NilValue);
        }

        for(n=0;n<(leng_dist-1);n++){ // find matching (1) mass distance -> (2) mass defect changes -> (3) RT change
            if((fabs(RVECTOR(dist_dist,(n+1)))-fabs(RVECTOR(dist_dist,n)))>(NUMERIC_VALUE(max_delmz))){ // pre-check on m/z uncertainty
                continue;
            }
            if(RVECTOR(dist_dist,n)<0){ // initiate a search at a negative distance ...
                min_mass=fabs(RVECTOR(dist_dist,n));
                min_mass_UB=(min_mass+RVECTOR(peaklist4,int(RVECTOR(dist_ID,n)-1))+delmz_use);
                for(m=(n+1);m<leng_dist;m++){
                    if((fabs(RVECTOR(dist_dist,m))-fabs(RVECTOR(dist_dist,n)))>(NUMERIC_VALUE(max_delmz))){
                        break;
                    }
                    if(RVECTOR(dist_dist,m)>0){
                        // (1) mass distance overlap?
                        max_mass=fabs(RVECTOR(dist_dist,m));
                        max_mass_LB=(max_mass-RVECTOR(peaklist4,int(RVECTOR(dist_ID,m)-1))-delmz_use);
                        if(max_mass_LB<min_mass_UB){
                            rtdif1=(RMATRIX(peaklist3,int(NUMERIC_VALUE(use)-1),1)-RMATRIX(peaklist3,int(RVECTOR(dist_ID,n)-1),1));
                            rtdif1_LB=(rtdif1-NUMERIC_VALUE(rttol));
                            rtdif2=(RMATRIX(peaklist3,int(RVECTOR(dist_ID,m)-1),1)-RMATRIX(peaklist3,int(NUMERIC_VALUE(use)-1),1));
                            rtdif2_UB=(rtdif2+NUMERIC_VALUE(rttol));
                            if(!(rtdif2_UB<rtdif1_LB)){
                                rtdif1_UB=(rtdif1+NUMERIC_VALUE(rttol));
                                rtdif2_LB=(rtdif2-NUMERIC_VALUE(rttol));
                                if(!(rtdif2_LB>rtdif1_UB)){
                                    // store IDs of peaks in triplet
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),0)=RVECTOR(dist_ID,n); // CHANGED to m below
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),1)=NUMERIC_VALUE(use);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),2)=RVECTOR(dist_ID,m); // CHANGED to n below
                                    // store mass bounds
                                    min_mass_LB=(min_mass-RVECTOR(peaklist4,int(RVECTOR(dist_ID,n)-1))-delmz_use);
                                    max_mass_UB=(max_mass+RVECTOR(peaklist4,int(RVECTOR(dist_ID,m)-1))+delmz_use);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),3)=std::max(min_mass_LB,max_mass_LB);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),4)=std::min(min_mass_UB,max_mass_UB);
                                    // store RT bounds
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),5)=std::min(rtdif1_LB,rtdif2_LB);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),6)=std::max(rtdif1_UB,rtdif2_UB);
                                    //
                                    RVECTOR(tupeldo,0)=(RVECTOR(tupeldo,0)+1);
                                    if(NUMERIC_VALUE(tupeldo)>=nrow){
                                        Rprintf("vec_size too small!");
                                        UNPROTECT(9);
                                        return(R_NilValue);
                                    }
                                }
                            }
                        }
                    }
                }
            }else{ // ... or a positive one
                min_mass=fabs(RVECTOR(dist_dist,n));
                min_mass_UB=(min_mass+RVECTOR(peaklist4,int(RVECTOR(dist_ID,n)-1))+delmz_use);
                for(m=(n+1);m<leng_dist;m++){
                    if(((RVECTOR(dist_dist,m))-fabs(RVECTOR(dist_dist,n)))>(NUMERIC_VALUE(max_delmz))){ // RVECTOR(dist_dist,n) is >0
                        break;
                    }
                    if(RVECTOR(dist_dist,m)<0){
                        max_mass=fabs(RVECTOR(dist_dist,m));
                        max_mass_LB=(max_mass-RVECTOR(peaklist4,int(RVECTOR(dist_ID,m)-1))-delmz_use);
                        if(max_mass_LB<min_mass_UB){ // (1) mass distance overlap? only check upper vs. lower bound, because sorted by increasing m/z distance
                            // (3) RT change within tolerance?
                            rtdif1=(RMATRIX(peaklist3,int(RVECTOR(dist_ID,n)-1),1)-RMATRIX(peaklist3,int(NUMERIC_VALUE(use)-1),1));
                            rtdif1_LB=(rtdif1-NUMERIC_VALUE(rttol));
                            rtdif2=(RMATRIX(peaklist3,int(NUMERIC_VALUE(use)-1),1)-RMATRIX(peaklist3,int(RVECTOR(dist_ID,m)-1),1));
                            rtdif2_UB=(rtdif2+NUMERIC_VALUE(rttol));
                            if(!(rtdif2_UB<rtdif1_LB)){
                                rtdif1_UB=(rtdif1+NUMERIC_VALUE(rttol));
                                rtdif2_LB=(rtdif2-NUMERIC_VALUE(rttol));
                                if(!(rtdif2_LB>rtdif1_UB)){
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),0)=RVECTOR(dist_ID,m);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),1)=NUMERIC_VALUE(use);
                                   RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),2)=RVECTOR(dist_ID,n);
                                    // store mass bounds
                                    min_mass_LB=(min_mass-RVECTOR(peaklist4,int(RVECTOR(dist_ID,n)-1))-delmz_use);
                                    max_mass_UB=(max_mass+RVECTOR(peaklist4,int(RVECTOR(dist_ID,m)-1))+delmz_use);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),3)=std::max(min_mass_LB,max_mass_LB);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),4)=std::min(min_mass_UB,max_mass_UB);
                                    // store RT bounds
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),5)=std::min(rtdif1_LB,rtdif2_LB);
                                    RMATRIX(triplets,int(NUMERIC_VALUE(tupeldo)-1),6)=std::max(rtdif1_UB,rtdif2_UB);
                                    RVECTOR(tupeldo,0)=(RVECTOR(tupeldo,0)+1);
                                    if(NUMERIC_VALUE(tupeldo)>=nrow){
                                        Rprintf("vec_size too small!");
                                        UNPROTECT(9);
                                        return(R_NilValue);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // for all dist_dist

        UNPROTECT(9);
        return(R_NilValue);

   }


/******************************************************************************/
/* build all feasible 2-combinations from (n-1) to n tuples of homologues *****/
/* tuples & bounds must be sorted by increasing LB of first column in bounds **/
/* here, e.g.: homologue mass differences *************************************/
/* tuples may have variable lengths *******************************************/
/******************************************************************************/

    SEXP combine_tuple(
		SEXP tuples,
		SEXP bounds,
		SEXP keeper,
		SEXP mat_size
    ){

        PROTECT(tuples = AS_NUMERIC(tuples));
        PROTECT(bounds = AS_NUMERIC(bounds));
        PROTECT(keeper = AS_NUMERIC(keeper));
        PROTECT(mat_size = AS_NUMERIC(mat_size));

        int n,m,nrow,ncol,ncol2,k,skip_a=0,skip_b=0,doing=0,siz;
        nrow=RRow(tuples);
        ncol=RCol(tuples);
        ncol2=RCol(bounds);

        siz=((ncol+1)+ncol2+2);
        SEXP merged_tuples;
        PROTECT(merged_tuples = allocMatrix(REALSXP, (nrow*NUMERIC_VALUE(mat_size)), siz));
        for(n=0;n<(nrow*NUMERIC_VALUE(mat_size));n++){
            for(m=0;m<siz;m++){
                RMATRIX(merged_tuples,n,m)=0;
            }
        }

        for(n=0;n<(nrow-1);n++){
            for(m=(n+1);m<nrow;m++){
                if( RMATRIX(bounds,n,1) < RMATRIX(bounds,m,0) ){ // check bounds[,1:2]: UB<LB?
                    break;
                }
                if( !(RMATRIX(bounds,m,3) < RMATRIX(bounds,n,2)) ){ // check bounds[,3:4]: UB<LB?
                    if( !(RMATRIX(bounds,m,2) > RMATRIX(bounds,n,3)) ){ // check bounds[,3:4]: LB>UB?
                        skip_a=0; // peak matching - downward
                        for(k=(ncol-2);k>=0;k--){ // check peak matching - downward
                            if(RMATRIX(tuples,n,(k+1))!=RMATRIX(tuples,m,k)){
                                skip_a=1;
                                break;
                            }
                        }
                        if(skip_a==1){ // check peak matching - upward
                            skip_b=0; // peak matching - upward
                            for(k=1;k<ncol;k++){
                                if(RMATRIX(tuples,n,(k-1))!=RMATRIX(tuples,m,k)){
                                    skip_b=1;
                                    break;
                                }
                            }
                        }else{
                            skip_b=1;
                        }
                        if((skip_a==1)&(skip_b==1)){ // peaks don`t match, neither down- nor upward
                            continue;
                        }
                        // create new n+1 tuple:
                        if(skip_a==0){ // either downward
                            RMATRIX(merged_tuples,doing,0)=RMATRIX(tuples,n,0);
                            for(k=0;k<ncol;k++){
                                RMATRIX(merged_tuples,doing,(k+1))=RMATRIX(tuples,m,k);
                            }
                        }
                        if(skip_b==0){ // ... or upward
                            RMATRIX(merged_tuples,doing,ncol)=RMATRIX(tuples,n,(ncol-1));
                            for(k=0;k<ncol;k++){
                                RMATRIX(merged_tuples,doing,k)=RMATRIX(tuples,m,k);
                            }
                        }
                        if(doing>=(nrow*NUMERIC_VALUE(mat_size))){ // more combination results than length of original matrix? then extend:
                            Rprintf("ERROR: matrix out of bounds! increase mat_size");
                            break;
                        }
                        // update bounds
                        // m/z differences
                        RMATRIX(merged_tuples,doing,(ncol+1))=std::max(RMATRIX(bounds,n,0),RMATRIX(bounds,m,0));
                        RMATRIX(merged_tuples,doing,(ncol+2))=std::min(RMATRIX(bounds,n,1),RMATRIX(bounds,m,1));
                        // RT differences
                        RMATRIX(merged_tuples,doing,(ncol+3))=std::min(RMATRIX(bounds,n,2),RMATRIX(bounds,m,2));
                        RMATRIX(merged_tuples,doing,(ncol+4))=std::max(RMATRIX(bounds,n,3),RMATRIX(bounds,m,3));
                        // mark the (n-1)-tupels involved
                        RMATRIX(merged_tuples,doing,(siz-2))=(n+1);
                        RMATRIX(merged_tuples,doing,(siz-1))=(m+1);
                        doing++;
                        // add counter
                        RVECTOR(keeper,n)=(RVECTOR(keeper,n)+1);
                        RVECTOR(keeper,m)=(RVECTOR(keeper,m)+1);
                    }
                }
            }
        }

        UNPROTECT(5);
        return(merged_tuples);
    }


}
