#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <vector>
#include <algorithm>
#include <stack>

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RMATRIX2(m,i,j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RRow(m) (INTEGER(GET_DIM(m))[0])
#define RCol(m) (INTEGER(GET_DIM(m))[1])


struct part{
    size_t from;
    size_t to;
    size_t dimension;
    std::vector<double> UB;
    std::vector<double> LB;
    std::vector<int> ended;
};

struct buck{
    size_t counts;
    double distort;
    std::vector<double> UB;
    std::vector<double> LB;
    std::vector<double> centro;
};

double *qua_a;
size_t qua_b=0;
bool smaller (size_t i,size_t j) {return(qua_a[qua_b+i]<qua_a[qua_b+j]);}

extern "C"{


/******************************************************************************/
/* assign buckets *************************************************************/
/******************************************************************************/

    SEXP bucket(
        SEXP quant,
        SEXP sized,
        SEXP pBar
    ){

            PROTECT(quant = AS_NUMERIC(quant));
            qua_a = NUMERIC_POINTER(quant);
            PROTECT(sized = AS_NUMERIC(sized));
            double *sized2;
            sized2 = NUMERIC_POINTER(sized);
            size_t n=0,m=0,d=0,k=0,s=0,ncol=0,nrow=0,counter=0;
            int done;
            double gapsize,min_value,max_value;
            ncol=RCol(quant);
            nrow=RRow(quant);
            if(nrow==1){
                Rprintf("Nothing to partition \n");
                return 0;
            };
            int *here;
            here = new int[nrow];
            for(n=0;n<nrow;n++){here[n]=n;}

            std::stack<part> partits;
            part old_part1; /* store intermediates for lower partition */
            part old_part2; /* store intermediates for higher partition */
            std::stack<buck> buckits;

            /* initialize incl. bounds */
            partits.push(part());
            partits.top().from=0;
            partits.top().to=nrow;
            partits.top().dimension=0;
            for(n=0;n<ncol;n++){
                min_value=RMATRIX(quant,0,n);
                max_value=RMATRIX(quant,0,n);
                for(m=1;m<nrow;m++){
                    if(RMATRIX(quant,m,n)<min_value){
                        min_value=RMATRIX(quant,m,n);
                    }
                    if(RMATRIX(quant,m,n)>max_value){
                        max_value=RMATRIX(quant,m,n);
                    }
                }
                partits.top().LB.push_back(min_value);
                partits.top().UB.push_back(max_value);
                if((max_value-min_value)>*(sized2+n)){
                    partits.top().ended.push_back(0);
                }else{
                    partits.top().ended.push_back(1);
                }
            }

            SEXP utilsPackage; /* definitions for the progres bar */
            PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
            SEXP percentComplete;
            PROTECT(percentComplete = NEW_NUMERIC(1));
            double *rPercentComplete;
            rPercentComplete = NUMERIC_POINTER(percentComplete);

            while(partits.size()>0){

                /* partition *************************************************************************/
                qua_b=(partits.top().dimension*nrow);
                std::sort(&here[partits.top().from],&here[partits.top().to],smaller);
                if((partits.top().to-partits.top().from)>1){ /* >1 data point */
                    m=(partits.top().from+1);
                    if((partits.top().to-partits.top().from)>2){ /* >2 data point */
                        gapsize=fabs(RMATRIX(quant,here[partits.top().from+1],partits.top().dimension)-RMATRIX(quant,here[partits.top().from],partits.top().dimension));
                        for(n=(partits.top().from+2);n<partits.top().to;n++){ /* find best split - largest gap */
                            if(fabs(RMATRIX(quant,here[n],partits.top().dimension)-RMATRIX(quant,here[n-1],partits.top().dimension))>gapsize){
                                gapsize=fabs(RMATRIX(quant,here[n],partits.top().dimension)-RMATRIX(quant,here[n-1],partits.top().dimension));
                                m=n;
                            }
                        }
                    }
                }else{ /* 1 data point */
                    m=(partits.top().from);
                }

                /* remove old partition **************************************************************/
                old_part1=partits.top();
                old_part2=partits.top();
                partits.pop();

                /* new, lower partition **************************************************************/
                old_part1.LB[old_part1.dimension]=RMATRIX(quant,here[old_part1.from],old_part1.dimension);
                if(m>old_part1.from){d=(m-1);}else{d=old_part1.from;};
                old_part1.UB[old_part1.dimension]=RMATRIX(quant,here[d],old_part1.dimension);
                if((RMATRIX(quant,here[d],old_part1.dimension)-RMATRIX(quant,here[old_part1.from],old_part1.dimension))<=*(sized2+old_part1.dimension)){
                    old_part1.ended[old_part1.dimension]=1;
                }
                for(n=0;n<ncol;n++){
                    if(n!=old_part1.dimension){
                        min_value=RMATRIX(quant,here[old_part1.from],n);
                        max_value=RMATRIX(quant,here[old_part1.from],n);
                        for(s=old_part1.from;s<=d;s++){
                            if(RMATRIX(quant,here[s],n)<min_value){
                                min_value=RMATRIX(quant,here[s],n);
                            }
                            if(RMATRIX(quant,here[s],n)>max_value){
                                max_value=RMATRIX(quant,here[s],n);
                            }
                        }
                        old_part1.LB[n]=min_value;
                        old_part1.UB[n]=max_value;
                        if((max_value-min_value)>*(sized2+n)){
                            old_part1.ended[n]=0;
                        }else{
                            old_part1.ended[n]=1;
                        }
                    }
                }
                done=1; /* required bucket sizes reached? */
                if(d!=old_part1.from){ /* only one instance left*/
                    for(n=0;n<ncol;n++){
                        if(old_part1.ended[n]==0){
                            done=2;
                            break;
                        }
                    }
                }
                if(done==2){ /* create new partition */
                    k=old_part1.dimension+1;
                    for(n=0;n<ncol;n++){
                        if(k>=ncol){
                            k=0;
                        }
                        if(old_part1.ended[k]==0){
                            break;
                        }else{
                            k++;
                        }
                    }
                    partits.push(part());
                    partits.top().from=old_part1.from;
                    partits.top().to=m;
                    partits.top().dimension=k;
                    partits.top().LB=old_part1.LB;
                    partits.top().UB=old_part1.UB;
                    partits.top().ended=old_part1.ended;
               }else{ /* create bucket */
                    buckits.push(buck());
                    buckits.top().counts=(d-old_part1.from+1);
                    buckits.top().distort=0;
                    buckits.top().UB=old_part1.UB;
                    buckits.top().LB=old_part1.LB;
                    for(n=0;n<ncol;n++){
                        buckits.top().centro.push_back(buckits.top().LB[n]+((buckits.top().UB[n]-buckits.top().LB[n])/2));
                    }
                    counter=(counter+buckits.top().counts);
                    *rPercentComplete = counter;
                    eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
                }
                /* new, higher partition *************************************************************/
                old_part2.LB[old_part2.dimension]=RMATRIX(quant,here[m],old_part2.dimension);
                old_part2.UB[old_part2.dimension]=RMATRIX(quant,here[(old_part2.to-1)],old_part2.dimension);
                if((RMATRIX(quant,here[(old_part2.to-1)],old_part2.dimension)-RMATRIX(quant,here[m],old_part2.dimension))<=*(sized2+old_part2.dimension)){
                    old_part2.ended[old_part2.dimension]=1;
                }
                for(n=0;n<ncol;n++){
                  /*  if(n!=old_part2.dimension){ */
                        min_value=RMATRIX(quant,here[m],n);
                        max_value=RMATRIX(quant,here[m],n);
                        for(s=m;s<old_part2.to;s++){
                            if(RMATRIX(quant,here[s],n)<min_value){
                                min_value=RMATRIX(quant,here[s],n);
                            }
                            if(RMATRIX(quant,here[s],n)>max_value){
                                max_value=RMATRIX(quant,here[s],n);
                            }
                        }
                        old_part2.LB[n]=min_value;
                        old_part2.UB[n]=max_value;
                        if((max_value-min_value)>*(sized2+n)){
                            old_part2.ended[n]=0;
                        }else{
                            old_part2.ended[n]=1;
                        }
                   /* }*/
                }
                done=1; /* required bucket sizes reached? */
                if((old_part2.to-1)>m){ /* only one instance left */
                    for(n=0;n<ncol;n++){
                        if(old_part2.ended[n]==0){
                            done=2;
                            break;
                        }
                    }
                }
                if(done==2){ /* create new partition */
                    k=old_part2.dimension+1;
                    for(n=0;n<ncol;n++){
                        if(k>=ncol){
                            k=0;
                        }
                        if(old_part2.ended[k]==0){
                            break;
                        }else{
                            k++;
                        }
                    }
                    partits.push(part());
                    partits.top().from=m;
                    partits.top().to=old_part2.to;
                    partits.top().dimension=k;
                    partits.top().LB=old_part2.LB;
                    partits.top().UB=old_part2.UB;
                    partits.top().ended=old_part2.ended;
               }else{ /* create bucket */
                    buckits.push(buck());
                    buckits.top().counts=(old_part2.to-m);
                    buckits.top().distort=0;
                    buckits.top().UB=old_part2.UB;
                    buckits.top().LB=old_part2.LB;
                    for(n=0;n<ncol;n++){
                        buckits.top().centro.push_back(buckits.top().LB[n]+((buckits.top().UB[n]-buckits.top().LB[n])/2));
                    }
                    counter=(counter+buckits.top().counts);
                    *rPercentComplete = counter;
                    eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
                }
            }

            SEXP results;
            k=buckits.size();
            PROTECT(results = allocMatrix(REALSXP, k, ncol*3+1));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<k;n++){
                for(m=0;m<ncol;m++){
                    results2[m*k+n]=buckits.top().centro[m];
                    results2[ncol*k+m*k+n]=buckits.top().LB[m];
                    results2[2*ncol*k+m*k+n]=buckits.top().UB[m];
                }
                results2[3*ncol*k+n]=buckits.top().counts;
                buckits.pop();
            }
            delete[] here;
            UNPROTECT(5);
            return(results);

    }


/******************************************************************************/
/* PNN Algorithm **************************************************************/
/******************************************************************************/

    SEXP PNNA(
        SEXP delmass,
        SEXP intrat,
        SEXP mass,
        SEXP markermass,
        SEXP size_deltamass,
        SEXP size_mass,
        SEXP size_int,
        SEXP scale_deltamass,
        SEXP scale_mass,
        SEXP scale_int
    ){

            PROTECT(delmass = AS_NUMERIC(delmass));
            double *delmass2;
            delmass2 = NUMERIC_POINTER(delmass);
            PROTECT(intrat = AS_NUMERIC(intrat));
            double *intrat2;
            intrat2 = NUMERIC_POINTER(intrat);
            PROTECT(mass = AS_NUMERIC(mass));
            double *mass2;
            mass2 = NUMERIC_POINTER(mass);
            PROTECT(markermass = AS_NUMERIC(markermass));
            double *markermass2;
            markermass2 = NUMERIC_POINTER(markermass);
            PROTECT(size_deltamass = AS_NUMERIC(size_deltamass));
            double *size_deltamass2;
            size_deltamass2 = NUMERIC_POINTER(size_deltamass);
            PROTECT(size_mass = AS_NUMERIC(size_mass));
            double *size_mass2;
            size_mass2 = NUMERIC_POINTER(size_mass);
            PROTECT(size_int = AS_NUMERIC(size_int));
            double *size_int2;
            size_int2 = NUMERIC_POINTER(size_int);
            PROTECT(scale_deltamass = AS_NUMERIC(scale_deltamass));
            double *scale_deltamass2;
            scale_deltamass2 = NUMERIC_POINTER(scale_deltamass);
            PROTECT(scale_mass = AS_NUMERIC(scale_mass));
            double *scale_mass2;
            scale_mass2 = NUMERIC_POINTER(scale_mass);
            PROTECT(scale_int = AS_NUMERIC(scale_int));
            double *scale_int2;
            scale_int2 = NUMERIC_POINTER(scale_int);

            size_t leng = LENGTH(delmass);
            size_t n=0,m=0,k=0,leng2=0,leng3=0;
            int index_dist_to=0,index_at=0,index_temp=0;
            double mindist,mindist_2,mindist_3,sumcount;

            /* initialize cluster container */
            std::vector<double> low_delmass (leng);
            std::vector<double> high_delmass (leng);
            std::vector<double> low_intrat (leng);
            std::vector<double> high_intrat (leng);
            std::vector<double> low_mass (leng);
            std::vector<double> high_mass (leng);
            std::vector<double> low_marker (leng);
            std::vector<double> high_marker (leng);
            std::vector<double> med_delmass (leng);
            std::vector<double> med_mass (leng);
            std::vector<double> med_intrat (leng);
            std::vector<double> med_marker (leng);
            std::vector<double> counts (leng);
            std::vector<int> ID (leng);
            for(n=0;n<leng;n++){
                low_delmass[n]=*(delmass2+n);
                high_delmass[n]=*(delmass2+n);
                low_intrat[n]=*(intrat2+n);
                high_intrat[n]=*(intrat2+n);
                low_mass[n]=*(mass2+n);
                high_mass[n]=*(mass2+n);
                low_marker[n]=*(markermass2+n);
                high_marker[n]=*(markermass2+n);
                med_delmass[n]=*(delmass2+n);
                med_mass[n]=*(mass2+n);
                med_intrat[n]=*(intrat2+n);
                med_marker[n]=*(markermass2+n);
                counts[n]=1;
                ID[n]=n; /* contains shrinking pointers to original data - erasion for each merge or dead end centroid*/
            }

            /* initialize all distances */
            std::vector<double> dist (leng);
            std::vector<int> dist_to (leng);
            mindist= R_PosInf;
            index_at=-1;
            leng2=0;
            k=0;
            leng3=leng;
            for(n=0;n<leng;n++){
                void R_CheckUserInterrupt(void);
                mindist_2= R_PosInf;
                mindist_3= R_PosInf;
                index_dist_to=-1;
                if(n>0){
                    for(m=(n-1);m>=0;m--){
                        if(fabs(*(mass2+n)-*(mass2+m))<*size_mass2){
                            if(fabs(*(intrat2+n)-*(intrat2+m))<*size_int2){
                                if(fabs(*(delmass2+n)-*(delmass2+m))<*size_deltamass2){
                                    if(fabs(*(markermass2+n)-*(markermass2+m))<*size_deltamass2){
                                        mindist_3=((
                                            pow((fabs(*(delmass2+n)-*(delmass2+m))/ *scale_deltamass2),2)+
                                            pow((fabs(*(intrat2+n)-*(intrat2+m))/ *scale_int2),2)+
                                            pow((fabs(*(mass2+n)-*(mass2+m))/ *scale_mass2),2)+
                                            pow((fabs(*(markermass2+n)-*(markermass2+m))/ *scale_deltamass2),2)
                                        )/4);
                                        if(mindist_3<=mindist_2){
                                            mindist_2=mindist_3;
                                            index_dist_to=(signed)m;
                                        };
                                    };
                                };
                            };
                        }else{
                            break;
                        }
                    }
                }
                if(n<leng){
                    for(m=(n+1);m<leng;m++){
                       if(fabs(*(mass2+n)-*(mass2+m))<*size_mass2){
                            if(fabs(*(intrat2+n)-*(intrat2+m))<*size_int2){
                                if(fabs(*(delmass2+n)-*(delmass2+m))<*size_deltamass2){
                                    if(fabs(*(markermass2+n)-*(markermass2+m))<*size_deltamass2){
                                        mindist_3=((
                                            pow((fabs(*(delmass2+n)-*(delmass2+m))/ *scale_deltamass2),2)+
                                            pow((fabs(*(intrat2+n)-*(intrat2+m))/ *scale_int2),2)+
                                            pow((fabs(*(mass2+n)-*(mass2+m))/ *scale_mass2),2)+
                                            pow((fabs(*(markermass2+n)-*(markermass2+m))/ *scale_deltamass2),2)
                                        )/4);
                                        if(mindist_3<=mindist_2){
                                            mindist_2=mindist_3;
                                            index_dist_to=(signed)m;
                                        };
                                    };
                                };
                            };
                        }else{
                            break;
                        }
                    }
                }
                if(index_dist_to!=-1){
                    dist[n]=mindist_2;
                    dist_to[n]=index_dist_to;
                    if(mindist_2<=mindist){
                        mindist=mindist_2;
                        index_at=signed(n-k); /* here: index_at = pointer on ID pointing to data */
                    }
                }else{
                    dist_to[n]=-1;
                    leng2++;
                    ID.erase (ID.begin()+n-k);
                    k++;
                    leng3--;
                }
            }

            /* run */
            while(index_at!=-1){

                void R_CheckUserInterrupt(void);
                /* merge *******************************************************************/
                if(dist_to[dist_to[ID[index_at]]]>-1){
                    sumcount = (counts[dist_to[ID[index_at]]] + counts[ID[index_at]]);
                    med_delmass[dist_to[ID[index_at]]]=(((med_delmass[dist_to[ID[index_at]]]*counts[dist_to[ID[index_at]]])+(med_delmass[ID[index_at]]*counts[ID[index_at]]))/sumcount);
                    med_mass[dist_to[ID[index_at]]]=(((med_mass[dist_to[ID[index_at]]]*counts[dist_to[ID[index_at]]])+(med_mass[ID[index_at]]*counts[ID[index_at]]))/sumcount);
                    med_intrat[dist_to[ID[index_at]]]=(((med_intrat[dist_to[ID[index_at]]]*counts[dist_to[ID[index_at]]])+(med_intrat[ID[index_at]]*counts[ID[index_at]]))/sumcount);
                    med_marker[dist_to[ID[index_at]]]=(((med_marker[dist_to[ID[index_at]]]*counts[dist_to[ID[index_at]]])+(med_marker[ID[index_at]]*counts[ID[index_at]]))/sumcount);
                    low_delmass[dist_to[ID[index_at]]]=std::min(low_delmass[dist_to[ID[index_at]]],low_delmass[ID[index_at]]);
                    high_delmass[dist_to[ID[index_at]]]=std::max(high_delmass[dist_to[ID[index_at]]],high_delmass[ID[index_at]]);
                    low_intrat[dist_to[ID[index_at]]]=std::min(low_intrat[dist_to[ID[index_at]]],low_intrat[ID[index_at]]);
                    high_intrat[dist_to[ID[index_at]]]=std::max(high_intrat[dist_to[ID[index_at]]],high_intrat[ID[index_at]]);
                    low_mass[dist_to[ID[index_at]]]=std::min(low_mass[dist_to[ID[index_at]]],low_mass[ID[index_at]]);
                    high_mass[dist_to[ID[index_at]]]=std::max(high_mass[dist_to[ID[index_at]]],high_mass[ID[index_at]]);
                    low_marker[dist_to[ID[index_at]]]=std::min(low_marker[dist_to[ID[index_at]]],low_marker[ID[index_at]]);
                    high_marker[dist_to[ID[index_at]]]=std::max(high_marker[dist_to[ID[index_at]]],high_marker[ID[index_at]]);
                    counts[dist_to[ID[index_at]]]=sumcount;
                    index_temp=dist_to[ID[index_at]];
                    dist_to[ID[index_at]]=-2;
                    ID.erase (ID.begin()+index_at);
                    leng3--;
                }else{
                    index_temp=ID[index_at];
                }
                /* recalculate distance ****************************************************/
                mindist_2= R_PosInf;
                index_dist_to=-1;
                for(n=0;n<leng3;n++){
                    if(dist_to[ID[n]]>-1){
                        if(ID[n]!=index_temp){
                            if(fabs(*(mass2+ID[n])-*(mass2+index_temp))<*size_mass2){
                                if(fabs(*(intrat2+ID[n])-*(intrat2+index_temp))<*size_int2){
                                    if(fabs(*(delmass2+ID[n])-*(delmass2+index_temp))<*size_deltamass2){
                                        if(fabs(*(markermass2+ID[n])-*(markermass2+index_temp))<*size_deltamass2){
                                            mindist_3=((
                                                pow((fabs(*(delmass2+ID[n])-*(delmass2+index_temp))/ *scale_deltamass2),2)+
                                                pow((fabs(*(intrat2+ID[n])-*(intrat2+index_temp))/ *scale_int2),2)+
                                                pow((fabs(*(mass2+ID[n])-*(mass2+index_temp))/ *scale_mass2),2)+
                                                pow((fabs(*(markermass2+ID[n])-*(markermass2+index_temp))/ *scale_deltamass2),2)
                                            )/4);
                                            if(mindist_3<=mindist_2){
                                                mindist_2=mindist_3;
                                                index_dist_to=ID[n];
                                            };
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(index_dist_to!=-1){
                    dist[index_temp]=mindist_2;
                    dist_to[index_temp]=index_dist_to;
                }else{
                    dist_to[index_temp]=-1;
                    leng2++;
                    for(n=0;n<leng3;n++){
                        if(ID[n]==index_temp){
                            ID.erase (ID.begin()+n);
                            break;
                        }
                    }
                    leng3--;
                }
                /* find minimum distance ****************************************************/
                index_at=-1;
                mindist=R_PosInf;
                for(n=0;n<leng3;n++){
                    if(dist_to[ID[n]]>-1){
                        if(dist[ID[n]]<mindist){
                            mindist=dist[ID[n]];
                            index_at=n;
                        }
                    }
                }

            }

            SEXP quantiz;
            PROTECT(quantiz = allocMatrix(REALSXP, leng2, 13));
            double *quantiz2;
            quantiz2 = REAL(quantiz);
            m=0;
            for(n=0;n<leng;n++){
                if(dist_to[n]==-1){
                    quantiz2[m]=low_delmass[n];
                    quantiz2[(1*leng2)+m]=high_delmass[n];
                    quantiz2[(2*leng2)+m]=low_mass[n];
                    quantiz2[(3*leng2)+m]=high_mass[n];
                    quantiz2[(4*leng2)+m]=low_intrat[n];
                    quantiz2[(5*leng2)+m]=high_intrat[n];
                    quantiz2[(6*leng2)+m]=low_marker[n];
                    quantiz2[(7*leng2)+m]=high_marker[n];
                    quantiz2[(8*leng2)+m]=med_delmass[n];
                    quantiz2[(9*leng2)+m]=med_mass[n];
                    quantiz2[(10*leng2)+m]=med_intrat[n];
                    quantiz2[(11*leng2)+m]=med_marker[n];
                    quantiz2[(12*leng2)+m]=counts[n];
                    m++;
                }
            }

            UNPROTECT(11);
            return(quantiz);

    }

}

