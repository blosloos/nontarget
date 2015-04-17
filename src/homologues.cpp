#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <vector>

extern "C"{

/******************************************************************************/
/* filter homologue peak-to-peak differences **********************************/
/******************************************************************************/

    SEXP homol_1(
        SEXP storage,
        SEXP mass,
        SEXP RT,
        SEXP ord,
        SEXP lowbound,
        SEXP highbound,
        SEXP minrt,
        SEXP maxrt,
        SEXP delmz,
        SEXP ppm,
        SEXP rttol,
        SEXP minleng,
        SEXP inter, // run in interactive mode?
        SEXP pBar
    ){

            PROTECT(storage = AS_INTEGER(storage));
            size_t stor;
            stor = INTEGER_VALUE(storage);
            PROTECT(mass = AS_NUMERIC(mass));
            double *mass2;
            mass2 = NUMERIC_POINTER(mass);
            PROTECT(RT = AS_NUMERIC(RT));
            double *RT2;
            RT2 = NUMERIC_POINTER(RT);
            PROTECT(ord = AS_INTEGER(ord));
            int *ord2;
            ord2 = INTEGER_POINTER(ord);
            PROTECT(lowbound = AS_NUMERIC(lowbound));
            double *lowbound2;
            lowbound2 = NUMERIC_POINTER(lowbound);
            PROTECT(highbound = AS_NUMERIC(highbound));
            double *highbound2;
            highbound2 = NUMERIC_POINTER(highbound);
            PROTECT(minrt = AS_NUMERIC(minrt));
            double *minrt2;
            minrt2 = NUMERIC_POINTER(minrt);
            PROTECT(maxrt = AS_NUMERIC(maxrt));
            double *maxrt2;
            maxrt2 = NUMERIC_POINTER(maxrt);
            PROTECT(delmz = AS_NUMERIC(delmz));
            double *delmz2;
            delmz2 = NUMERIC_POINTER(delmz);
            PROTECT(ppm= AS_INTEGER(ppm));
            int *ppm2;
            ppm2 = INTEGER_POINTER(ppm);
            PROTECT(rttol = AS_NUMERIC(rttol));
            double *rttol2;
            rttol2 = NUMERIC_POINTER(rttol);
            PROTECT(minleng = AS_INTEGER(minleng));
            size_t minleng2;
            minleng2 = INTEGER_VALUE(minleng);

            PROTECT(inter = AS_NUMERIC(inter));
            int intera = INTEGER_VALUE(inter);
            SEXP utilsPackage; /* definitions for progress bar */
            PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
            SEXP percentComplete;
            PROTECT(percentComplete = NEW_NUMERIC(1));
            double *rPercentComplete;
            rPercentComplete = NUMERIC_POINTER(percentComplete);

            size_t k,m,n,x,z,homnum,countit,a,b,c,d,g,anew;
            double delmz3;
            unsigned int leng = LENGTH(mass);
            unsigned int lengbound = LENGTH(highbound);

            SEXP from;
            PROTECT(from = NEW_INTEGER(stor));
            int *atfrom;
            atfrom = INTEGER_POINTER(from);
            for(n=0;n<stor;n++){*(atfrom+n) = 0;}
            SEXP to;
            PROTECT(to = NEW_INTEGER(stor));
            int *atto;
            atto = INTEGER_POINTER(to);
            for(n=0;n<stor;n++){*(atto+n) = 0;};
            SEXP dmass;
            PROTECT(dmass = NEW_NUMERIC(stor));
            double *atdmass;
            atdmass = NUMERIC_POINTER(dmass);
            for(n=0;n<stor;n++){*(atdmass+n) = 0;};
            SEXP usemass;
            PROTECT(usemass = NEW_NUMERIC(stor));
            double *atusemass;
            atusemass = NUMERIC_POINTER(usemass);
            for(n=0;n<stor;n++){*(atusemass+n) = 0;};
            SEXP dRT;
            PROTECT(dRT = NEW_NUMERIC(stor));
            double *atdRT;
            atdRT = NUMERIC_POINTER(dRT);
            for(n=0;n<stor;n++){*(atdRT+n) = 0;};
            SEXP ordit;
            PROTECT(ordit = NEW_INTEGER(stor));
            int *atordit;
            atordit = INTEGER_POINTER(ordit);
            for(n=0;n<stor;n++){*(atordit+n) = n;};
            SEXP fertig;
            PROTECT(fertig = NEW_INTEGER(stor));
            int *fer;
            fer = INTEGER_POINTER(fertig);
            for(n=0;n<stor;n++){*(fer+n) = 0;};
            SEXP atk;
            PROTECT(atk = NEW_INTEGER(leng));
            int *atk2;
            atk2 = INTEGER_POINTER(atk);
            for(n=0;n<leng;n++){*(atk2+n) = n;};

            /* intermediate output vectors */
            std::vector<double> out_from;
            std::vector<double> out_to;
            std::vector<double> out_dmass;
            std::vector<double> out_dRT;
            std::vector<double> out_HS;

            /* HS screening vectors */
            std::vector<double> do_from;
            std::vector<double> do_to;
            std::vector<double> do_dmz;
            std::vector<double> do_tolmz;
            std::vector<double> do_dRT;
            std::vector<unsigned int> do_index;

            x=0;
            homnum=1;
            /* (1) over all mass windows ... */
            for( z=0; z<lengbound; z++ ){
                // Rprintf("*");
                if(intera==1){
                    *rPercentComplete = z;
                    eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
                }
                m=0;
                void R_CheckUserInterrupt(void);
                /* (2) ... find relevant mass differences ... */
                for(n=0;n<(leng-1);n++){
                    if( m >= ((stor)-1) ){ Rprintf("increase vec_size!"); break; }
                    for(k=*(atk2+n);k<leng;k++){
                        if( *(mass2+k) > ( *(mass2+n) + *(highbound2+z) )){
                            *(atk2+n)=k;
                            break;
                        }
                        if( m >= ((stor)-1) ){ Rprintf("increase vec_size!"); break; }
                        if( (*(RT2+k)>=(*(RT2+n)+*minrt2)) && (*(RT2+k)<=(*(RT2+n)+*maxrt2)) ){
                            if(
                               ((*(mass2+k)-*(mass2+n))<=*(highbound2+z))&&
                               ((*(mass2+k)-*(mass2+n))>=*(lowbound2+z))
                            ){
                                atfrom[m]=*(ord2+n);
                                atto[m]=*(ord2+k);
                                atusemass[m]=((*(mass2+k) + *(mass2+n))/2);
                                atdmass[m]=(*(mass2+k) - *(mass2+n));
                                atdRT[m]=(*(RT2+k) - *(RT2+n));
                                m++;
                            }
                        }
                    }
                }
                /* (3) ... sort index ... */
                R_orderVector(atordit,m,Rf_lang1(dmass),FALSE,FALSE);
                /* (4) ... and screen those mass differences: */
                void R_CheckUserInterrupt(void);
                for(n=0;n<stor;n++){*(fer+n) = 0;};
                for(n=0;n<m;n++){
                    if( *(fer+n)==0 ){
                        *(fer+n) = 1;
                        do_from.clear();
                        do_to.clear();
                        do_dmz.clear();
                        do_tolmz.clear();
                        do_dRT.clear();
                        do_index.clear();
                        do_from.push_back(*(atfrom+*(atordit+n)));
                        do_to.push_back(*(atto+*(atordit+n)));
                        do_dmz.push_back(*(atdmass+*(atordit+n)));
                        if(*ppm2==1){
                            delmz3=(2*(*delmz2 * *(atusemass+*(atordit+n)) /1E6));
                            do_tolmz.push_back(delmz3);
                        }else{
                            delmz3=(2*(*delmz2));
                            do_tolmz.push_back(delmz3);
                        }
                        do_dRT.push_back(*(atdRT+*(atordit+n)));
                        do_index.push_back(n);
                        countit=1;
                        a=0;
                        b=1;
                        anew=1;
                        while( anew==1 ){
                            anew=0;
                            g=b;
                            for(c=a;c<b;c++){
                                /* search forward **********************************************************************/
                                if(do_index[c]<(m-1)){
                                    for(d=(do_index[c]+1);d<m;d++){
                                        if(fabs(*(atdmass+*(atordit+d)) - do_dmz[c]) > do_tolmz[c] ){
                                            break;
                                        }
                                        if( *(fer+d)==0 ){
                                            if( (do_to[c]==*(atfrom+*(atordit+d)) ) || (do_from[c]==*(atto+*(atordit+d))) ){
                                                if( fabs(*(atdRT+*(atordit+d)) - do_dRT[c]) <= *rttol2 ){               /* shift in RT ok?*/
                                                    if( fabs(*(atdmass+*(atordit+d)) - do_dmz[c]) <= do_tolmz[c] ){     /* within mz-tolerance? */
                                                        do_from.push_back(*(atfrom+*(atordit+d)));
                                                        do_to.push_back(*(atto+*(atordit+d)));
                                                        do_dmz.push_back(*(atdmass+*(atordit+d)));
                                                        if(*ppm2==1){
                                                            delmz3=(2*(*delmz2 * *(atusemass+*(atordit+d)) /1E6));
                                                            do_tolmz.push_back(delmz3);
                                                        }else{
                                                            delmz3=(2*(*delmz2));
                                                            do_tolmz.push_back(delmz3);
                                                        }
                                                        do_dRT.push_back(*(atdRT+*(atordit+d)));
                                                        do_index.push_back(d);
                                                        *(fer+d)=1;
                                                        countit++;
                                                        anew=1;
                                                        g++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                /* search backward **********************************************************************/
                                if(do_index[c]>n){
                                    for(d=(do_index[c]);d>0;d--){
                                        if(fabs(*(atdmass+*(atordit+d-1)) - do_dmz[c]) > do_tolmz[c] ){
                                            break;
                                        }
                                        if( *(fer+d-1)==0 ){
                                            if( (do_to[c]==*(atfrom+*(atordit+d-1)) ) || (do_from[c]==*(atto+*(atordit+d-1))) ){
                                                if( fabs(*(atdRT+*(atordit+d-1)) - do_dRT[c]) <= *rttol2 ){               /* shift in RT ok?*/
                                                    if( fabs(*(atdmass+*(atordit+d-1)) - do_dmz[c]) <= do_tolmz[c] ){     /* within mz-tolerance? */
                                                        do_from.push_back(*(atfrom+*(atordit+d-1)));
                                                        do_to.push_back(*(atto+*(atordit+d-1)));
                                                        do_dmz.push_back(*(atdmass+*(atordit+d-1)));
                                                        if(*ppm2==1){
                                                            delmz3=(2*(*delmz2 * *(atusemass+*(atordit+d-1)) /1E6));
                                                            do_tolmz.push_back(delmz3);
                                                        }else{
                                                            delmz3=(2*(*delmz2));
                                                            do_tolmz.push_back(delmz3);
                                                        }
                                                        do_dRT.push_back(*(atdRT+*(atordit+d-1)));
                                                        do_index.push_back(d-1);
                                                        *(fer+d-1)=1;
                                                        countit++;
                                                        anew=1;
                                                        g++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            a=b;
                            b=g;
                        }
                        if(countit>=minleng2){
                            for(c=0;c<countit;c++){
                                out_from.push_back (do_from[c]);
                                out_to.push_back (do_to[c]);
                                out_dmass.push_back (do_dmz[c]);
                                out_dRT.push_back (do_dRT[c]);
                                out_HS.push_back (homnum);
                                x++;
                            }
                            homnum++;
                        }
                    }
                } /*(4)*/

            } /*(1)*/

            SEXP homologues;
            PROTECT(homologues = allocMatrix(REALSXP, x, 5));
            double *hom;
            hom = REAL(homologues);
            for(m=0;m<x;m++){
                hom[m]=out_from[m];
                hom[x+m]=out_to[m];
                hom[(2*x)+m]=out_dmass[m];
                hom[(3*x)+m]=out_dRT[m];
                hom[(4*x)+m]=out_HS[m];
            }

            UNPROTECT(24);
            return(homologues);

    }

}
