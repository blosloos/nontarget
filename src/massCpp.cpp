#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R.h>

using namespace std;

struct dat{
       vector<double> mass;
       vector<double> retent;
       bool operator() (size_t i, size_t j) { return (mass[i]<mass[j]);}
       };

extern "C"
{
void mass(  double *mass,
            double *retent,
            int *a,
            double *masstol,
            double *massfrac,
            double *rttollow,
            double *rttolup,
            int *manyisos,
            double *isomat1,
            int *isomat3,
            double *maxmass1,
            int *isomat4,
            int *entry,
            int *ppm2,
            int *getit1,
            int *getit2,
            int *getit4,
            int *getit5,
            int *getit6
        )
    {

    size_t i=0,k=0,j=0,do_at=0,l=0,upcount=0,lowcount=0,rtup=0,rtlow=0,howmany=0;
    int entry2=*entry;
    double uptol, lowtol, thismasslow, thismasslow2, thismassup, thismassup2;
    //generate index vectors:
    vector<int> index;
    vector<int> index2;

    SEXP getit1bstore;
    PROTECT(getit1bstore = NEW_INTEGER(*a+1));
    int *getit1b;
    getit1b = INTEGER_POINTER(getit1bstore);
    for(k=0;k<(*a+1);k++){*(getit1b+k) = 0;}

    SEXP getit2bstore;
    PROTECT(getit2bstore = NEW_INTEGER(*a+1));
    int *getit2b;
    getit2b = INTEGER_POINTER(getit2bstore);
    for(k=0;k<(*a+1);k++){*(getit2b+k) = 0;}

    SEXP getit4bstore;
    PROTECT(getit4bstore = NEW_INTEGER(*a+1));
    int *getit4b;
    getit4b = INTEGER_POINTER(getit4bstore);
    for(k=0;k<(*a+1);k++){*(getit4b+k) = 0;}

    SEXP getit5bstore;
    PROTECT(getit5bstore = NEW_INTEGER(*a+1));
    int *getit5b;
    getit5b = INTEGER_POINTER(getit5bstore);
    for(k=0;k<(*a+1);k++){*(getit5b+k) = 0;}

    SEXP getit6bstore;
    PROTECT(getit6bstore = NEW_INTEGER(*a+1));
    int *getit6b;
    getit6b = INTEGER_POINTER(getit6bstore);
    for(k=0;k<(*a+1);k++){*(getit6b+k) = 0;}

    // read in data: /////////////////////////////////////////////////////////
    dat dat1;
    // (a) m/z
    for(i=0;i<(unsigned)*a;i++){dat1.mass.push_back(mass[i]);}
    // (b) retention time
    for(i=0;i<(unsigned)*a;i++){dat1.retent.push_back(retent[i]);}
    //////////////////////////////////////////////////////////////////////////

    // run search: ///////////////////////////////////////////////////////////
    for(i=0;i<(unsigned)*a;i++){
        if((i==0)||(dat1.retent[i]!=dat1.retent[i-1])){ // build new index vector = all peaks within dRT
            uptol=dat1.retent[i]+*rttolup;
            lowtol=dat1.retent[i]+*rttollow;
            while(((rtup+1)<(unsigned)*a) && (dat1.retent[rtup+1]<=uptol)) {rtup++;};
            while(((rtlow+1)<(unsigned)*a) && (dat1.retent[rtlow]<lowtol)) {rtlow++;};
            index.erase (index.begin(),index.end());
            for(j=rtlow; j<=rtup; j++){
                    index.push_back(j);
            };
            sort(index.begin(),index.end(),dat1); // sort vector by mass
            do_at=0;
            while( (((do_at+1)<index.size()) && (dat1.mass[index[do_at]]<=dat1.mass[i]))) {do_at++;};
        }else{ // old index vector valid = sorted by mass, searched by increasing masses only
            if(dat1.mass[i]<dat1.mass[i-1]){ // reset
                do_at=0;
                while ((((do_at+1) < index.size() ) && (dat1.mass[index[do_at]] <= dat1.mass[i]))) {do_at++;};
            }else{ // increase
                 while ((((do_at+1) < index.size() ) && (dat1.mass[index[do_at]] <= dat1.mass[i]))) {do_at++;};
            }
        }; // if1



        howmany = index.size(); // within RT-window?
        if(howmany>0){ // sth in vector and not at its mass end
            if(*ppm2==1){
                thismasslow=(dat1.mass[i]-(dat1.mass[i]**masstol/1E6));
                thismasslow2=(dat1.mass[i]-(((dat1.mass[i]**masstol/1E6))**massfrac));
                thismassup=(dat1.mass[i]+(dat1.mass[i]**masstol/1E6));
                thismassup2=(dat1.mass[i]+(((dat1.mass[i]**masstol/1E6))**massfrac));
            }else{
                thismasslow=(dat1.mass[i]-*masstol);
                thismasslow2=(dat1.mass[i]-(*masstol**massfrac));
                thismassup=(dat1.mass[i]+*masstol);
                thismassup2=(dat1.mass[i]+(*masstol**massfrac));
            };
            upcount=do_at;
            lowcount=do_at;
            for(k = 0; k<(unsigned)*manyisos; k++){
                while(((upcount+1)<howmany ) && (dat1.mass[index[upcount+1]]<=(thismassup + isomat1[k]))) {
                    upcount++;
                };
                while(((lowcount+1)<howmany) && (dat1.mass[index[lowcount]]<(thismasslow + isomat1[k]))) {
                    lowcount++;
                };//set mass window
                for(l=lowcount;l<=upcount;l++){
                    if( (dat1.mass[index[l]]<=(thismassup + isomat1[k])) && (dat1.mass[index[l]]>=(thismasslow + isomat1[k]))){
                        if(*(getit2b+index[l])<(*entry+1)){ // from?
                            getit2[(index[l]**entry)+(*(getit2b+index[l]))]=(i+1);
                            *(getit2b+index[l]) = (*(getit2b+index[l])+1);
                        }
                        if( *(getit4b+i)<(*entry+1) ){ // to?
                            getit4[i**entry + *(getit4b+i)]=(index[l]+1);
                            *(getit4b+i) = (*(getit4b+i)+1);
                        }
                        if(*(getit1b+i)<(*entry+1)){ // which isotope?
                            getit1[i**entry+*(getit1b+i)]=(k+1);
                            *(getit1b+i) = (*(getit1b+i)+1);
                        }
                        if(*(getit6b+i)<(*entry+1)){ // which charge level?
                            getit6[i**entry+*(getit6b+i)]=(isomat4[k]);
                            *(getit6b+i) = (*(getit6b+i)+1);
                        }
                        // large or small mass tolerance?
                        if( (dat1.mass[index[l]]<=(thismassup2 + isomat1[k])) && (dat1.mass[index[l]]>=(thismasslow2 + isomat1[k]))){
                            getit5[i**entry+*(getit5b+i)]=1; // 1 = small
                            *(getit5b+i) = (*(getit5b+i)+1);
                        }else{
                            getit5[i**entry+*(getit5b+i)]=1; // 2 = large
                            *(getit5b+i) = (*(getit5b+i)+1);
                        };
                        isomat3[k]=isomat3[k]+1;
                    } // if within mass-window
                } // for l
            } // for k
        }  // if howmany>0
    } // for i
    //////////////////////////////////////////////////////////////////////////

    // check if entry has reached limit //////////////////////////////////////
    for(j=0;j<(*a-1);j++){if(*(getit1b+j)>entry2){*entry=*(getit1b+j);};};
    for(j=0;j<(*a-1);j++){if(*(getit2b+j)>entry2){*entry=*(getit2b+j);};};
    for(j=0;j<(*a-1);j++){if(*(getit4b+j)>entry2){*entry=*(getit4b+j);};};
    for(j=0;j<(*a-1);j++){if(*(getit5b+j)>entry2){*entry=*(getit5b+j);};};

    UNPROTECT(5);

    } // main
} // extern "C"









