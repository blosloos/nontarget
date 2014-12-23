#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <Rdefines.h>
#include <R.h>

using namespace std;

struct dat{
       vector<double> mass;
       vector<double> retent;
       bool operator() (size_t i, size_t j) { return (mass[i]<mass[j]);}
       };

extern "C"
{
void mass(             double *mass, double *retent, int *a, double *masstol, double *massfrac, double *rttollow, double *rttolup,
                       int *manyisos, double *isomat1, int *isomat3, double *maxmass1, int *isomat4, int *entry, int *ppm2,
                       int *getit1, int *getit2, int *getit4, int *getit5, int *getit6
                       )
     {

      size_t i=0,k=0,j=0,m=0,l=0,upcount=0,lowcount=0,rtup=0,rtlow=0,howmany=0;
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
      for(i=0;i<(unsigned)*a;i++)
      {
       if((i==0) || (dat1.retent[i]!=dat1.retent[i-1]))
           {
            uptol=dat1.retent[i]+*rttolup;
            lowtol=dat1.retent[i]+*rttollow;
            while((rtup <(unsigned)*a - 1) && (dat1.retent[rtup]<=uptol)) {rtup=(rtup+1);};
            rtup=(rtup-1);
            while((rtlow <(unsigned)*a - 1) && (dat1.retent[rtlow]<lowtol)) {rtlow=(rtlow+1);};
            index.erase (index.begin(),index.end());
            for(j=rtlow; (j<=rtup) && (j<(unsigned)*a); j++){
                         if(dat1.mass[j]>dat1.mass[i]){
                                          index.push_back(j);
                         };
            };
            sort(index.begin(),index.end(),dat1);
      }else { // if1
            j=0;
            if (index.size()) {
                while (((j < index.size() -1 ) && (dat1.mass[index[j]] <= dat1.mass[i]))){
                     j++;
                };
                index.erase(index.begin(), index.begin() + j);
            }
      }; // if1

      howmany = index.size();

      if(howmany>0){
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

                       upcount=0;
                       lowcount=0;

                       for(k = 0; k<(unsigned)*manyisos; k++){

                             if(upcount<=(howmany-1)){

                             while((upcount < howmany - 1)
                                   && (dat1.mass[index[upcount]] <= (thismassup + isomat1[k]))) {
                                 upcount = upcount + 1;
                             };
                             while((lowcount < howmany - 1)
                                   && (dat1.mass[index[lowcount]]<(thismasslow + isomat1[k]))) {
                                  lowcount =  lowcount + 1;
                             };

                           if(lowcount==howmany){
                              if( (dat1.mass[index[lowcount]]<=(thismassup + isomat1[k])) &&
                                 (dat1.mass[index[lowcount]]>=(thismasslow + isomat1[k]))
                             ){
                                upcount=upcount+1;
                              }
                              else {
                                upcount=lowcount;
                              }
                             };

                             if((upcount-lowcount)>0.0){
                               index2.erase (index2.begin(),index2.end());
                               for(l=lowcount;l<upcount;l++){index2.push_back(index[l]);};
                               for(m=0;m<index2.size();m++){

                                if(*(getit2b+index2[m])<(*entry+1)){              // from?
                                    getit2[(index2[m]**entry)+(*(getit2b+index2[m]))]=(i+1);
                                    *(getit2b+index2[m]) = (*(getit2b+index2[m])+1);
                                }

                                if( *(getit4b+i)<(*entry+1) ){                               // to?
                                    getit4[i**entry + *(getit4b+i)]=(index2[m]+1);
                                    *(getit4b+i) = (*(getit4b+i)+1);
                                }

                                if(*(getit1b+i)<(*entry+1)){                              // which isotope?
                                    getit1[i**entry+*(getit1b+i)]=(k+1);
                                    *(getit1b+i) = (*(getit1b+i)+1);
                                }

                                if(*(getit6b+i)<(*entry+1)){                              // which charge level?
                                    getit6[i**entry+*(getit6b+i)]=(isomat4[k]);
                                    *(getit6b+i) = (*(getit6b+i)+1);
                                }

                               // in which mass tolerance?
                               // large or small?
                               if((dat1.mass[index2[m]]<=(isomat1[k]+thismassup2)) && (dat1.mass[index2[m]]>=(isomat1[k]+thismasslow2))){

                                     getit5[i**entry+*(getit5b+i)]=1;                       // 1 = small
                                     *(getit5b+i) = (*(getit5b+i)+1);

                                }else{

                                     getit5[i**entry+*(getit5b+i)]=1;                       // 2 = large
                                     *(getit5b+i) = (*(getit5b+i)+1);

                                };

                              isomat3[k]=isomat3[k]+1;

                               } // for m
                              if(upcount>howmany){upcount=lowcount;};
                             } // if
                            } // if
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









