#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP R2_C(SEXP RinMatrix, SEXP EINSEN){

SEXP ret = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;

Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

double value1;
double value2;

Rvalue           = coerceVector(RinMatrix, REALSXP);
double *Rval     = REAL(Rvalue);
double *einsen   = REAL(EINSEN);

double  ones;
double  zeros;
double     id;
double     id2;
double     count;

double  freqsite1;
double  freqsite2;
double  freqsite1_low;
double  freqsite2_low;
double  site_length = I;
double  d_raw;
double  r2;
int     fill;

PROTECT(ret = allocVector(REALSXP,J*J));

// Init ret

for(int i=0; i< J*J; i++){
REAL(ret)[i]=0; // default monomorph
}

//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}


fill = 0;
for (int m=0; m < J-1; m++){
 
 // sites 1
    ones  = einsen[m];
    zeros = site_length - ones;
    if(ones>=zeros){
       freqsite1 = ones/site_length;
       id = 1;
    }else{
       freqsite1 = zeros/site_length;
       id = 0;
    } 
 // End of sites 1
 
 for (int i = m+1; i < J; i++){

 // sites 2
    ones  = einsen[m];
    zeros = I - ones;
    if(ones>=zeros){
       freqsite2 = ones/site_length;
       id2 = 1;
    }else{
       freqsite2 = zeros/site_length;
       id2 = 0;
    } 
 // End of sites 2
    
   count = 0;
   for (int j = 0; j < I; j++){
   
    value1 = Rval[j +I*m];
    value2 = Rval[j +I*i];
    
     if(value1==id && value2 ==id2){
       count ++;
     }
    
   }

   d_raw  = (double)count/(double)site_length - freqsite1*freqsite2;

   freqsite1_low = 1-freqsite1;
   freqsite2_low = 1-freqsite2;
   r2 = (d_raw*d_raw)/(freqsite1*freqsite1_low*freqsite2*freqsite2_low);
   REAL(ret)[fill] = r2;
   fill ++;

   //printf("%f",value1);
  
   // if(value2==6){break;}
   // if(value2==5){break;}  
   // if(value1!=value2){
        //printf("Treffer");
   //	INTEGER(ret)[i]=1; // polymorphic
   //     break;      
     
   
  
 }
}

//

UNPROTECT(1);

return ret;

}




/*SEXP Ccompare(SEXP Rvector1, SEXP Rvector2)
{

SEXP val = R_NilValue;
PROTECT(val= Rf_allocVector(INTSXP,1));
INTEGER(val)[0] = 1;

// Rvector1 = coerceVector(Rvector1, INTSXP);
// Rvector2 = coerceVector(Rvector2, INTSXP);

int size;
const double *vec1 = REAL(Rvector1);
const double *vec2 = REAL(Rvector2);
size            = length(Rvector1);

for(int i=0; i < size; i++) {

	if(vec1[i]!=vec2[i]){          
	   INTEGER(val)[0] = 0;
           break;
	}
}

UNPROTECT(1);
return val;


}
*/
