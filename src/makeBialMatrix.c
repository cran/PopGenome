#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>

SEXP makeBialMatrix(SEXP RinMatrix){

SEXP list = R_NilValue;

int I;
int J;
SEXP Rdim;
SEXP Rvalue;
SEXP BACK;
SEXP TRANS;
SEXP SUBST;


Rdim = getAttrib(RinMatrix, R_DimSymbol);
I    = INTEGER(Rdim)[0]; // Reihen 
J    = INTEGER(Rdim)[1]; // Spalten

double value1=-1;
double value2=-1;
double count1;
double count2;
double minor=0;
double major=0;

PROTECT(BACK = allocMatrix(INTSXP, I, J));

Rvalue           = coerceVector(RinMatrix, INTSXP);
int *Rval        = INTEGER(Rvalue);

PROTECT(TRANS = allocVector(INTSXP,J));
PROTECT(SUBST = allocMatrix(INTSXP,2,J));

// Init SUBST
for (int i = 0; i < J; i++){
 for (int j = 0; j < 2; j++){
   INTEGER(SUBST)[j +2*i]=0;
 }
}

// Init TRANS
for(int i=0; i< J; i++){
 INTEGER(TRANS)[i]=0; // default monomorph
}

// Init Matrix
for (int i = 0; i < J; i++){
 for (int j = 0; j < I; j++){
   INTEGER(BACK)[j +I*i]=0;
 }
}


//for(int i=0; i< J*I; i++){
//value2 = Rval[i];
//printf("%f",value2);
//}

for (int i = 0; i < J; i++){
   // value1 = Rval[I*i];
 for (int j = 0; j < I; j++){
  
  if(j==0){
   value1 = Rval[j +I*i];
   count1 = 0;
  for (int k = 0; k < I; k++){
      if(Rval[k + I*i]==value1){
        count1++;
      }else{
        value2 = Rval[k + I*i];
      }        
  }
   count2 = 0;
  for (int m = 0; m < I; m++){
      if(Rval[m + I*i]==value2){
        count2++; 
      } 
  }
 // } // end if j==0

// Transitions
  if((value1==1 && value2==2) || (value1==2 && value2==1)){INTEGER(TRANS)[i]=1;}
  if((value1==3 && value2==4) || (value1==4 && value2==3)){INTEGER(TRANS)[i]=1;}


 
  if(count1<=count2){
    minor = value1;
    major = value2;
  }else{
    minor = value2;
    major = value1;	
  }
  
 //fill SUBST
 INTEGER(SUBST)[0+2*i]= (int)minor;
 INTEGER(SUBST)[1+2*i]= (int)major;

 } // end if j==0
  
  if(minor==Rval[j +I*i]){
    INTEGER(BACK)[j + I*i] = 1;
  }

 
  //if(major==Rval[j +I*i]){
  //  REAL(BACK)[j + I*i] = 0.0;
  //}

  



 } 
}


   //printf("%f",value1);
  
//   if(value2==6){break;}
//   if(value2==5){break;}  
//   if(value1!=value2){
        //printf("Treffer");
//	INTEGER(ret)[i]=1; // polymorphic
//        break;      
//   } 
   
// }
//}


    // Creating a list with 3 vector elements:
   PROTECT(list = allocVector(VECSXP, 3)); 
     // attaching myint vector to list:
   SET_VECTOR_ELT(list, 0, BACK); 
     // attaching mydouble vector to list:
   SET_VECTOR_ELT(list, 1, TRANS); 
     // and attaching the vector names:
     //setAttrib(list, R_NamesSymbol, list_names); 
   SET_VECTOR_ELT(list, 2, SUBST);
   UNPROTECT(4);
   return list;

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
