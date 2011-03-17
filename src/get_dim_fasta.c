#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <string.h>

SEXP get_dim_fasta (SEXP RRfilename) {

  SEXP ret = R_NilValue;
  
  PROTECT(ret = allocVector(INTSXP,2));

  char *file;
  int  ch;
  int  n_ind=0;

  SEXP Rfilename;
  Rfilename   = STRING_ELT(RRfilename,0);
  file        = (char*)CHAR(Rfilename);


  FILE *fp;
  fp = fopen( file ,"r");

  if(fp==NULL) {
    Rprintf("Cannot open file.\n");
    UNPROTECT(1);
    return R_NilValue;
  }


  // check rows of ff matrix
  // ##################################
  while((ch = fgetc( fp )) != EOF) {
    
   // printf("%c", ch);

    if(ch=='>'){
      n_ind ++;
      
    }

  }

   INTEGER(ret)[0] = n_ind;
  // ##################################


// SAVING NAMES
// ######################################

 rewind(fp);
 SEXP row_names;

 PROTECT(row_names = allocVector(STRSXP,n_ind));
 
  char str[100];
  char temp[100];

  int count;
  int count_ind=-1;

  while((ch = fgetc( fp )) != EOF) {
    
    if(ch=='>'){
         count_ind ++;
         count = 0;
         while(1){
            ch = fgetc( fp );              
           if(ch=='\n'|| ch=='\t'){            
              break;
           }else{
             temp[count] = ch;
             count ++;
             strcat(str,temp);
             SET_STRING_ELT(row_names,count_ind,mkChar(str));
             strcpy(str, "");
           }
         }
        
    }

  }
//##################################################



  

// check columns of ff matrix
// ############################################
  rewind(fp); 
  int n_nucs  = 0;
 
   // Get first nucleotide
    
   while(1){
     ch = fgetc( fp );
     //printf("%c", ch);
     if(ch=='\n'|| ch=='\t'){
        break;
     }
   }
   
  // Count nucleotides
  
   while(1){
   ch = fgetc( fp );

     if(ch=='>' || ch==EOF){
      break;
     }
 
     if(ch!='\n'){      
        n_nucs ++;
     }
    
  // printf("%c", ch);
   } 

  INTEGER(ret)[1] = n_nucs; 
  fclose(fp);
//###################################################

SEXP list;
PROTECT(list = allocVector(VECSXP, 2));
SET_VECTOR_ELT(list, 0, ret);    
SET_VECTOR_ELT(list, 1, row_names); 

UNPROTECT(3);
return list;

}


