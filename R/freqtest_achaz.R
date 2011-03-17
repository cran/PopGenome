# /*Frequency tests from Achaz Genetics 2009: outgroup*/

freqtesto_achaz <- function(sample_size,fr,S,w1,w2)  #/* nomes outgroup */

# fr = sfreq[pop]
# S = Segregating Sites

{
    Test <- 0 
    fr <- fr[2:(length(fr)-1)] # ignorier monomorphe Frequenz Stellen 
    
    if(sample_size < 2 || S == 0){ return(NaN)}

    Th1   <- 0
	  sumw1 <- 0

    for(i in 1:(sample_size-1)) {
		 Th1   <- Th1 + fr[i]*(i*w1[i])
		 sumw1 <- sumw1 + w1[i]
	  }
     Th1 <- Th1/sumw1
    
    Th2   <- 0
	  sumw2 <- 0

    for(i in 1:(sample_size-1)) {
		 Th2   <- Th2   + fr[i]*(i*w2[i])
		 sumw2 <- sumw2 + w2[i]
	  }
	   
    Th2 <-  Th2/sumw2
  
  alfan  <- 0
	
	for(i in 1:(sample_size-1)){
	
		omi   <- omegai(sample_size,i,w1,w2)  # --------- Function omegai
		alfan <- alfan + i*(omi*omi)
	}
	
  betan <- 0
 	
	for(i in 1:(sample_size-1)) {
		omi   <- omegai(sample_size,i,w1,w2)
		sigmaii_val <- sigmaii(sample_size,i)
		betan <- betan + i*i * omi*omi * sigmaii_val
  if((i+1)<=(sample_size-1)){ # for Schleifen in R !!! oh man !
    for(j in (i+1):(sample_size-1)){
      sigmaij_val <- sigmaij(sample_size,j,i)
     	betan <- betan + 2.0 * i*j * omegai(sample_size,i,w1,w2) * omegai(sample_size,j,w1,w2) * sigmaij_val
    }
  }
  }
  
	Test <- (Th1 - Th2)/(sqrt(alfan* S/an(sample_size) + betan*S*(S-1)/(an(sample_size)*an(sample_size)+a2n(sample_size))))
  
  if(is.na(Test)){return(NaN)}
	if (abs(Test) < 1.0E-15){return(0)}

    return(Test)
}
###################################################################################
# /*Frequency tests from Achaz Genetics 2009: NO outgroup*/
###################################################################################

freqtestn_achaz <- function(sample_size,fr,S,w1,w2) #/* NO outgroup */

{

    if(sample_size < 2 || S == 0) {return(NaN)}

    Th1   <- 0
	  sumw1 <- 0
	
    for(i in 1:floor(sample_size/2)) {
		  Th1   <-  Th1 + (fr[i])*w1[i]/(psii(sample_size,i))
		  sumw1 <-  sumw1 + w1[i]
    }

    Th1 <- Th1/sumw1

	  
    Th2   <- 0
	  sumw2 <- 0
	
    for(i in 1:floor(sample_size/2)) {
		  Th2   <-  Th2 + (fr[i])*w2[i]/(psii(sample_size,i))
		  sumw2 <- sumw2 + w2[i]
	  }

    Th2 <- Th2/sumw2
  
	alfan <- 0
	
	for(i in 1:floor(sample_size/2)) {
		omi   <- omegain(sample_size,i,w1,w2)
		psi   <- psii(sample_size,i)
		alfan <- alfan + (omi*omi)/psi
	}

  
 	betan <- 0
 	
	for(i in 1:floor(sample_size/2)) {
	
		omi <- omegain(sample_size,i,w1,w2)
		psi <- psii(sample_size,i)
		betan <- betan + omi/psi * omi/psi * rhoii(sample_size,i)

    if((i+1)<=floor(sample_size/2)){ # Wegen R Schleifen
     for(j in (i+1):floor(sample_size/2)) {
			omj <- omegain(sample_size,j,w1,w2)
			psj <- psii(sample_size,j)
			betan <- betan + 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i)
		 }
		}
	}

	Test <- (Th1 - Th2)/(sqrt(alfan*S/an(sample_size) + betan*S*(S-1.0)/(an(sample_size)*an(sample_size)+a2n(sample_size))))

  if (is.na(Test)){return(NaN)}
	if (abs(Test) < 1.0E-15){return(0)}


    return(Test)
    
}
########## SUBFUNCTIONS ########################################################
an <- function(n)
{
  an <- as.integer(n)
	ani <- 0
	if(n>1){
   for(i in 1:(n-1)){
		 ani <- ani + 1/i
   }
  }
  return (ani)
}

a2n <- function(n)
{
	a2ni <- 0

  if(n>1){ # wegen der Schleife ! oh man
	 for(i in 1:(n-1)){
		 a2ni <- a2ni + (1.0/(i*i))
   }
  }
  
  return (a2ni)
}

bn <- function (n,i)
{
  n <- as.integer(n)
  i <- as.integer(i)
	bni <- (2.0*n)/((n-i+1)*(n-i)) * (an(n+1) - an(i)) - 2.0/((n-i))
	return (bni)
}

sigmaii <- function(n,i)

{
  n <- as.integer(n)
  i <- as.integer(i)

  if((2*i)<n){sii <- bn(n,i+1)}
		
	if((2*i)==n){sii <- 2.0 * (an(n)-an(i))/(n-i) - 1/(i*i)}
	
	if((2*i)>n){sii <- bn(n,i) - 1.0/(i*i)}
	
  return(sii)
}

sigmaij <- function(n,i,j)
{
  n <- as.integer(n)
  i <- as.integer(i)
  j <- as.integer(j)

	if(i < j) {
		ii <- j
		jj <- i
	}
	else {
		if(i==j) {
			return(sigmaii(n,i))
		}
		else {
			ii <- i
			jj <- j
		}
	}
	if((ii+jj)<n){sij <- (bn(n,ii+1) - bn(n,ii))/2.0}
	
	if((ii+jj)==n){sij <- (an(n)-an(ii))/((n-ii)) + (an(n)-an(jj))/((n-jj)) - (bn(n,ii) + bn(n,jj+1))/2.0 - 1.0/(ii*jj)}
	
	if((ii+jj)>n){sij <- (bn(n,jj) - bn(n,jj+1))/2.0 - 1/(ii*jj)}
	
	return(sij)
	
}


omegai <- function(n,i,w1,w2)

{
  n <- as.integer(n)
  i <- as.integer(i)
  
	sumw1 <- 0
	
	for(x in 1:(n-1)){
     sumw1 <- sumw1 + w1[x]
  }
  sumw2 <- 0
  
   for(x in 1:(n-1)){
      sumw2 <- sumw2 + w2[x]
   } 
     
 omi <- w1[i]/sumw1 - w2[i]/sumw2
   
  return(omi)
  
}

psii <- function(n,i)
{
	if(i==n-i) {krond <- 1}
	else{krond <- 0}

  psi <- n/((1 + krond)*i*(n-i))
	return(psi)
	
}

rhoii <- function (n,i)
{

	if(i==(n-i)){krond <- 1}
	else{krond <- 0}

	rhoi  <- (sigmaii(n,i)+sigmaii(n,n-i)+2.0*sigmaij(n,i,n-i))
	rhoi <- rhoi/ ((1.0 + krond) * (1.0 + krond))

  return (rhoi)
  
}

rhoij <- function(n,i,j)
{

	if(i==(n-i)){krondi <- 1}
	else{ krondi <- 0}
	if(j==(n-j)){ krondj <- 1}
	else{ krondj <- 0}

	rhoj  <- (sigmaij(n,i,j)+sigmaij(n,i,n-j)+sigmaij(n,n-i,j)+sigmaij(n,n-i,n-j))
	rhoj  <- rhoj /((1.0 + krondi) * (1.0 + krondj))
	return(rhoj)
}


omegain <- function(n,i,w1,w2)
 
 {

	sumw1 <- 0
	for(x in 1:floor(n/2)){
      sumw1 <- sumw1 + w1[x]
	    sumw2 <- 0
	  for(x in 1:floor(n/2)){ 
        sumw2 <- sumw2 + w2[x]
   }
  } 
	omi <- w1[i]/sumw1 - w2[i]/sumw2

  return(omi)
}



