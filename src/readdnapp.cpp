/*
**
**		ReadDNAplusplus
**
**		replacement for read.dna for use in PopGen
**
**
**
**
**
**
**
*/

//*
//*			INCLUDES
//*

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

//*
//*			DEFINES
//*

	//
	//
#define			PLATFORM_BITS		32

#	define		_FILE_OFFSET_BITS PLATFORM_BITS


#ifdef			DEBUG
#	define			ONDBG			if( true )
#else
#	define			ONDBG			
#endif

//*
//*			STRUCTS
//*



//*
//*			CLASSES
//*


//!
class dynstorage			//
{
public:

	//
	dynstorage(unsigned int inblocksize ){
		blocksize=inblocksize;
		currentblock = 0;
		currentblockpointer=0;
		currentblockindex=0;
		top=0;
		enlarge();
	}
	
	//! ONLY CALL IF CURRENT BLOCK IS EXHAUSTED - OTHERWISE REMAINING BYTES ARE INVALID DATA AND WILL BE CONSIDERED VALID!
	void	enlarge( void )
	{
		chunkheader * newptr;
		if( currentblock == 0 || currentblock->next == 0 )
		{

			newptr = (chunkheader*)malloc( blocksize+sizeof(chunkheader)-1 );
											//-1 because chunkheader contains 1 databyte
			if( newptr == 0 )
				throw "dynstorage : failed malloc to enlarge !\n";
			if( currentblock )
				currentblock->next = newptr;
			newptr->next = 0;
		}
		else
			newptr = currentblock->next;
		
		//switch to next block
		//
		currentblock = newptr;
		currentblockpointer = &currentblock->databegin[0];
		currentblockindex = 0;
		
		//
		//
		if( top == 0 )
			top = newptr;

		//
	}//enlarge
	
	//
	struct chunkheader {
		chunkheader		*next;
		char			databegin[1];
	};

	
	//!
	inline	bool	get( int index, char&c )
	{
		int numchunk = index/blocksize;
		int	numbyte = index%blocksize;
		chunkheader * p = top;
//		printf("get[%d:c%d/b%d]",index,numchunk,numbyte);
		while( p && numchunk-- )
		{
			p = p->next;
		}
		
		//
		if( p == 0 || (p==currentblock && numbyte >= currentblockindex) )
		{
			c=0;
			return false;
		}
		
		//
		c = p->databegin[numbyte];
		
		//
		return true;
		
		//
	}

	//
	inline void		serialise( char * p, int numbytes=1 )
	{
	
		//
		//
		while( numbytes > 0 )
		{

			// determine the maximum number of bytes we can copy
			int numbytescopyable = numbytes>(blocksize-currentblockindex) ? (blocksize-currentblockindex) : numbytes;
			
			//
			memcpy(currentblockpointer,p, numbytescopyable);
			
			//
			currentblockpointer += numbytescopyable;
			currentblockindex += numbytescopyable;
			p += numbytescopyable;
			numbytes -= numbytescopyable;
			
			if( currentblockindex >= blocksize )
				enlarge();
		
		}
	
		
		//
	}//serialisE( buffer )
	
	//
	//
	inline	void	reset(void)
	{
		chunkheader * p = top;
		while( p )
		{
			memset( &p->databegin[0], 0, blocksize );
			p = p->next;
		}
		currentblock=top;
		currentblockpointer=&currentblock->databegin[0];
		currentblockindex=0;
	}
	
	
	//!
	inline	void	serialise( char c )
	{
//		printf("ser[%02d:%d,%d] -> '%s'\n",c,currentblockindex,blocksize,&top->databegin[0]);
		if( currentblockindex >= blocksize )
			enlarge();
		*currentblockpointer = c;
		currentblockpointer++;
		*currentblockpointer=0;
		currentblockindex++;
	}//serialise( char )
	
	//
	chunkheader			*top;
	chunkheader			*currentblock;
	char				*currentblockpointer;
	unsigned int		currentblockindex;
	
	unsigned int		blocksize;
};


//*
//*			DATA
//*

int		returncode = 0;

char		*filename=(char*)"ALN2/10102075.fas";
//long long	filebytelength=0;
uint64_t	filebytelength=0;
char		*filedatabuffer = 0;

/*

		ids <- c("T" ,"t",	"U","u",	"C","c",	"G","g",	"A","a",	"N","n","?",	"-")
		nuks <- c(1,1,		1,1,		2,2,    	3,3,		4,4,    	 5,5,5,			6)
*/

//!	For quick mapping of nukleotide character codes into numeric constants for the biallelic matrix
char	nucleotide_mapping[] = {
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 5	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 16	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,6,5,5,				// 32	: !"#$%&'()*+’-./
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 48	:5123456789:;<=>?
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 64	:@ABCDEFGHIJKLMNO
	5,5,5,5,	1,5,5,5,	5,5,5,5,	5,5,5,5,				// 80	:PQRSTUVWXYZ[\]^_
	5,4,5,2,	5,5,5,3,	5,5,5,5,	5,5,5,5,				// 96	:`abcdefghijklmno
	5,5,5,5,	1,5,5,5,	5,5,5,5,	5,5,5,5,				// 112	:pqrstuvwxyz{|}~

	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 128	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 144	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 160	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 176	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 192	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 208	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5,				// 224	: - mostly nonprintable - 
	5,5,5,5,	5,5,5,5,	5,5,5,5,	5,5,5,5					// 240	: - mostly nonprintable - 
};

//!
dynstorage		rownames(2048);

//*
//*			EXTERNS
//*



//*
//*			CODE
//*






		//************************
		//************************
		//************************
		//************************

	//
	//
	//
#ifndef LINUX
#define      ftello               ftell
#define      fseeko                   fseek
#endif


//!	Returns the filesize of the given file
//long long	filesize( FILE * fh )
//{

uint64_t   filesize( FILE * fh )
{
	uint64_t res = 0;
	uint64_t curpos = ftello(fh);
	
        
        if( fseeko(fh,0,SEEK_END) == 0 )
	      	res = ftello(fh);
	fseeko(fh,curpos,SEEK_SET);
	return res;
}



	//
	//
	//

bool	setMatrixRownames( SEXP mat, int num_rows )
{
	SEXP		dimnamesvec,
				rownamesvec
					;
	
	//	allocate vectors to set rownames
	//
	PROTECT(dimnamesvec = allocVector(VECSXP, 2));
	PROTECT(rownamesvec = allocVector(STRSXP, num_rows));

	//	set rownames as vector elements
	//
	char	rownambuf[512];
	int		rownambufidx=0;
	int		smplnum=0;
	char	cc=0;
	int		i=0;
	while( rownames.get(i,cc) && num_rows > smplnum )
	{
		//
		if( cc )
		{
			if( rownambufidx < sizeof(rownambuf)-1 )
				rownambuf[rownambufidx++] = cc;
		}
		else
		{
			rownambuf[rownambufidx++] = 0;
			SET_STRING_ELT(rownamesvec,smplnum,mkChar(&rownambuf[0]));
			rownambufidx=0;
			smplnum++;
		}
		
		i++;
	}
	
	//	set dimnames
	//
	SET_VECTOR_ELT(dimnamesvec, 0, rownamesvec);
	setAttrib(mat, R_DimNamesSymbol, dimnamesvec);
	
	UNPROTECT(2);
	
	return true;
	
	//
}

	//
	//
	//

SEXP	copyToMatrix( const char *codemat, int numrows, int numcolumns )
{

	//
	SEXP		currentfilename,
				ans = R_NilValue
					;
	int*		rans;
	int			i,j;


	//	allocate matrix
	//
	//
	PROTECT(ans = allocMatrix(INTSXP, numrows, numcolumns));

	//	copy data
	//
	rans = INTEGER(ans);
	for(i = 0; i < numrows; i++)
	{
		for(j = 0; j < numcolumns; j++)
			rans[i + numrows*j] = *codemat++;
	}
	
	//
	//
	setMatrixRownames(ans, numrows);

	//
	UNPROTECT(1);
	
	return ans;
}




	//
	//
	//
	
	

//!
SEXP	processAlignmentFasta( char * fastabuffer, int bufferlength, int &numnucleotides, int& numsamples )
{
	char	*bufferbegin = fastabuffer;					//remember start of buffer
	char	*bufferend = bufferbegin + bufferlength;	//to find out when we reached end of buffer
	char	*translationpointer = bufferbegin;			//pointer to the address where the translated characters are stored
	char	translatedchar;								//result of translation from ACTGN?- to 123456
	int		translationpos = 0;							//
	int		parsestate=0;								//parsing the file has 3 different states - in-between, >headerline and DNA
	int		alignmentlinelength = -1,					//expected length of DNA sequences in the file based on the first one - all lines must have same length!
			currentlinelength = -1;						//length of the currently parsed DNA sequence
	int		numsamplescontained=0;						//number of samples found in the file
	SEXP	res = R_NilValue;
	
	//
	//
	while( fastabuffer < bufferend )
	{
		char c = *fastabuffer++;
		switch( parsestate )
		{
			//
			//	parsestate 0 = neither DNA/RNA nor >headerline
			//
			case 0:
				if( c == '>' )
					parsestate=1;
				break;
			//
			//	parsestate 1 = >headerline , store the names somewhere
			//
			case 1:
				if( c == '\n' )
				{
					currentlinelength=0;
					numsamplescontained++;
					parsestate=2;
					rownames.serialise( (char)0 );
				}
				else if( c != '>' )
					rownames.serialise(c);

				break;
			//
			//	parsestate 2 = DNA/RNA
			//
			case 2:
			
				//
				//
				if( c == '>' )
				{
					parsestate = 1;
					if( alignmentlinelength == -1 )
						alignmentlinelength = currentlinelength;
					else if( alignmentlinelength != currentlinelength )
					{
						Rprintf("ERROR! alignment line length %d mismatch current line length %d!\n",alignmentlinelength, currentlinelength );return(R_NilValue);
					}
				}
				else if( c != '\n' )
				{
					currentlinelength++;
					translatedchar = nucleotide_mapping[c];
					*translationpointer++ = translatedchar;
					translationpos++;
				}
				break;
		}//switch parsestate

	}//while( not at end of fasta data )
	
	//
	int		translationsize = translationpointer - bufferbegin;
	
	numnucleotides = alignmentlinelength;
	numsamples = numsamplescontained;
	
	//
	if( numsamples > 0 && numnucleotides > 0 )
		res = copyToMatrix( filedatabuffer, numsamples, numnucleotides );

	//
	return res;

	//
}



	//
	//
	//
	


//!
SEXP	processAlignmentPhylip( unsigned char * phylipbuffer, int bufferlength, int &numnucleotides, int& numsamples )
{

	//
	//
	unsigned char	*bufferbegin = phylipbuffer;					//remember start of buffer
	unsigned char	*bufferend = bufferbegin + bufferlength;	//to find out when we reached end of buffer
	unsigned char	*translationpointer = bufferbegin;			//pointer to the address where the translated characters are stored
	char			translatedchar;								//result of translation from ACTGN?- to 123456

	int				translationpos = 0,							//
					oldtranslationpos,							//start of translation for the current block of <numsamples> alignment lines
					matrixstride								// = sizeof(int)*numnucleotides
						;
	SEXP			ans = R_NilValue;
	int				*rans=0;
	int				*samplerans;
	int				numpos=0;
	char			lasttrans=0;
	int				lastoffs=0;
	int				nucpos=0;

	
	//	1.	parse first line with numsaples / alignmentlength ints
	//
	//			E.g.	" 34 1743\n"
	//
	if( bufferlength <= 8  )
	{
		//ONDBG printf("	File  length of %d <= 8 OR first character '%c' is not ' '\n",numsamples, bufferlength,*phylipbuffer );
		return false;
								//NOTE: expect 8 to be the minimum size a (albeit useless) phylip file can have :
								//	" 1 1\n"  = 5 bytes (one gene with one nucleotide)
								//	"Q A"	= 3 bytes (a gene named Q with a single nucleotide A)
	}
	
	if( *phylipbuffer == ' ' )
		phylipbuffer++;
									
		//
		//
	numsamples=0;
	while( *phylipbuffer != ' ' )
	{
		if( phylipbuffer >= bufferend )		//fail if too short
		{
			//ONDBG printf("	too short2\n" );
			return ans;
		}
		if( *phylipbuffer >= '0' && *phylipbuffer <= '9' )
			numsamples = (10*numsamples) + ( *phylipbuffer - '0' );
		else
		{
			//ONDBG printf("	got '%c' not 0-9\n",*phylipbuffer );
			return ans;
		}
		phylipbuffer++;
	}
	
	phylipbuffer++;	//skip ' ' between numbers
	
		//
		//
	numnucleotides=0;
	while( *phylipbuffer != '\n' )		//FIXME \n always valid assumption
	{
		if( phylipbuffer >= bufferend )		//fail if too short
		{
			//ONDBG printf("	too short2\n" );
			return ans;
		}
		if( *phylipbuffer >= '0' && *phylipbuffer <= '9' )
			numnucleotides = (10*numnucleotides) + ( *phylipbuffer - '0' );
		else
		{
			//ONDBG printf("	got '%c' not 0-9\n",*phylipbuffer );
			return ans;
		}
		phylipbuffer++;
	}
	
	//
	//printf("	File has %d samples and alignment length of %d\n",numsamples, numnucleotides );
	
	//	skip any kind of newline sequence
	if( *phylipbuffer == '\r' )
		phylipbuffer++;
	if( *phylipbuffer == '\n' )
		phylipbuffer++;
	
	//	allocate Matrix to store result in
	//
	PROTECT(ans = allocMatrix(INTSXP, numsamples, numnucleotides));
	rans = INTEGER(ans);		//get int-pointer to matrix data so that we can fill right now
	
	
	//
	matrixstride = numnucleotides;

	//	Prevent infinite loops
	//
//#ifdef DEBUG
#	define		COUNTGUARD		if( countguard++ > maxcount )		break;
#	define		COUNTRESET		countguard=0;

	int		countguard=0,maxcount=256;
//#else
//#	define		COUNTGUARD
//#	define		COUNTRESET
//#endif
	
	//	#numsamples runs to get a line "<samplename>   <sequence>" and store the sample names.
	//
	//
	for( int  i=0; (i < numsamples) && (phylipbuffer < bufferend) ; i++ )
	{

		//
		//
		while( (*phylipbuffer != ' ') && (*phylipbuffer != '\t') )
		{
			COUNTGUARD;
			rownames.serialise((char*)phylipbuffer,1);		
			phylipbuffer++;
		}
		rownames.serialise((char)0);
		
		while( (*phylipbuffer == ' ') && (*phylipbuffer == ' ') )
			phylipbuffer++;
			;

		
		//
		//
		COUNTRESET;
		
		//
		samplerans = &rans[ 0 + i ];
		
		//
		while( *phylipbuffer != '\n' )
		{
			COUNTGUARD;
			if( *phylipbuffer!=' ')
			{
				*samplerans = (int)(nucleotide_mapping[*phylipbuffer]);
				samplerans+=numsamples;
				numpos++;
				if( i == 0 )
					nucpos++;
			}
			phylipbuffer++;
		}
		
		
		//	skip newline
		phylipbuffer++;
		
		COUNTRESET;

		//
	}//...for all samples

	//
	translationpos = numpos;
	
	//	
	//
	char	        ring[30];
	int		ringpos=0;
	
	//printf("	rans %x\n",rans);
	
	while( phylipbuffer < bufferend )
	{
		//	expect a newline	(empty line)
		//
		if( *phylipbuffer++ != '\n' )
			break;
	
		//
		numpos=0;
	//	printf("	at p %d\n",nucpos);
		for( int  i=0; (i < numsamples) && (phylipbuffer < bufferend) ; i++ )
		{

			while( (*phylipbuffer == ' ') || (*phylipbuffer == '\t') )
			{
				//TODO phylipbuffer ptr check
				if( phylipbuffer >=  bufferend )
					return R_NilValue;
				COUNTGUARD;
				phylipbuffer++;
			}

			//
			//
			COUNTRESET;
			samplerans = &rans[ translationpos + i ];

//			printf("%d[",i);
			for(; *phylipbuffer != '\n' ;  )
			{
			
				//
				if( phylipbuffer >=  bufferend )
					return R_NilValue;
				if( (samplerans - rans) > (numsamples*numnucleotides) )
				{
				//	printf("	at [%d == %x] for spec %d pos %d\n",
// (samplerans - rans),(samplerans - rans),i,nucpos);
					return R_NilValue;
				}
					

				COUNTGUARD;
				if( *phylipbuffer!=' ')
				{
					lasttrans=*phylipbuffer;

					ring[ringpos++]=lasttrans;ringpos%=sizeof(ring);
					
					*samplerans = nucleotide_mapping[*phylipbuffer];
					samplerans+=numsamples;
					numpos++;
					if( i == 0 )
						nucpos++;
				}
				phylipbuffer++;
			}
	//		printf("]\n");
		
			//	skip newline
			phylipbuffer++;
		
			COUNTRESET;

			//
		}//for all samples

		translationpos += numpos;

	}//while not finished

	
	//printf("	translationpos %d, last=%c lastoffs=%d\n",translationpos,lasttrans,lastoffs);
	
	//printf("RING:\n");
	//for( int i=0; i < sizeof(ring);i++)
	//	printf("%c",ring[i%sizeof(ring)]);
	//printf("\n");

	//
	setMatrixRownames( ans , numsamples );
	
	//
	

	//
	//
	UNPROTECT(1);
	return ans;

	//
}



	//
	//
	//
	
SEXP	processAlignmentAny( char * databuffer, int bufferlength, int &numnucleotides, int& numsamples )
{

	//	'>Samplename something bla\n'
	//
	if( *databuffer == '>' )
		return processAlignmentFasta( databuffer, bufferlength, numnucleotides, numsamples );

	//	' 39 1029\n'
	//
	if( *databuffer == ' ' || ( (*databuffer>='0' ) && (*databuffer<='9') ) )
		return processAlignmentPhylip( (unsigned char*)databuffer, bufferlength, numnucleotides, numsamples );


	//	'#Mega\n'
	//
//	if( *databuffer == '#' )
//		return processAlignmentMega( databuffer, bufferlength, numnucleotides, numsamples );
	
	
	//
	return R_NilValue;
}





	//
	//
	//


extern "C"  SEXP readdna(SEXP filenames )
 {
 
 	//
 	//
	int			numfilenames = length(filenames),
				i=0,
				j,
				numsamples=0,
				numbases=0,
				*rx,
				*ry
					;
	int			*rans,
				*cptr
					;
					
	bool		translationsuccess = false;
	
	SEXP		ans = R_NilValue;
	SEXP		currentfilename,
				dimnamesvec,
				rownamesvec
					;
	
	
	//
	currentfilename = STRING_ELT(filenames,i);
	filename = (char*)CHAR(currentfilename);

	//	open file
	//
	FILE * fh = fopen( filename ,"rt");
	if( fh )
	{
	
		//	load complete file
		//
		filebytelength = filesize( fh );

		//
		//
		if( filebytelength <= 0 )
			return(ans);

		//
		//
		filedatabuffer = (char*)malloc( filebytelength );
		if( filedatabuffer==0 )
		{
			Rprintf("(!!) Failed to allocate %lld bytes to load file into memory!\n",filebytelength);
			return(ans);
		}

		//
		//		
		int numbytesloaded = fread( filedatabuffer, 1, (size_t)filebytelength, fh );
		if( numbytesloaded < filebytelength )
		{
			Rprintf("(!!) Only %dL bytes of %dL could be read!\n",numbytesloaded,filebytelength);
			free( filedatabuffer );
			return( ans );
		}
		
		//
		//
		rownames.reset();
	
		//
		//
		ans = processAlignmentAny( filedatabuffer, filebytelength, numbases, numsamples );
		
		if( ans == R_NilValue )
			Rprintf("	Translation failed in file %s!\n",filename);
		else
		{
			//	set the 'path' attribute
			//
			PROTECT(currentfilename = allocVector(STRSXP, 1));
			SET_STRING_ELT(currentfilename,0,mkChar(filename));
			setAttrib(ans, install("path"), currentfilename);
			UNPROTECT(1);
		}

		//
		free( filedatabuffer );
		filedatabuffer = 0;

		//
		fclose( fh );

	}//if ( fh )
	else
		Rprintf("(!!) Could not open file for reading: '%s'\n",filename);


	//
	return(ans);
	
	//
}


		//************************
		//************************
		//************************
		//************************

//!
R_CallMethodDef callMethods[]  = {
	{"readdna", (DL_FUNC) &readdna, 1 },
	{NULL, NULL, 0}
};

//!
void	R_init_readdnapp(DllInfo *info)
{
	/* Register routines, allocate resources. */
	//printf("Loading DLL readdna++ !!\n");
	
	R_registerRoutines(info,0,callMethods,0,0);
}

//!
void	R_unload_readdnapp(DllInfo *info)
{
	/* Release resources. */
	//printf("unloading DLL readdna++ !!\n");
}



