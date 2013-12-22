#include <iostream>
#include <cstdlib> 	
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "csm.h"    	// for using the csm call
#include "sm.h"    	// for using the sm call
#include "defs.h"    	// for some definitions
#include "lfs.h"    	// for reading the input ala divsufsort

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

int main(int argc, char **argv)
{
        
  	FILE * fp;
  	const char * fname;
  	unsigned char * x;
  	unsigned char * t;
  	LFS_OFF_T n;
  	LFS_OFF_T m;
  	int err, needclose = 1;
	double 	start = 0, end = 0;

	unsigned int num_of_occ = 0;
	unsigned int * Occ = NULL;

  	/* Check arguments */
  	if(argc != 4) 
	{
    		fprintf(stderr, " Usage: %s <text-file> <pattern-file> <k-mismatches >= 0> \n", argv[0]);
    		return 0;
  	}

  	/* Open a file for reading: the text*/
  	if(strcmp(argv[1], "-") != 0) {
	#if defined(_MSC_VER) && _MSC_VER >= 1400
    	if(fopen_s(&fp, fname = argv[1], "rb") != 0) 
	{
		#else
    		if((fp = fopen(fname = argv[1], "rb")) == NULL) 
		{
	#endif
      			fprintf( stderr, " Error: Cannot open file `%s'.\n", fname);
      			perror(NULL);
      			exit(EXIT_FAILURE);
    		}
  	} 
	else 
	{
		#if defined(_WIN32)
    		_setmode(_fileno(stdin), _O_BINARY);
		#endif
    		fp = stdin;
    		fname = "stdin";
    		needclose = 0;
  	}

  	/* Get the file size */
  	if(LFS_FSEEK(fp, 0, SEEK_END) == 0) 
	{
    		n = LFS_FTELL(fp);
    		rewind(fp);
    		if(n < 0) 
		{
      			fprintf(stderr, " Error: Cannot ftell `%s'.\n", fname);
      			perror(NULL);
      			exit(EXIT_FAILURE);
    		}
    		if(0x7fffffff <= n) 
		{
      			fprintf(stderr, " Error: Input file `%s' is too big.\n", fname);
      			exit(EXIT_FAILURE);
    		}
  	} 
	else 
	{
    		fprintf(stderr, " Error: Cannot fseek `%s'.\n", fname);
    		perror(NULL);
    		exit(EXIT_FAILURE);
  	}

  	/* Allocate 5n bytes of memory */
  	t = ( unsigned char * ) malloc( ( size_t ) n * sizeof( unsigned char ) );
  	if( ( t == NULL ) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit(EXIT_FAILURE);
  	}

  	/* Read n bytes of data */
  	if(fread( t, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) 
	{
    		fprintf(stderr, " Error: %s `%s'.\n",
		argv[0],
      		(ferror(fp) || (! feof(fp))) ? "Cannot read from":"Unexpected EOF in",
      		fname);
    		perror(NULL);
   	 	exit(EXIT_FAILURE);
  	}
  	if(needclose & 1) { fclose(fp); }

  	/* Open a file for reading: the pattern*/
  	if(strcmp(argv[2], "-") != 0) {
	#if defined(_MSC_VER) && _MSC_VER >= 1400
    	if(fopen_s(&fp, fname = argv[2], "rb") != 0) 
	{
		#else
    		if((fp = fopen(fname = argv[2], "rb")) == NULL) 
		{
	#endif
      			fprintf(stderr, " Error: Cannot open file `%s'.\n", fname);
      			perror(NULL);
      			exit(EXIT_FAILURE);
    		}
  	} 
	else 
	{
		#if defined(_WIN32)
    		_setmode(_fileno(stdin), _O_BINARY);
		#endif
    		fp = stdin;
    		fname = "stdin";
    		needclose = 0;
  	}

  	/* Get the file size */
  	if(LFS_FSEEK(fp, 0, SEEK_END) == 0) 
	{
    		m = LFS_FTELL(fp);
    		rewind(fp);
    		if(m < 0) 
		{
      			fprintf(stderr, " Error: Cannot ftell `%s'.\n", argv[0], fname);
      			perror(NULL);
      			exit(EXIT_FAILURE);
    		}
    		if(0x7fffffff <= m) 
		{
      			fprintf(stderr, " Error: Input file `%s' is too big.\n", fname);
      			exit(EXIT_FAILURE);
    		}
  	} 
	else 
	{
    		fprintf( stderr, " Error: Cannot fseek `%s'.\n", fname );
    		perror( NULL );
    		exit( EXIT_FAILURE );
  	}

  	/* Allocate 5n bytes of memory */
  	x = ( unsigned char * ) malloc( ( size_t ) m * sizeof( unsigned char ) );
  	if( ( x == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}

  	/* Read n bytes of data */
  	if( fread( x, sizeof( unsigned char ), ( size_t ) m, fp) != ( size_t ) m ) 
	{
    		fprintf(stderr, " Error: %s `%s'\n.",
      		( ferror( fp ) || ( ! feof( fp ) ) ) ? "Cannot read from":"Unexpected EOF in",
      		fname);
    		perror(NULL);
   	 	exit(EXIT_FAILURE);
  	}
  	if( needclose & 1 ) { fclose( fp ); }
       
   	int k = atoi ( argv[3] );
	x[m - 1] = '\0';
	t[n - 1] = '\0';
	n--;
	m--;

	/* You may change the ALLOC_SIZE definition in file defs.h --- it is resized automatically anyway */
  	Occ = ( unsigned int * ) realloc ( Occ,   ( ALLOC_SIZE ) * sizeof ( unsigned int ) );
  	if( ( Occ == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}

	start = gettime();
        int block_size = 1048576; //1Mb block will be used
	if ( ! ( acsmf_simple_ms ( x, t, k, &Occ, &num_of_occ, block_size ) ) )
	{
		fprintf(stderr, " Error: acsmf_simple_ms() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of acsmf_simple_ms: %lf\n", ( end - start ));

	#if 0
	start = gettime();
	if ( ! ( acsmf_simple ( x, t, k, &Occ, &num_of_occ ) ) )
	{
		fprintf(stderr, " Error: acsmf_simple() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of acsmf_simple: %lf\n", ( end - start ));

	start = gettime();
	if ( ! ( acsmf ( x, t, k, &Occ, &num_of_occ ) ) )
	{
		fprintf(stderr, " Error: acsmf() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of acsmf: %lf\n", ( end - start ));


	start = gettime();
        int block_size = 5; // block will be used
	if ( ! ( fpt_simple_ms ( x, t, k, &Occ, &num_of_occ, block_size ) ) )
	{
		fprintf(stderr, " Error: fpt_simple_ms() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of fpt_simple_ms: %lf\n", ( end - start ));

	start = gettime();
	if ( ! ( fpt_simple ( x, t, k, &Occ, &num_of_occ ) ) )
	{
		fprintf(stderr, " Error: fpt_simple() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of fpt_simple: %lf\n", ( end - start ));

	start = gettime();
	if ( ! ( fpt ( x, t, k, &Occ, &num_of_occ ) ) )
	{
		fprintf(stderr, " Error: fpt() failed.\n" );
		exit(EXIT_FAILURE);
	}
	end = gettime();
	fprintf( stderr, "Elapsed time of fpt: %lf\n", ( end - start ));
	#endif

	fprintf( stderr, "Occ: %d\n", num_of_occ );
	for ( int i = 0; i < num_of_occ; i++ )
		fprintf( stderr, "%d\n", Occ[i] );

  	free ( Occ );
  	free ( t );
  	free ( x );
       
  	return ( 0 );
}
