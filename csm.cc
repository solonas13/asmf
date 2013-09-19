#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include "aca.h" 						//include header for Aho Corasick automaton
#include <divsufsort.h> 					//include header for suffix sort
#include <sdsl/rmq_support.hpp>					//include header for range minimum queries
#include "lcp.h"						//include header for LCP array
#include "trie.h"						//include header for trie
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "csm.h"

using namespace sdsl;
using namespace std;

void fragments ( int i, int f, unsigned int m, int * mf, int * ind  )
{
	unsigned int modulo = m % f;
	double nf = (double) ( m ) / f ;
	int first;
	int last;

	if( i < modulo )
	{
		first = i * ( ceil ( nf ) );
		last  = first + ceil ( nf ) - 1;
	}
	else
	{
		first = i * ( floor ( nf ) ) + modulo;
		last  = first + floor ( nf ) - 1;
	}
	ind[i] = first;
	mf[i] = last - first + 1;
}

unsigned int extract_dups ( char * xx, unsigned int mm, unsigned int f, int * mf, int * ind, int * dups )
{
	AlphaMap *      alphabet = NULL;
        Trie *          trie = NULL;

        /* Create an empty alphabet */
        alphabet = alpha_map_new();

        /* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

        /* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	unsigned int uniq = 0;

        for( int i = 0; i < f; ++i )
        {
		TrieData data;
                int Nf   = mf[i];

		AlphaChar * pf;    
		pf = (AlphaChar *) calloc ( ( Nf + 1 ) , sizeof( AlphaChar ) );
		for ( int j = 0; j < Nf; j ++ )
		{
			pf[j] = ( AlphaChar ) xx[ ind[i] + j ];
		}
                pf[Nf] = '\0';

		if ( trie_retrieve ( trie, pf, &data ) != TRUE )
		{
			trie_store ( trie, pf, i );
			dups[ i ] = -1;		//first occurrence
			uniq ++;
		}
		else
		{
			dups[ i ] = data;
		}
 
		free ( pf );
	}
	
	trie_free ( trie );
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	return ( uniq );
}

unsigned int ecsmf ( unsigned char * x, unsigned char * t, unsigned int ** Occ, unsigned int * num_of_occ )
{
	/* Lengths */ 
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );

	unsigned char *xx;
	char ** seqs;
	unsigned int mm = 2 * m - 1;

	unsigned char *T;
	unsigned char *Tr;
	unsigned char *txx;
	unsigned int N = mm + n;

  	int * SA;
  	int * invSA;
  	int * LCP;

  	int * SAr;
  	int * invSAr;
  	int * LCPr;

	/* Fragment variables */
	unsigned int f = 4;
	int   * mf;
        int   * ind;

	/*dups[i] stores the id of the next fragment that is equal to itself*/
	int   	* dups;
	unsigned int uniq;
	int   	* d_occ;
	int   	* l_occ;

	/* Occurrences */
	int matches;
	int * F;
	int * P;
	int * M;
	unsigned int max_alloc = ALLOC_SIZE;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	xx = ( unsigned char * ) malloc( ( size_t ) ( mm + 1 ) * sizeof( unsigned char ) );
  	if( ( xx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &xx[0], x, m );
        memmove ( &xx[m], x, m - 1 );
        xx[mm] = '\0';

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], xx, mm );
        memmove ( &T[mm], t, n );
        T[N] = '\0';

  	txx = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( txx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &txx[0], t, n );
        memmove ( &txx[n], xx, mm );
        txx[N] = '\0';

  	Tr = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( Tr == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	for( int i = 0; i < N; ++i )
  		Tr[i] = txx[ N - 1 - i ];
        Tr[N] = '\0';

	free ( txx );

	F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( F == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	P  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( P == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i ++ )
        	fragments ( i, f, mm, mf, ind  );

	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
	uniq = extract_dups ( ( char * ) xx, mm, f, mf, ind, dups );

        l_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
        d_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f )
	{
               	seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for ( int i = 0; i < f; i++ )	
		{
			/* Add the fragment once */
			if ( dups[i] < 0 )
			{
                		seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 		memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                		seqs[i][mf[i]] = '\0';
			}
			else //add nothing since it is already added 
			{
                		seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                		seqs[i][0] = '\0';

				if ( l_occ[ dups[i] ] == -1 )
					d_occ[ dups[i] ] = i;
				else
					d_occ[ l_occ[ dups[i] ] ] = i;
				l_occ[ dups[i] ] = i;
			}		
		}	
	}
	else //add all the fragments since there exist no duplicates
	{
                seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for( int i = 0; i < f; ++i )
        	{
                	seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 	memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                	seqs[i][mf[i]] = '\0';
        	}
	}

  	free ( xx );

	/* Compute the fragment's occurrences using Aho Corasick Automaton */
	filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, P );

	for( int i = 0; i < f; ++i )
		free ( seqs [i] );
	free ( seqs );

  	/* Compute the suffix array of T */
  	SA = (int *) malloc( (size_t) ( N ) * sizeof(int));
  	if( ( SA == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}
	if( divsufsort( T, SA, (saidx_t) ( N ) ) != 0 ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}
        
	/* Compute the inverse SA array of T */
  	invSA = (int *) calloc((size_t) ( N ) , sizeof(int));
  	if( ( invSA == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		invSA [SA[i]] = i;
	}

	/* Compute the LCP array of T */
        LCP = LCParray( T, N, SA, invSA );

	free ( T );
	free ( SA );
 
	int_vector<> v( N, 0 ); // create a vector of length n and initialize it with 0s

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		v[i] = LCP[i];
	}

	rmq_succinct_sct<> rmq(&v);

	util::clear(v);

	/* Compute the SA array of Tr */

  	SAr = (int *) malloc((size_t) ( N ) * sizeof(int));
  	if( ( SAr == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}
	if( divsufsort( Tr, SAr, (saidx_t) ( N ) ) != 0 ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}
        
	/* Compute the inverse SA array of Tr */
  	invSAr = (int *) calloc((size_t) ( N ) , sizeof(int));
  	if( ( invSAr == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		invSAr [SAr[i]] = i;
	}

	/* Compute the LCP array of Tr */
        LCPr = LCParray( Tr, N, SAr, invSAr );
 
	free ( Tr );	
	free ( SAr );

	int_vector<> vr( N, 0 ); // create a vector of length n and initialize it with 0s

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		vr[i] = LCPr[i];
	}

	rmq_succinct_sct<> rmqr(&vr);

	util::clear(vr);

	/* An array to mark the valid occurrences */
	M  = ( int * ) calloc( ( n ) , sizeof( int ) );
  	if( ( M == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < matches; i ++ )
	{
		unsigned int jj = F[i]; 	//this is the id of the fragment
		
		do 
		{

		unsigned int Er = 0;	
		unsigned int ii = P[i]; 	//this is its starting position in t

		unsigned int mina = min ( ind[jj] + mf[jj], mm - 1 );
		unsigned int minb = min ( mm + ii + mf[jj], N - 1 );

		unsigned int l = min ( invSA[ mina ], invSA[ minb ] );
		unsigned int r = max ( invSA[ mina ], invSA[ minb ] );
		Er = min ( LCP[rmq ( l + 1, r ) ], m - mf[jj] ); // to ensure that the extension concenrs only fragment jj
		Er = min ( Er, mm - ( ind[jj] + mf[jj] ) );      // to ensure that we do not extend beyond xx[mm-1]
		
		unsigned int El = 0;
		unsigned int iir = n - P[i] - mf[jj];

		unsigned int minar = min ( mm - ind[jj], mm - 1 );
		unsigned int minbr = min ( mm + iir + mf[jj], N - 1 );

		unsigned int lr = min ( invSAr[ minar ], invSAr[ minbr ] ); 
		unsigned int rr = max ( invSAr[ minar ], invSAr[ minbr ] );
		El = min ( LCPr[rmqr ( lr + 1, rr ) ], m - mf[jj] ); // to ensure that the extension concenrs only fragment jj
		El = min ( El, ind[jj] ); 			     // to ensure that we do not extend beyond xx[mm-1]
		
		/* If the full extension is greater than or equal to m */
		if ( Er + El + mf[jj] >= m )
		{
			for ( int j = max ( P[i] - El, P[i] + mf[jj] - m); j <= min ( P[i] + mf[jj] - m + Er, P[i]); j ++ )
			{
				if ( M[j] == 0 )
				{ 
					M[j] = 1;
				}
			}
		}

		jj = d_occ[jj];	//get the next frament that is equal to fragment jj

		} while ( jj != -1 );
	}

	for ( int i = 0; i <= n - m + 1; i++ )
	{
		if ( ( * num_of_occ ) >= max_alloc )
		{
			( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
			max_alloc += ALLOC_SIZE;
		}
		if ( M[i] )
		{ 
			( *Occ )[ ( * num_of_occ ) ] = i;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}

	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}
	
  	free ( invSA );
  	free ( LCP );
  	free ( invSAr );
  	free ( LCPr );

	free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

  	free ( F );
  	free ( P );
  	free ( M );

	return ( 1 );
}

unsigned int ecsmf_simple ( unsigned char * x, unsigned char * t, unsigned int ** Occ, unsigned int * num_of_occ )
{
	/* Lengths */ 
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );

	unsigned char *xx;
	char ** seqs;
	unsigned int mm = 2 * m - 1;

	unsigned char *T;
	unsigned char *Tr;
	unsigned char *txx;
	unsigned int N = mm + n;

	/* Fragment variables */
	unsigned int f = 4;
	int   * mf;
        int   * ind;

	/*dups[i] stores the id of the next fragment that is equal to itself*/
	int   	* dups;
	unsigned int uniq;
	int   	* d_occ;
	int   	* l_occ;

	/* Occurrences */
	int matches;
	int * F;
	int * P;
	int * M;
	unsigned int max_alloc = ALLOC_SIZE;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	xx = ( unsigned char * ) malloc( ( size_t ) ( mm + 1 ) * sizeof( unsigned char ) );
  	if( ( xx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &xx[0], x, m );
        memmove ( &xx[m], x, m - 1 );
        xx[mm] = '\0';

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], xx, mm );
        memmove ( &T[mm], t, n );
        T[N] = '\0';

  	txx = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( txx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &txx[0], t, n );
        memmove ( &txx[n], xx, mm );
        txx[N] = '\0';

  	Tr = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( Tr == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	for( int i = 0; i < N; ++i )
  		Tr[i] = txx[ N - 1 - i ];
        Tr[N] = '\0';

	free ( txx );

	F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( F == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	P  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( P == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i ++ )
        	fragments ( i, f, mm, mf, ind  );

	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
	uniq = extract_dups ( ( char * ) xx, mm, f, mf, ind, dups );

        l_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
        d_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f )
	{
               	seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for ( int i = 0; i < f; i++ )	
		{
			/* Add the fragment once */
			if ( dups[i] < 0 )
			{
                		seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 		memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                		seqs[i][mf[i]] = '\0';
			}
			else //add nothing since it is already added 
			{
                		seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                		seqs[i][0] = '\0';

				if ( l_occ[ dups[i] ] == -1 )
					d_occ[ dups[i] ] = i;
				else
					d_occ[ l_occ[ dups[i] ] ] = i;
				l_occ[ dups[i] ] = i;
			}		
		}	
	}
	else //add all the fragments since there exist no duplicates
	{
                seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for( int i = 0; i < f; ++i )
        	{
                	seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 	memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                	seqs[i][mf[i]] = '\0';
        	}
	}


	/* Compute the fragment's occurrences using Aho Corasick Automaton */
	filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, P );

	for( int i = 0; i < f; ++i )
		free ( seqs [i] );
	free ( seqs );

	/* An array to mark the valid occurrences */
	M  = ( int * ) calloc( ( n ) , sizeof( int ) );
  	if( ( M == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < matches; i ++ )
	{
		unsigned int jj = F[i]; 	//this is the id of the fragment

		do
		{

		unsigned int ii = P[i]; 	//this is its starting position in t
		unsigned int Er = 0;	

		unsigned int mina = min ( ind[jj] + mf[jj], mm - 1 );
		unsigned int minb = min ( mm + ii + mf[jj], N - 1 );

		while ( T[mina++] == T[minb++] )  
			Er ++;
		Er = min ( Er, m - mf[jj] ); 			 // to ensure that the extension concenrs only fragment jj
		Er = min ( Er, mm - ( ind[jj] + mf[jj] ) );      // to ensure that we do not extend beyond xx[mm-1]
		
		unsigned int El = 0;
		unsigned int iir = n - P[i] - mf[jj];

		unsigned int minar = min ( mm - ind[jj], mm - 1 );
		unsigned int minbr = min ( mm + iir + mf[jj], N - 1 );

		while ( Tr[minar++] == Tr[minbr++] )  
			El ++;
		El = min ( El, m - mf[jj] ); 			   // to ensure that the extension concenrs only fragment jj
		El = min ( El, ind[jj] ); 			   // to ensure that we do not extend beyond xx[mm-1]

		/* If the full extension is greater than or equal to m */
		if ( Er + El + mf[jj] >= m )
		{
			for ( int j = max ( P[i] - El, P[i] + mf[jj] - m); j <= min ( P[i] + mf[jj] - m + Er, P[i]); j ++ )
			{
				if ( M[j] == 0 )
				{ 
					M[j] = 1;
				}
			}
		}

		jj = d_occ[jj];	//get the next frament that is equal to fragment jj

		} while ( jj != -1 );
	}

	for ( int i = 0; i <= n - m + 1; i++ )
	{
		if ( ( * num_of_occ ) >= max_alloc )
		{
			( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
			max_alloc += ALLOC_SIZE;
		}
		if ( M[i] )
		{ 
			( *Occ )[ ( * num_of_occ ) ] = i;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}

	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}

	free ( T );
	free ( Tr );	
  	free ( xx );
	
	free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

  	free ( F );
  	free ( P );
  	free ( M );

	return ( 1 );
}

unsigned int acsmf ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ )
{
	/* Lengths */
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );

	unsigned char *xx;
	unsigned int mm = 2 * m - 1;

	unsigned char *T;
	unsigned char *Tr;
	unsigned char *txx;
	unsigned int N = mm + n;

	char ** seqs;
  	int * SA;
  	int * invSA;
  	int * LCP;

  	int * SAr;
  	int * invSAr;
  	int * LCPr;

	/* Fragment variables */
	unsigned int f = 2 * k + 4;
	int   * mf;
        int   * ind;

	/* Duplicates variables */
	int   	* dups;		
	unsigned int uniq;
	int   	* d_occ;	// d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
	int   	* l_occ;

	/* Occurrences */
	int matches;
	int * F;
	int * P;
	int * M;
	int * mis;
	unsigned int max_alloc = m;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	xx = ( unsigned char * ) malloc( ( size_t ) ( mm + 1 ) * sizeof( unsigned char ) );
  	if( ( xx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &xx[0], x, m );
        memmove ( &xx[m], x, m - 1 );
        xx[mm] = '\0';

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], xx, mm );
        memmove ( &T[mm], t, n );
        T[N] = '\0';

  	txx = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( txx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &txx[0], t, n );
        memmove ( &txx[n], xx, mm );
        txx[N] = '\0';

  	Tr = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( Tr == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	for( int i = 0; i < N; ++i )
  		Tr[i] = txx[ N - 1 - i ];
        Tr[N] = '\0';

	free ( txx );

	/*F[i] stores the fragment id occuring as the i-th match of it in t*/
	F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( F == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/*P[i] stores the position that F[i] occurs in t*/
	P  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( P == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i ++ )
        	fragments ( i, f, mm, mf, ind  );

	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
	uniq = extract_dups ( ( char * ) xx, mm, f, mf, ind, dups );

        l_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
        d_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f )
	{
               	seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for ( int i = 0; i < f; i++ )	
		{
			/* Add the fragment once */
			if ( dups[i] < 0 )
			{
                		seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 		memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                		seqs[i][mf[i]] = '\0';
			}
			else //add nothing since it is already added 
			{
                		seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                		seqs[i][0] = '\0';

				if ( l_occ[ dups[i] ] == -1 )		//if it the first duplicated fragment
					d_occ[ dups[i] ] = i;
				else
					d_occ[ l_occ[ dups[i] ] ] = i;
				l_occ[ dups[i] ] = i;
			}		
		}	
	}
	else //add all the fragments since there exist no duplicated fragment
	{
                seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for( int i = 0; i < f; ++i )
        	{
                	seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 	memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                	seqs[i][mf[i]] = '\0';
        	}
	}

  	free ( xx );

	/* Compute the fragment's occurrences using Aho Corasick Automaton */
	filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, P );

	for( int i = 0; i < f; ++i )
		free ( seqs [i] );
	free ( seqs );

  	/* Compute the suffix array of T */
  	SA = (int *) malloc( ( size_t ) ( N ) * sizeof( int ) );
  	if( SA == NULL ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}
	if( divsufsort( T, SA, ( saidx_t ) ( N ) ) != 0 ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}
        
	/* Compute the inverse SA array of T */
  	invSA = ( int * ) calloc( ( size_t ) ( N ) , sizeof( int ) );
  	if( ( invSA == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		invSA [SA[i]] = i;
	}

	/* Compute the LCP array of T */
        LCP = LCParray( T, N, SA, invSA );

	free ( T );
	free ( SA );
 
	int_vector<> v( N , 0 ); // create a vector of length n and initialize it with 0s

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		v[i] = LCP[i];
	}

	rmq_succinct_sct<> rmq(&v);

	util::clear(v);

	/* Compute the SA array of Tr */
  	SAr = (int *) malloc((size_t) ( N ) * sizeof(int));
  	if( ( SAr == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}
	if( divsufsort( Tr, SAr, (saidx_t) ( N ) ) != 0 ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}
        
	/* Compute the inverse SA array of Tr */
  	invSAr = (int *) calloc((size_t) ( N ) , sizeof(int));
  	if( ( invSAr == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		invSAr [SAr[i]] = i;
	}

	/* Compute the LCP array of Tr */
        LCPr = LCParray( Tr, N, SAr, invSAr );
 
	free ( Tr );	
	free ( SAr );

	int_vector<> vr( N , 0 ); // create a vector of length n and initialize it with 0s

//	#pragma omp parallel for
	for ( int i = 0; i < N; i ++ )
	{
		vr[i] = LCPr[i];
	}

	rmq_succinct_sct<> rmqr(&vr);

	util::clear(vr);

	/* An array to mark the valid occurrences */
	M  = ( int * ) calloc( ( n ) , sizeof( int ) );
  	if( ( M == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	/* For i = 0 .. num of matches - 1*/
	for ( int i = 0; i < matches; i ++ )
	{
		int jj = F[i]; 		//this is the id of the fragment

		do			//do this for all fragments that are equal to fragment jj
		{

		unsigned int Er = 0;
		unsigned int El = 0;
		unsigned int ii = P[i]; 		//this is its starting position in t
		unsigned int iir = n - P[i] - mf[jj];	//this is its starting position in rev(t)

	        /* An array to mark the mismatches for every matching fragment */
		mis  = ( int * ) calloc( mm , sizeof( int ) ); 		//This size could actually be ( 2 * m - mf[jj] ) --- doesn't matter
  		if( ( mis == NULL) ) 
		{
    			fprintf( stderr, " Error: Cannot allocate memory.\n" );
    			return ( 0 );
  		}

		/* For j = 0 .. k */
		for ( int j = 0; j <= k; j++ )
		{
			/* Find the j-th right extension */
			unsigned int mina = min ( ind[jj] + mf[jj] + Er, mm - 1 );
			unsigned int minb = min ( mm + ii + mf[jj] + Er, N - 1 );

			unsigned int l = min ( invSA[ mina ], invSA[ minb ] );
			unsigned int r = max ( invSA[ mina ], invSA[ minb ] );
			Er += min ( LCP[rmq ( l + 1, r ) ], m - mf[jj] );
			Er = min ( Er, m - mf[jj] ); 		     // to ensure that the extension concenrs only fragment jj
			Er = min ( Er, mm - ( ind[jj] + mf[jj] ) );  // to ensure that we do not extend beyond xx[mm-1]

			/* Mark the j-th mismatch in array mis */
			if ( ind[jj] + mf[jj] + Er < mm && j != k )	//we check whether it is within the limits of array mis
			{	
				mis[ind[jj] + mf[jj] + Er] = 1;		//here we mark the mismatch with `1'
			}

			/* Increment by one, that is the mismatch added in the extension; we check first IF we are still within the limits */
			if ( ind[jj] + mf[jj] + Er + 1 < mm && j != k )
				Er = Er + 1;				//add the mismatch

			/* Find the j-th left extension */
			unsigned int minar = min ( mm - ind[jj] + El, mm - 1 );
			unsigned int minbr = min ( mm + iir + mf[jj] + El, N - 1 );

			unsigned int lr = min ( invSAr[ minar ], invSAr[ minbr ] );
			unsigned int rr = max ( invSAr[ minar ], invSAr[ minbr ] );
			El += min ( LCPr[rmqr ( lr + 1, rr ) ], m - mf[jj] );
			El = min ( El, m - mf[jj] );              // to ensure that the extension concenrs only fragment jj
			El = min ( El, ind[jj] ); 		  // to ensure that we do not extend beyond xx[mm-1]

			/* Mark the j-th mismatch in array mis */
			if ( ( int ) ( ind[jj] - El - 1 ) >= 0 && j != k )	
			{
				mis[ ind[jj] - El - 1] = 1;		//here we mark the mismatch with `1'
			}

			/* Increment by one, that is the mismatch added in the extension; we check first IF we are still within the limits */
			if ( ind[jj] - El - 1 - 1 >= 0 && j != k )	
				El = El + 1;				//add the mismatch
		}

		/* If the full extension (right, left, and fragment's length) is greater than or equal to m we have a match */
		#if 1
		if ( Er + El + mf[jj] >= m )
		{
			unsigned int extension = Er + El + mf[jj];

			/* Here we just sum up the total num of mismatches for the first position */
			unsigned int total_mis = 0; 
			int start_pos_mis = ( int ) max ( ( int ) ( ind[jj] - El ), 0 );
			int end_pos_mis = ( int ) min ( ind[jj] - El + m - 1, mm - 1 );
			for ( int j = start_pos_mis; j <= end_pos_mis; j ++ )
				total_mis += mis[j];

			/* Here we compute the range of potential valid positions on t */
			int lpos = max ( int ( max ( P[i] - El, P[i] + mf[jj] - m) ), 0 );
			int rpos = min ( int ( min ( P[i] + mf[jj] - m + Er, P[i] ) ), n - m + 1);

			/* Then we go through all potential valid positions by summing up the mismatches for each position */
			for ( int j = 0; j <= extension - m; j ++ )
			{
				int pos = P[i] - El + j;

				/* If it is in the range */
				if ( pos >= lpos && pos <= rpos )
				{
					/* It is not reported and with less than k mismatches */
					if ( M[pos] == 0 && total_mis <= k )
						M[pos] = 1;
				}

				int lmp = start_pos_mis + j;	//chop the leftmost element
				int rmp = end_pos_mis + j + 1;  //add one element from the right

				if ( rmp >= mm ) break;			

				/* Count the num of mismatches for the next position */
				total_mis = total_mis - mis[lmp] + mis[rmp];	
			}
		}
		#endif

		free ( mis );

		jj = d_occ[jj];	//get the next frament that is equal to fragment jj

		} while ( jj != -1 );
	}

	for ( int i = 0; i <= n - m + 1; i++ )
	{
		if ( ( * num_of_occ ) >= max_alloc )
		{
			( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
			max_alloc += ALLOC_SIZE;
		}
		if ( M[i] )
		{ 
			( *Occ )[ ( * num_of_occ ) ] = i;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}
	
	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}
	
  	free ( invSA );
  	free ( LCP );
  	free ( invSAr );
  	free ( LCPr );

	free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

  	free ( F );
  	free ( P );
  	free ( M );

	return ( 1 );
}

unsigned int acsmf_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ )
{
	/* Lengths */
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );

	unsigned char *xx;
	unsigned int mm = 2 * m - 1;

	unsigned char *T;
	unsigned char *Tr;
	unsigned char *txx;
	unsigned int N = mm + n;

	char ** seqs;

	/* Fragment variables */
	unsigned int f = 2 * k + 4;
	int   * mf;
        int   * ind;

	/* Duplicates variables */
	int   	* dups;		
	unsigned int uniq;
	int   	* d_occ;	// d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
	int   	* l_occ;

	/* Occurrences */
	int matches;
	int * F;
	int * P;
	int * M;
	int * mis;
	unsigned int max_alloc = m;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	xx = ( unsigned char * ) malloc( ( size_t ) ( mm + 1 ) * sizeof( unsigned char ) );
  	if( ( xx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &xx[0], x, m );
        memmove ( &xx[m], x, m - 1 );
        xx[mm] = '\0';

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], xx, mm );
        memmove ( &T[mm], t, n );
        T[N] = '\0';

  	txx = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( txx == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &txx[0], t, n );
        memmove ( &txx[n], xx, mm );
        txx[N] = '\0';

  	Tr = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( Tr == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	for( int i = 0; i < N; ++i )
  		Tr[i] = txx[ N - 1 - i ];
        Tr[N] = '\0';

	free ( txx );

	/*F[i] stores the fragment id occuring as the i-th match of it in t*/
	F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( F == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/*P[i] stores the position that F[i] occurs in t*/
	P  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
  	if( ( P == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i ++ )
        	fragments ( i, f, mm, mf, ind  );

	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
	uniq = extract_dups ( ( char * ) xx, mm, f, mf, ind, dups );

        l_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
        d_occ = ( int * ) calloc( ( f ) , sizeof( int ) );
	for ( int i = 0; i < f; i++ )	
	{
		d_occ[i] = -1;
		l_occ[i] = -1;
	}

	/* In case there exist duplicated fragmnents */
	if ( uniq < f )
	{
               	seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for ( int i = 0; i < f; i++ )	
		{
			/* Add the fragment once */
			if ( dups[i] < 0 )
			{
                		seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 		memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                		seqs[i][mf[i]] = '\0';
			}
			else //add nothing since it is already added 
			{
                		seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                		seqs[i][0] = '\0';

				if ( l_occ[ dups[i] ] == -1 )		//if it the first duplicated fragment
					d_occ[ dups[i] ] = i;
				else
					d_occ[ l_occ[ dups[i] ] ] = i;
				l_occ[ dups[i] ] = i;
			}		
		}	
	}
	else //add all the fragments since there exist no duplicated fragment
	{
                seqs = ( char ** ) malloc( ( f ) * sizeof( char * ) );
		for( int i = 0; i < f; ++i )
        	{
                	seqs[i] = ( char * ) malloc( ( mf[i] + 1 ) * sizeof( char ) );
               	 	memmove ( &seqs[i][0], &xx[ ind[i] ], mf[i] );
                	seqs[i][mf[i]] = '\0';
        	}
	}


	/* Compute the fragment's occurrences using Aho Corasick Automaton */
	filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, P );

	for( int i = 0; i < f; ++i )
		free ( seqs [i] );
	free ( seqs );

	/* An array to mark the valid occurrences */
	M  = ( int * ) calloc( ( n ) , sizeof( int ) );
  	if( ( M == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

	/* For i = 0 .. num of matches - 1*/
	for ( int i = 0; i < matches; i ++ )
	{
		int jj = F[i]; 		//this is the id of the fragment

		do			//do this for all fragments that are equal to fragment jj
		{

		unsigned int Er = 0;
		unsigned int El = 0;
		unsigned int ii = P[i]; 		//this is its starting position in t
		unsigned int iir = n - P[i] - mf[jj];	//this is its starting position in rev(t)

	        /* An array to mark the mismatches for every matching fragment */
		mis  = ( int * ) calloc( mm , sizeof( int ) ); 		//This size could actually be ( 2 * m - mf[jj] ) --- doesn't matter
  		if( ( mis == NULL) ) 
		{
    			fprintf( stderr, " Error: Cannot allocate memory.\n" );
    			return ( 0 );
  		}

		/* For j = 0 .. k */
		for ( int j = 0; j <= k; j++ )
		{
			/* Find the j-th right extension */
			unsigned int mina = min ( ind[jj] + mf[jj] + Er, mm - 1 );
			unsigned int minb = min ( mm + ii + mf[jj] + Er, N - 1 );

			while ( T[mina++] == T[minb++] )  
				Er ++;
			Er = min ( Er, m - mf[jj] ); 		     // to ensure that the extension concenrs only fragment jj
			Er = min ( Er, mm - ( ind[jj] + mf[jj] ) );  // to ensure that we do not extend beyond xx[mm-1]

			/* Mark the j-th mismatch in array mis */
			if ( ind[jj] + mf[jj] + Er < mm && j != k )	//we check whether it is within the limits of array mis
			{	
				mis[ind[jj] + mf[jj] + Er] = 1;		//here we mark the mismatch with `1'
			}

			/* Increment by one, that is the mismatch added in the extension; we check first IF we are still within the limits */
			if ( ind[jj] + mf[jj] + Er + 1 < mm && j != k )
				Er = Er + 1;				//add the mismatch

			/* Find the j-th left extension */
			unsigned int minar = min ( mm - ind[jj] + El, mm - 1 );
			unsigned int minbr = min ( mm + iir + mf[jj] + El, N - 1 );

			/* If we could extend within pattern's length: no meaning to extend beyond that */
			while ( Tr[minar++] == Tr[minbr++] )  
				El ++;
			El = min ( El, m - mf[jj] );              // to ensure that the extension concenrs only fragment jj
			El = min ( El, ind[jj] ); 		  // to ensure that we do not extend beyond xx[mm-1]
			
			/* Mark the j-th mismatch in array mis */
			if ( ( int ) ( ind[jj] - El - 1 ) >= 0 && j != k )	
			{
				mis[ ind[jj] - El - 1] = 1;		//here we mark the mismatch with `1'
			}

			/* Increment by one, that is the mismatch added in the extension; we check first IF we are still within the limits */
			if ( ind[jj] - El - 1 - 1 >= 0 && j != k )	
				El = El + 1;				//add the mismatch
		}

		/* If the full extension (right, left, and fragment's length) is greater than or equal to m we have a match */
		#if 1
		if ( Er + El + mf[jj] >= m )
		{
			unsigned int extension = Er + El + mf[jj];

			/* Here we just sum up the total num of mismatches for the first position */
			unsigned int total_mis = 0; 
			int start_pos_mis = ( int ) max ( ( int ) ( ind[jj] - El ), 0 );
			int end_pos_mis = ( int ) min ( ind[jj] - El + m - 1, mm - 1 );
			for ( int j = start_pos_mis; j <= end_pos_mis; j ++ )
				total_mis += mis[j];

			/* Here we compute the range of potential valid positions on t */
			int lpos = max ( int ( max ( P[i] - El, P[i] + mf[jj] - m) ), 0 );
			int rpos = min ( int ( min ( P[i] + mf[jj] - m + Er, P[i] ) ), n - m + 1);

			/* Then we go through all potential valid positions by summing up the mismatches for each position */
			for ( int j = 0; j <= extension - m; j ++ )
			{
				int pos = P[i] - El + j;

				/* If it is in the range */
				if ( pos >= lpos && pos <= rpos )
				{
					/* It is not reported and with less than k mismatches */
					if ( M[pos] == 0 && total_mis <= k )
						M[pos] = 1;
				}

				int lmp = start_pos_mis + j;	//chop the leftmost element
				int rmp = end_pos_mis + j + 1;  //add one element from the right

				if ( rmp >= mm ) break;			

				/* Count the num of mismatches for the next position */
				total_mis = total_mis - mis[lmp] + mis[rmp];	
			}
		}
		#endif

		free ( mis );

		jj = d_occ[jj];	//get the next frament that is equal to fragment jj

		} while ( jj != -1 );
	}

	for ( int i = 0; i <= n - m + 1; i++ )
	{
		if ( ( * num_of_occ ) >= max_alloc )
		{
			( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
			max_alloc += ALLOC_SIZE;
		}
		if ( M[i] )
		{ 
			( *Occ )[ ( * num_of_occ ) ] = i;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}
	
	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}
	
  	free ( xx );
  	free ( T );
  	free ( Tr );

	free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

  	free ( F );
  	free ( P );
  	free ( M );

	return ( 1 );
}
