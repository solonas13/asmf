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
#include "sm.h"
#include "filter.h"
#include "defs.h"

using namespace sdsl;
using namespace std;

int * PrefixArray ( unsigned char * s )
{
	unsigned int         i, g, n, f;
   	int       * P;

   	n = strlen ( ( char * ) s );
   	P = ( int * ) calloc ( n , sizeof ( int ) );
   	P[0] = 0;
   	g = f = 0;  
   	for ( i = 1; i < n; ++i )
    	{
      		if ( i < g && P[i - f] != g - i )
          		P[i] = min(P[i - f], g - i);
      		else 
        	{ 
          		g = max(g, i);
          		f = i;
          		while ( g < n && s[g] == s[g - f] )
            			g++;
          		P[i] = g - f;
		}
    	}                                
   	return ( P );
}

unsigned int fpt_simple_ms ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ, unsigned int block_size )
{ 
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );
	unsigned char * T;

	char ** seqs;

        /* Fragment variables */
        unsigned int f = k + 1;
        int   * mf;
        int   * ind;

        /* Duplicates variables */
        int     * dups;
        unsigned int uniq;
        int     * d_occ;        // d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
        int     * l_occ;

        /* Occurrences */
        int matches;
        int * F;
        int * Pos;

	int ** P;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

        if ( m > block_size )
        {
                fprintf( stderr, " Error: Invalid length of block size.\n" );
                return ( 0 );
        }

        if ( block_size > n - m + 1 )
        {
                block_size = n - m + 1;
        }

        /* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
        for ( int i = 0; i < f; i ++ )
                fragments ( i, f, m, mf, ind  );


	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
        uniq = extract_dups ( ( char * ) x, m, f, mf, ind, dups );

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
                                memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                                seqs[i][mf[i]] = '\0';
                        }
                        else //add nothing since it is already added 
                        {
                                seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                                seqs[i][0] = '\0';

                                if ( l_occ[ dups[i] ] == -1 )           //if it the first duplicated fragment
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
                        memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                        seqs[i][mf[i]] = '\0';
                }
        }

	/* Here we have to split the text into blocks */
        unsigned int b_size = block_size;
        unsigned int nb_blocks = ceil ( ( double ) n / b_size );
        unsigned int mod_size = n % b_size; //this is the occ_size of the last block

        for ( int b = 0; b < nb_blocks; b++ )
	{
                unsigned int txt_size = b_size + m - 1; // we add ( m - 1 ) as this is the last possible occurrence starting at position T[mm + b_size].
                unsigned int occ_size = b_size;
                unsigned int txt_index = b * b_size;

                /* This is the last block */
                if ( b == nb_blocks - 1 && mod_size > 0 )
                {
                        b_size = mod_size;
                        if  ( b_size < m ) //if the block size is less than the pattern length break.
                                break;
                        else
                        {
                                txt_size = mod_size;
                                occ_size = mod_size - m + 1;
                        }
                }

                unsigned char * txt;
                txt = ( unsigned char * ) malloc( ( size_t ) ( txt_size + 1 ) * sizeof( unsigned char ) );
                if( ( txt == NULL) )
                {
                        fprintf( stderr, " Error: Cannot allocate memory.\n" );
                        return ( 0 );
                }
                memmove ( txt, &t[txt_index], txt_size );
                txt[txt_size] = '\0';

                unsigned int N = m + txt_size;

		T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
		if( ( T == NULL ) ) 
		{
			fprintf( stderr, " Error: Cannot allocate memory.\n" );
			return ( 0 );
		}

		memmove ( T, x, m );
		memmove ( &T[m], txt, txt_size );
		T[N] = '\0';

		/*F[i] stores the fragment id occuring as the i-th match of it in t*/
		F  = ( int * ) calloc( ( f * occ_size ) , sizeof( int ) );
		if( ( F == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory.\n" );
			return ( 0 );
		}

		/*Pos[i] stores the position that F[i] occurs in t*/
		Pos  = ( int * ) calloc( ( f * occ_size ) , sizeof( int ) );
		if( ( Pos == NULL) )
		{
			fprintf( stderr, " Error: Cannot allocate memory.\n" );
			return ( 0 );
		}

        	/* Compute the fragment's occurrences using Aho Corasick Automaton */
        	filtering ( ( char * ) txt, txt_size, ( char ** ) seqs, f, &matches, F, Pos );

		P = ( int ** ) calloc( ( k + 1 ) , sizeof( int * ) );
		P[0] = PrefixArray ( T );	
		for( int i = 1; i < k + 1; ++i )
			P[i] = ( int * ) calloc( ( N ) , sizeof( int ) );

		/* For i = 0 .. num of matches - 1*/
		for ( int i = 0; i < matches; i ++ )
		{
			int jj = F[i];          	//this is the id of the fragment
			int ii = Pos[i] + m - ind[jj];	

			do                      //do this for all fragments that are equal to fragment jj
			{
				for ( int q = 1; q < k + 1; q ++ )
				{
					int delta = ii + P[q - 1][ii];
					if ( delta <= N )
					{
						do
							delta++;
						while ( delta <= N && T[delta - ii] == T[delta] );
					}
					P[q][ii] = min ( delta - ii, N - ii );
				}

				jj = d_occ[jj]; //get the next frament that is equal to fragment jj

			} while ( jj != -1 );
		}

                unsigned int max_alloc = ALLOC_SIZE;
		for ( int i = m; i < m + occ_size; i ++ )
		{
			if ( P[k][i] >= m )
			{
				if ( ( * num_of_occ ) >= max_alloc )
				{
					( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
					max_alloc += ALLOC_SIZE;
				} 
				( *Occ )[ ( * num_of_occ ) ] = i - m + txt_index;	//text index
				( * num_of_occ ) = ( * num_of_occ ) + 1;
			}
		}
		
		/* Deallocate memory */
		for( int i = 0; i < k + 1; ++i )
			free ( P[i] );
		free ( P );
		free ( T );

		free ( F );
		free ( Pos );
	}

	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}

        for( int i = 0; i < f; ++i )
               	free ( seqs [i] );
        free ( seqs );	
        free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

	return ( 1 );
}

unsigned int fpt_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ )
{ 
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );
	unsigned int N = m + n;
	unsigned char *T;

	char ** seqs;

        /* Fragment variables */
        unsigned int f = k + 1;
        int   * mf;
        int   * ind;

        /* Duplicates variables */
        int     * dups;
        unsigned int uniq;
        int     * d_occ;        // d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
        int     * l_occ;

        /* Occurrences */
        int matches;
        int * F;
        int * Pos;

	int ** P;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], x, m );
        memmove ( &T[m], t, n );
        T[N] = '\0';

        fprintf( stderr, "Input is read in memory.\n");

	 /*F[i] stores the fragment id occuring as the i-th match of it in t*/
        F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
        if( ( F == NULL) )
        {
                fprintf( stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }

        /*Pos[i] stores the position that F[i] occurs in t*/
        Pos  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
        if( ( Pos == NULL) )
        {
                fprintf( stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }

        /* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
        for ( int i = 0; i < f; i ++ )
                fragments ( i, f, m, mf, ind  );


	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
        uniq = extract_dups ( ( char * ) x, m, f, mf, ind, dups );

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
                                memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                                seqs[i][mf[i]] = '\0';
                        }
                        else //add nothing since it is already added 
                        {
                                seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                                seqs[i][0] = '\0';

                                if ( l_occ[ dups[i] ] == -1 )           //if it the first duplicated fragment
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
                        memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                        seqs[i][mf[i]] = '\0';
                }
        }

        /* Compute the fragment's occurrences using Aho Corasick Automaton */
        filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, Pos );

        for( int i = 0; i < f; ++i )
                free ( seqs [i] );
        free ( seqs );	

	P = ( int ** ) calloc( ( k + 1 ) , sizeof( int * ) );
        P[0] = PrefixArray ( T );	
	for( int i = 1; i < k + 1; ++i )
  		P[i] = ( int * ) calloc( ( N ) , sizeof( int ) );

        fprintf(stderr, "Prefix array is computed.\n");

	/* For i = 0 .. num of matches - 1*/
        for ( int i = 0; i < matches; i ++ )
        {
                int jj = F[i];          	//this is the id of the fragment
		int ii = Pos[i] + m - ind[jj];	

                do                      //do this for all fragments that are equal to fragment jj
                {
			for ( int q = 1; q < k + 1; q ++ )
			{
				int delta = ii + P[q - 1][ii];
				if ( delta <= N )
				{
					do
						delta++;
					while ( delta <= N && T[delta - ii] == T[delta] );
				}
				P[q][ii] = min ( delta - ii, N - ii );
			}

			jj = d_occ[jj]; //get the next frament that is equal to fragment jj

                } while ( jj != -1 );
	}
	fprintf( stderr, "All pattern's occurrences are found.\n" );

	unsigned int max_alloc = ALLOC_SIZE;
	for ( int i = m; i < N - m + 1; i ++ )
	{
		if ( P[k][i] >= m )
		{
			if ( ( * num_of_occ ) >= max_alloc )
			{
				( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
				max_alloc += ALLOC_SIZE;
			} 
			( *Occ )[ ( * num_of_occ ) ] = i - m;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}

	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}
	
  	/* Deallocate memory */
	for( int i = 0; i < k + 1; ++i )
    		free ( P[i] );
	free ( P );
  	free ( T );

  	free ( F );
  	free ( Pos );
        free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );

	return ( 1 );
}

unsigned int fpt ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ )
{ 
	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) t );
	unsigned int N = m + n;
	unsigned char *T;
	char ** seqs;

        /* Fragment variables */
        unsigned int f = k + 1;
        int   * mf;
        int   * ind;

        /* Duplicates variables */
        int     * dups;
        unsigned int uniq;
        int     * d_occ;        // d_occ[i] stores the id of the next fragment that is equal to itself a la linked-list manner
        int     * l_occ;

        /* Occurrences */
        int matches;
        int * F;
        int * Pos;

	int ** P;
  	int * SA;
  	int * invSA;
  	int * LCP;

	if ( m > n ) 
	{
		fprintf( stderr, " Error: Invalid length of pattern.\n" );
    		return ( 0 );
	}

  	T = ( unsigned char * ) malloc( ( size_t ) ( N + 1 ) * sizeof( unsigned char ) );
  	if( ( T == NULL) ) 
	{
    		fprintf( stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

        memmove ( &T[0], x, m );
        memmove ( &T[m], t, n );
        T[N] = '\0';

        fprintf( stderr, "Input is read in memory.\n");

	/*F[i] stores the fragment id occuring as the i-th match of it in t*/
        F  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
        if( ( F == NULL) )
        {
                fprintf( stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }

        /*Pos[i] stores the position that F[i] occurs in t*/
        Pos  = ( int * ) calloc( ( f * n ) , sizeof( int ) );
        if( ( Pos == NULL) )
        {
                fprintf( stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }

        /* Partition the pattern in fragments */
        mf  = ( int * ) calloc( ( f ) , sizeof( int ) );
        ind = ( int * ) calloc( ( f ) , sizeof( int ) );
        for ( int i = 0; i < f; i ++ )
                fragments ( i, f, m, mf, ind  );


	/* Check whether there exist duplicated fragments */
        dups  = ( int * ) calloc( ( f ) , sizeof( int ) );
        uniq = extract_dups ( ( char * ) x, m, f, mf, ind, dups );

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
                                memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                                seqs[i][mf[i]] = '\0';
                        }
                        else //add nothing since it is already added 
                        {
                                seqs[i] = ( char * ) malloc( ( 1 ) * sizeof( char ) );
                                seqs[i][0] = '\0';

                                if ( l_occ[ dups[i] ] == -1 )           //if it the first duplicated fragment
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
                        memmove ( &seqs[i][0], &x[ ind[i] ], mf[i] );
                        seqs[i][mf[i]] = '\0';
                }
        }

        /* Compute the fragment's occurrences using Aho Corasick Automaton */
        filtering ( ( char * ) t, n, ( char ** ) seqs, f, &matches, F, Pos );

        for( int i = 0; i < f; ++i )
                free ( seqs [i] );
        free ( seqs );	

	P = ( int ** ) calloc( ( k + 1 ) , sizeof( int * ) );
        P[0] = PrefixArray ( T );	
	for( int i = 1; i < k + 1; ++i )
  		P[i] = ( int * ) calloc( ( N ) , sizeof( int ) );

        fprintf(stderr, "Prefix array is computed.\n");

  	/* Compute the suffix array */
  	SA = (int *) malloc((size_t) ( N + 1 ) * sizeof(int));
  	if( ( SA == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}
	if( divsufsort( T, SA, (saidx_t) ( N + 1 ) ) != 0 ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		exit( EXIT_FAILURE );
  	}
        
        fprintf(stderr, "SA is computed.\n");

	/*Compute the inverse SA array */
  	invSA = (int *) calloc((size_t) ( N + 1 ) , sizeof(int));
  	if( ( invSA == NULL) ) 
	{
    		fprintf(stderr, " Error: Cannot allocate memory.\n" );
    		return ( 0 );
  	}

//	#pragma omp parallel for
	for ( int i = 0; i < N + 1; i ++ )
	{
		invSA [SA[i]] = i;
	}

        fprintf(stderr, "insSA is computed.\n");

	/* Compute the LCP array */
        LCP = LCParray( T, N + 1, SA, invSA );
 
	int_vector<> v( N + 1 , 0 ); // create a vector of length n and initialize it with 0s

//	#pragma omp parallel for
	for ( int i = 0; i < N + 1; i ++ )
	{
		v[i] = LCP[i];
	}

        fprintf(stderr, "LCP is computed.\n");

	rmq_succinct_sct<> rmq(&v);

	util::clear(v);

        fprintf(stderr, "LCP is preprocessed for RMQs.\n");

	/* For i = 0 .. num of matches - 1*/
        for ( int i = 0; i < matches; i ++ )
        {
                int jj = F[i];          	//this is the id of the fragment
		int ii = Pos[i] + m - ind[jj];	

                do                      //do this for all fragments that are equal to fragment jj
                {
			for ( int q = 1; q < k + 1; q ++ )
			{
	                        unsigned int mina = min ( ii + P[q - 1][ii] + 1, N - 1 );
	                        unsigned int minb = min ( P[q - 1][ii] + 1, N - 1 );
                        	if ( mina == minb ) 
                        	{  
                               		P[q][ii] = min ( minb, m ); 
                                	continue; 
                        	} 
				unsigned int l = min ( invSA[ mina ], invSA[ minb ] );
				unsigned int r = max ( invSA[ mina ], invSA[ minb ] );
				P[q][ii]        = min ( minb + LCP[rmq ( l + 1, r ) ], m );
			}

			jj = d_occ[jj]; //get the next frament that is equal to fragment jj

                } while ( jj != -1 );
	}
	fprintf( stderr, "All pattern's occurrences are found.\n" );

	unsigned int max_alloc = ALLOC_SIZE;
	for ( int i = m; i < N - m + 1; i ++ )
	{
		if ( P[k][i] == m )
		{
			if ( ( * num_of_occ ) >= max_alloc )
			{
				( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( max_alloc + ALLOC_SIZE ) * sizeof ( unsigned int ) );
				max_alloc += ALLOC_SIZE;
			} 
			( *Occ )[ ( * num_of_occ ) ] = i - m;	//text index
			( * num_of_occ ) = ( * num_of_occ ) + 1;
		}
	}

	if ( ( * num_of_occ ) > 0  )
	{
		( * Occ )   = ( unsigned int * ) realloc ( ( * Occ ), ( ( * num_of_occ ) ) * sizeof ( unsigned int ) );
		fprintf( stderr, "The occurrences vector is resized to %d.\n", ( * num_of_occ ) );
	}
	
  	/* Deallocate memory */
	for( int i = 0; i < k + 1; ++i )
    		free ( P[i] );
	free ( P );

  	free ( T );
  	free ( SA );
  	free ( invSA );
  	free ( LCP );

  	free ( F );
  	free ( Pos );
        free ( mf );
        free ( ind );
        free ( dups );
        free ( d_occ );
        free ( l_occ );


	return ( 1 );
}
