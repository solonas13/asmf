#include <math.h>
#include <stdlib.h>
#include "trie.h"						//include header for trie

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
