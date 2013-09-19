#include <stdlib.h>
#include <string.h>
#include "lcp.h"

int * LCParray ( unsigned char *text, int n, int * SA, int * ISA )
{										
	int i=0, j=0;

	int * LCP = ( int * ) calloc  ( n, sizeof( int ) );

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}
	return ( LCP );
}
