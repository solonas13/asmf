asmf
====

Approximate String Matching via Filtering

Copyright (C) 2013 Solon P. Pissis.

The following functions for approximate string matching with k-mismatches, under the Hamming distance model, are currently supported:

	unsigned int fpt ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int fpt_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int fpt_simple_ms ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ, unsigned int block_size );


The following functions for approximate circular string matching with k-mismatches, under the Hamming distance model, are currently supported:

	unsigned int acsmf ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int acsmf_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int acsmf_simple_ms ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ, unsigned int block_size );


All functions take as input parameters the pattern x, the text t, and the integer threshold; and then return the starting positions of the occurrences of the rotations of x in t with k-mismatches as output. Functions ending with `ms' take an additional parameter block_size which specifies the length of the non-overlapping fragments into which the text is split. This is useful in the case when the text is too large to fit in memory.
