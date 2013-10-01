asmf
====

Approximate String Matching via Filtering

Copyright (C) 2013 Solon P. Pissis.

The following functions for approximate string matching with k-mismatches, under the Hamming distance model, are currently supported:

	unsigned int fpt ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int fpt_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );


The following functions for approximate circular string matching with k-mismatches, under the Hamming distance model, are currently supported:

	unsigned int acsmf ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
	unsigned int acsmf_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );


All functions take as input arguments the pattern x, the text t, and the integer threshold k<m; and then return the starting positions of the occurrences of the rotations of x in t with k-mismatches as output. 
