#define min(x,y) ( ( (x) < (y) ) ? (x) : (y) ) 
#define max(x,y) ( ( (x) > (y) ) ? (x) : (y) )
#define ALLOC_SIZE 1000
void fragments ( int i, int f, unsigned int m, int * mf, int * ind  );
unsigned int extract_dups ( char * xx, unsigned int mm, unsigned int f, int * mf, int * ind, int * dups );

unsigned int ecsmf ( unsigned char * x, unsigned char * t, unsigned int ** Occ, unsigned int * num_of_occ );
unsigned int ecsmf_simple ( unsigned char * x, unsigned char * t, unsigned int ** Occ, unsigned int * num_of_occ );
unsigned int acsmf ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );
unsigned int acsmf_simple ( unsigned char * x, unsigned char * t, unsigned int k, unsigned int ** Occ, unsigned int * num_of_occ );

