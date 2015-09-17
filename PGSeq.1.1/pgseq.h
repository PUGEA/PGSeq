#include "pgseq_constants.h"

typedef struct 
{
	int ExonNo;
	int Rep;
	int AlphaNo;
	int IsoNo;	
	
	double NormConst[MAX_NUM_Rep];
	
	double GeneData[MAX_NUM_ExonNo][MAX_NUM_Rep]; 
 	double GeneLen[MAX_NUM_ExonNo];
	double GeneMap[MAX_NUM_ExonNo][MAX_NUM_AlphaNo];
	    
    double Grad_a[MAX_NUM_AlphaNo];   
    double parameters[MAX_NUM_AlphaNo]; 
			
} expparam;





