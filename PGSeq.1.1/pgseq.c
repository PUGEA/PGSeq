/* **************************************************************************** */
/*        PGSeq C implementation                                          */
/* **************************************************************************** */
#include<Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "o8para.h"
#include "pgseq.h"
#include "pgseq_constants.h"

#include <malloc.h>

static expparam in_param;            
void donlp2(void);

typedef void (*func_void_void_t)(void);
typedef void (*func_void_type_liste_donlp_x_err_t)(IINTEGER type, IINTEGER liste[],
				DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
typedef void (*func_void_liste_shift_donlp_x_grad_t)(IINTEGER liste[], IINTEGER shift ,
				DDOUBLE donlp2_x[], DDOUBLE **grad);
typedef void (*func_void_donlp_x_fx_t)(DDOUBLE donlp2_x[],DDOUBLE *fx);
typedef void (*func_void_donlp_x_gradf_t)(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
typedef void (*func_void_mode_t)(IINTEGER mode);
typedef void (*func_void_t)();

extern func_void_type_liste_donlp_x_err_t econ;
extern func_void_liste_shift_donlp_x_grad_t econgrad;
extern func_void_donlp_x_fx_t ef;
extern func_void_donlp_x_gradf_t egradf;
extern func_void_mode_t eval_extern;
extern func_void_t freemem;
extern func_void_t initialparams;
extern func_void_void_t setup;
extern func_void_void_t solchk;
extern func_void_void_t user_init;
extern func_void_void_t user_init_size;
extern func_void_t allocatemem;

void econ_poga(IINTEGER type, IINTEGER liste[],DDOUBLE donlp2_x[], DDOUBLE con[], LLOGICAL err[]);
void econgrad_poga(IINTEGER liste[], IINTEGER shift ,DDOUBLE donlp2_x[], DDOUBLE **grad);
void ef_poga(DDOUBLE donlp2_x[],DDOUBLE *fx);
void egradf_poga(DDOUBLE donlp2_x[],DDOUBLE gradf[]);
void eval_extern_poga(IINTEGER mode);
void setup_poga();
void solchk_poga();
void user_init_poga(void);
void user_init_size_poga(void);
void initialparams_poga();
void allocatemem_poga();
void freemem_poga();


//* Copyright (c) 1995-2004 by Radford M. Neal   
double digamma (double x)
{
  double r, f, t;

  r = 0;
  while (x<=5) 
  { r -= 1/x;
    x += 1;
  }

  f = 1/(x*x);
  t = f*(-1/12.0 + f*(1/120.0 + f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0 + f*(691/32760.0 + f*(-1/12.0 + f*3617/8160.0)))))));

  return r + log(x) - 0.5/x + t;
}

	
void initialparams_poga()
{
	econ = econ_poga;
	econgrad = econgrad_poga;
	ef = ef_poga;
	egradf = egradf_poga;
	eval_extern = eval_extern_poga;
	initialparams = initialparams_poga;
	setup = setup_poga;
	solchk = solchk_poga;
	user_init = user_init_poga;
	user_init_size = user_init_size_poga;
	allocatemem = allocatemem_poga;
}
	
void allocatemem_poga()
{
	int i;
	for(i=0;i<in_param.AlphaNo+2;i++)
		in_param.parameters[i]=0.0;
}

void getNormc(PyObject* Normc)
{
	int size;
	int i,x;
	size=PyList_Size(Normc);
	in_param.Rep=size;
	
	for (i=0;i<size;i++)
	{
		in_param.NormConst[i]=0.0;
	}
	for (i=0;i<size;i++)
	{
    	in_param.NormConst[i]=PyFloat_AsDouble(PyList_GetItem(Normc,i));
	}
	return;
}

void getgenedata(char *GeneDataP,char *GeneMapP,char *GeneLenP,int IsoNo1,int ExonNo1)//outputpath,pwd
{
	int i,j,k;
	
	FILE *fData,*fMap,*fLen;   //read file *
	fData=fopen(GeneDataP,"r");
	fMap=fopen(GeneMapP,"r");
	fLen=fopen(GeneLenP,"r");
	
	in_param.IsoNo = IsoNo1;
	in_param.ExonNo = ExonNo1;
	
	in_param.AlphaNo = in_param.Rep * in_param.IsoNo;
	
	if(fData==NULL||fMap==NULL||fLen==NULL)
	{
		printf("can not open!");
	}
		
	for(i=0; i<in_param.ExonNo;i++)
		for(j=0;j<in_param.Rep;j++)
		in_param.GeneData[i][j]=0.0;
	for(i=0;i<in_param.ExonNo;i++)
		for(j=0;j<in_param.IsoNo;j++)
		in_param.GeneMap[i][j]=0.0;
	for(i=0;i<in_param.ExonNo;i++)
		in_param.GeneLen[i]=0.0;	
		
	for(i=0; i<in_param.ExonNo;i++)
		for(j=0;j<in_param.Rep;j++)
		fscanf(fData,"%lf",&in_param.GeneData[i][j]);

	for(i=0;i<in_param.ExonNo;i++)
		for(j=0;j<in_param.IsoNo;j++)
		fscanf(fMap,"%lf",&in_param.GeneMap[i][j]);
					
	for(i=0;i<in_param.ExonNo;i++)    //copy
		for(j=in_param.IsoNo;j<in_param.AlphaNo;j++)
		in_param.GeneMap[i][j] = in_param.GeneMap[i][j%in_param.IsoNo]; 
		             
	for(i=0;i<in_param.ExonNo;i++)
		{
		fscanf(fLen,"%lf",&in_param.GeneLen[i]);   
		in_param.GeneLen[i]=(in_param.GeneLen[i]+0.0)/1000;
		}	
	fclose(fData);
	fclose(fMap);
	fclose(fLen);	    		
	return;
}

void  poga_calparameters() 
{
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"

	int i;			                        
	
	initialparams_poga();
	allocatemem_poga();
	
	donlp2();
	// printf("optie %d \n",(int)optite+11);
	return ;
}

double * calexpression(int IsoNo,int Rep)
{
	int i,j,k;
	double c,d;
	double l1[in_param.AlphaNo],l2[in_param.AlphaNo],SigmaF[in_param.AlphaNo];
	double IsoData[in_param.AlphaNo];
	double MapMean[in_param.AlphaNo],MapVar[in_param.AlphaNo],MeanTG[in_param.AlphaNo],VarTG[in_param.AlphaNo];
	double TransMean[in_param.AlphaNo],TransVar[in_param.AlphaNo];
	double GeneMean[in_param.Rep],GeneVar[in_param.Rep],GeneMeanTG[in_param.Rep],GeneVarTG[in_param.Rep];
	double Const,Par5,Par6;
	double pi=3.141592653;
	double GeneLength;
	
	double Par2,Par3;
	double Q[in_param.Rep];
	double alphaExon[in_param.AlphaNo];
	double Z[in_param.AlphaNo];
	double ParArray1[in_param.ExonNo][in_param.Rep],ParArray2[in_param.ExonNo][in_param.Rep];
	
	double *expression;
	expression=(double *)malloc(sizeof(double) * (2*in_param.AlphaNo+2*in_param.Rep));
	
	poga_calparameters();

	// for(i=0;i<in_param.AlphaNo;i++)
	// {
 //             printf("%f\n", in_param.parameters[i]);
	// }

	
	GeneLength=0.0;
	for(j=0;j<in_param.ExonNo;j++)
	{
		GeneLength += in_param.GeneLen[j];	
	}
	
	for(i=0;i<in_param.AlphaNo;i++)
	{
		l1[i]=0.0;   //Grad_a
		l2[i]=0.0;   //
		IsoData[i]=0.0;
		SigmaF[i]=0.0;
		MapMean[i]=0.0;
		MapVar[i]=0.0;
		MeanTG[i]=0.0;
		VarTG[i]=0.0;
		TransMean[i]=0.0;
		TransVar[i]=0.0;
	}
	c=in_param.parameters[in_param.AlphaNo];
	d=in_param.parameters[in_param.AlphaNo+1];

	/**** l1    l2 ************************************************/
	
	for(j=0;j<in_param.AlphaNo;j++)
		alphaExon[j]=0.0;	
			
	for(i=0;i<in_param.Rep;i++)
	{
		Z[i]=0.0;
		Q[j]=0.0;
	}		
		
	for(i=0;i<in_param.ExonNo;i++)
		for(j=0;j<in_param.Rep;j++)
		{
			ParArray1[i][j]=0.0;
			ParArray2[i][j]=0.0;
		}
				
	for(i=0;i<in_param.ExonNo;i++)
	{
		Par2 = 0.0;	
  		Par3 = 0.0;		
		for(j=0;j<in_param.AlphaNo;j++)
			alphaExon[j] = in_param.parameters[j]*in_param.GeneMap[i][j];
				
		for(j=0;j<in_param.Rep;j++)
		{
			Q[j] = 0.0;
			k = in_param.IsoNo*j;       // k start point
			while(k<in_param.IsoNo*(j+1)) // k end point
			{
				Q[j] += alphaExon[k];
				k++; 	
			}
			Z[j] = in_param.NormConst[j] * in_param.GeneLen[i] * Q[j];
		}	
		for(j=0;j<in_param.Rep;j++)
		{
			Par2 += in_param.GeneData[i][j]; 
			Par3 += Z[j];	
		}
		Par2 += c;	
		Par3 += d;	    
	    for(j=0;j<in_param.Rep;j++)
		{
			ParArray1[i][j] = in_param.NormConst[j] * in_param.GeneLen[i] *(in_param.GeneData[i][j] / Z[j] - Par2 / Par3);//l1
			ParArray2[i][j] = in_param.NormConst[j] * in_param.NormConst[j] * in_param.GeneLen[i] *in_param.GeneLen[i]*( (c/d)/Par3 - (c/d)/Z[j]);//l2
		}		
	}
	for(i=0;i<in_param.Rep;i++)			
		for(j=0;j<in_param.ExonNo;j++)
		   for(k=0;k<in_param.IsoNo;k++)
		   {
		   	 l1[i * in_param.IsoNo + k] += ParArray1[j][i] *in_param.GeneMap[j][i * in_param.IsoNo + k];
		   	 l2[i * in_param.IsoNo + k] += ParArray2[j][i] *in_param.GeneMap[j][i * in_param.IsoNo + k]*in_param.GeneMap[j][i * in_param.IsoNo + k];		   	 
		   }
	
	/****l1 l2 end *************************************************************************************/	
	for(i=0;i<in_param.AlphaNo;i++)   //SigmaF
	{
		SigmaF[i]=-1* 1.0 / l2[i];
		MapVar[i]=SigmaF[i];             
		MapMean[i]=l1[i]*in_param.parameters[i]*SigmaF[i]+in_param.parameters[i];  //in_param.parameters[i]=alpha
	} 	
	for(i=0;i<in_param.AlphaNo;i++) 
	{
		k=i/in_param.IsoNo;  //floor
		for(j=0;j<in_param.ExonNo;j++)
		{
			IsoData[i]+=in_param.GeneData[j][k]*in_param.GeneMap[j][i];	
		}
	}	
	//Truncated Gaussian
	for(i=0;i<in_param.AlphaNo;i++) 
	{
		Const=2/(1+erf(MapMean[i]/sqrt(2*MapVar[i])));
		Par5=sqrt(MapVar[i]/(2*pi))*exp((-MapMean[i]*MapMean[i])/(2*MapVar[i]));
		Par6=1-erf(-MapMean[i]/sqrt(2*MapVar[i]));
		
		MeanTG[i]=Const*(Par5+(MapMean[i]/2)*Par6);
		VarTG[i]=Const*((MapVar[i]+(MapMean[i]-MeanTG[i])*(MapMean[i]-MeanTG[i]))*Par6/2.0+(MapMean[i]-2*MeanTG[i])*Par5);
		
		if(fpclassify(Const)==FP_INFINITE)  //INF
		{
			MeanTG[i] = in_param.parameters[i];
			VarTG[i] = SigmaF[i];
		}
	}
	
	//the mean and variance of transcript
	for(i=0;i<in_param.AlphaNo;i++)
	{
		TransMean[i] = (c/d) *MeanTG[i];
		TransVar[i] = (c/d)*(c/d)*VarTG[i];
	}
	
	//the mean and variance of gene
	for(i=0;i<in_param.Rep;i++)
	{
		GeneMean[i]=0.0;
		GeneMeanTG[i]=0.0;
		GeneVar[i]=0.0;
		GeneVarTG[i]=0.0;
	}
	if(in_param.IsoNo==1)
	{
		for(i=0;i<in_param.Rep;i++)
		{
			GeneMean[i] = TransMean[i];
			GeneVar[i] = GeneLength * in_param.NormConst[i] * TransVar[i];        
		}
	}
	else
	{
		for(i=0;i<in_param.Rep;i++)
		{
		   k=in_param.IsoNo*i;       
			while(k<in_param.IsoNo*(i+1)) 
			{
				GeneMeanTG[i] += MeanTG[k];
				GeneVarTG[i] += VarTG[k];
				k++; 	
			}
			GeneVarTG[i] = GeneLength*in_param.NormConst[i]*GeneVarTG[i];  					
		}
		
		for(i=0;i<in_param.Rep;i++)
		{
			GeneMean[i] = GeneMeanTG[i]*c/d;
			GeneVar[i] = (c/d)*(c/d)*GeneVarTG[i];
		}
	}
	
	for(i=0;i<in_param.AlphaNo;i++)
		expression[i]=TransMean[i];
	for(i=0;i<in_param.AlphaNo;i++)
		expression[in_param.AlphaNo+i]=TransVar[i];
	for(i=0;i<in_param.Rep;i++)
		expression[2*in_param.AlphaNo+i]=GeneMean[i];
	for(i=0;i<in_param.Rep;i++)
		expression[2*in_param.AlphaNo+in_param.Rep+i]=GeneVar[i];
	return expression;	
	
}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup_poga(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk_poga(void) {
    #define  X extern
    #include "o8comm.h"
    #undef   X
    #include "o8cons.h"
    
    int i;  
    for(i=0;i<in_param.AlphaNo+2;i++)
    	in_param.parameters[i]=donlp2_x[i+1];

    return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */
void ef_poga(DDOUBLE donlp2_x[],DDOUBLE *fx) {
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    
    int i,j,k;
    double c,d; 
    double Par1,Par2,Par3;  
    double alpha[in_param.AlphaNo];
    double Q[in_param.Rep];
    double alphaExon[in_param.AlphaNo];
    double Z[in_param.AlphaNo];
	
    *fx=0.0;
     
	for (i=0; i<in_param.AlphaNo; i++)	
	   alpha[i] = donlp2_x[i+1];
	      
    c = donlp2_x[in_param.AlphaNo+1];
	d = donlp2_x[in_param.AlphaNo+2];

	for(i=0;i<in_param.Rep;i++)
	   Z[i]=0.0;
	for(j=0;j<in_param.Rep;j++)
		Q[j]=0.0;
	
	for(j=0;j<in_param.AlphaNo;j++)
		alphaExon[j]=0.0;
			
	*fx = in_param.ExonNo * ( c * log(d) - lgamma(c) );	
		
	for(i=0;i<in_param.ExonNo;i++)
	{	
  		Par1 = 0.0;
  		Par2 = 0.0;
  		Par3 = 0.0;	
		
		for(j=0;j<in_param.AlphaNo;j++)
			alphaExon[j] = alpha[j] * in_param.GeneMap[i][j];
		
		for(j=0;j<in_param.Rep;j++)
		{
			Q[j]=0.0;
			k=in_param.IsoNo*j;       
			while(k<in_param.IsoNo*(j+1)) 
			{
				Q[j] += alphaExon[k];
				k++; 	
			}	
			Z[j]=in_param.NormConst[j] * in_param.GeneLen[i] * Q[j];
		}
				
		for(j=0;j<in_param.Rep;j++)
		{
			Par1 += in_param.GeneData[i][j] * log( Z[j] );
			Par2 += in_param.GeneData[i][j]; 
			Par3 += Z[j];	
		}

		Par2 += c;	
		Par3 += d;
		
		*fx = *fx + lgamma(Par2) - Par2*log(Par3) + Par1;  //similar Function
		
	}
		*fx = -*fx;

		// printf("fx:%f\n", *fx);

    return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf_poga(DDOUBLE donlp2_x[],DDOUBLE gradf[]) { 

    #define  X extern
    #include "o8fuco.h"
    #undef   X
    
    int i,j,k;
    double c,d;
    double Par1,Par2,Par3;
	double alpha[in_param.AlphaNo];
	double Q[in_param.Rep];
	double alphaExon[in_param.AlphaNo];
	double Z[in_param.AlphaNo];
	double ParArray[in_param.ExonNo][in_param.Rep];
    
	for (i=0; i<in_param.AlphaNo+2; i++)
	    gradf[i+1] = 0.0;
	
    for (i=0; i<in_param.AlphaNo; i++)
    {
      alpha[i] = donlp2_x[i+1];  
    }		     
    c = donlp2_x[in_param.AlphaNo+1];
	d = donlp2_x[in_param.AlphaNo+2];	
	
	for(j=0;j<in_param.AlphaNo;j++)
		alphaExon[j]=0.0;
	
	for(i=0;i<in_param.Rep;i++)
		Z[i]=0.0;
	for(j=0;j<in_param.Rep;j++)
		Q[j]=0.0;
		
	for(i=0;i<in_param.ExonNo;i++)
		for(j=0;j<in_param.Rep;j++)
			ParArray[i][j]=0.0;
	for(i=0; i<in_param.AlphaNo; i++)
	      in_param.Grad_a[i]=0.0;
            
	gradf[in_param.AlphaNo+1] = in_param.ExonNo * (log(d) - digamma(c)); //psi
	gradf[in_param.AlphaNo+2] = in_param.ExonNo * (c/d);
		
	for(i=0;i<in_param.ExonNo;i++)
	{
		
  		Par1 = 0.0;
  		Par2 = 0.0;
  		Par3 = 0.0;	
		
		for(j=0;j<in_param.AlphaNo;j++)
		{
			alphaExon[j] = alpha[j]*in_param.GeneMap[i][j];
		}
			
		for(j=0;j<in_param.Rep;j++)
		{
			Q[j] = 0.0;
			k = in_param.IsoNo*j;       // k start point
			while(k<in_param.IsoNo*(j+1)) // k end point
			{
				Q[j] += alphaExon[k];
				k++; 	
			}
	
			Z[j] = in_param.NormConst[j] * in_param.GeneLen[i] * Q[j];
		}
		
		for(j=0;j<in_param.Rep;j++)
		{
			Par1 += in_param.GeneData[i][j] * log( Z[j] );
			Par2 += in_param.GeneData[i][j]; 
			Par3 += Z[j];	
		}

		Par2 += c;	
		Par3 += d;

		gradf[in_param.AlphaNo+1] = gradf[in_param.AlphaNo+1] + digamma(Par2)-log(Par3);
	    gradf[in_param.AlphaNo+2] = gradf[in_param.AlphaNo+2] -(Par2 / Par3);
	    
	    for(j=0;j<in_param.Rep;j++)
		{
		ParArray[i][j] = in_param.NormConst[j] * in_param.GeneLen[i] *(in_param.GeneData[i][j] / Z[j] - Par2 / Par3);
		}		
	}

	for(i=0;i<in_param.Rep;i++)			
		for(j=0;j<in_param.ExonNo;j++)
		   for(k=0;k<in_param.IsoNo;k++)
			     in_param.Grad_a[i * in_param.IsoNo + k] += ParArray[j][i] *in_param.GeneMap[j][i * in_param.IsoNo + k];
			 	 
	for(i=0; i<in_param.AlphaNo; i++)
		gradf[i+1] = in_param.Grad_a[i];
	for (i=0; i<in_param.AlphaNo+2; i++)
	  {
	  	 gradf[i+1] = -gradf[i+1];

	  	 // printf("gradf %d\t%f\n",(i+1),gradf[1+i]);
	  }
    return;
}

/* **************************************************************************** */
/*                              donlp2-intv size initialization                 */
/* **************************************************************************** */
void user_init_size_poga(void) {

    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X

    n = in_param.AlphaNo+2;  
	nstep = 40;
    
    nlin   =  0;
    nonlin =  0;
    if(in_param.IsoNo>40)
    iterma = 5000;
    else
    iterma=4000;
}

/* **************************************************************************** */
/*                              donlp2 standard setup                           */
/* **************************************************************************** */
void user_init_poga(void) {

    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #include "o8cons.h"
    #undef   X    

    int i,j;
    silent = 1;
    big = 1.e18;
 
	/* initialise the parameters */
	for (i=1; i<=in_param.AlphaNo; i++)
	{	
		donlp2_x[i] = 10.0;
		low[i] = ALOW;
		up[i] = big;	
	}
	
	donlp2_x[in_param.AlphaNo+1] = 20.0;
	low[in_param.AlphaNo+1] = CLOW;
	up[in_param.AlphaNo+1] = big;
	donlp2_x[in_param.AlphaNo+2] = 20.0;
	low[in_param.AlphaNo+2] = DLOW;
	up[in_param.AlphaNo+2] = big;
    
    analyt = 1;
    // difftype=3;
    // epsfcn=1.e-16;
    epsdif = 1.e-16;  
   
    nreset = n; 
    del0 = 1.0e0;
    tau0 = 1.0e1;
    tau  = 0.1e0;
     
    return;
}
/* **************************************************************************** */
/*                        no nonlinear constraints                              */
/* **************************************************************************** */
void econ_poga(IINTEGER type, IINTEGER liste[], DDOUBLE donlp2_x[], DDOUBLE con[], 
              LLOGICAL err[]) {
             
    #define  X extern
    #include "o8fuco.h"
    #undef   X

    return;
}

/* **************************************************************************** */
/* **************************************************************************** */
void econgrad_poga(IINTEGER liste[], IINTEGER shift ,  DDOUBLE donlp2_x[],
               DDOUBLE **grad) {
              
    #define  X extern
    #include "o8fuco.h"
    #undef   X
    return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern_poga(IINTEGER mode) {

    #define  X extern
    #include "o8comm.h"
    #include "o8fint.h"
    #undef   X
    #include "o8cons.h"

    return;
}

