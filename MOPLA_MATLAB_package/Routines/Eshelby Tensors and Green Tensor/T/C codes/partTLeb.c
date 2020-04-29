#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include <mex.h>
// serial code for T
// indices start from 0 not 1
// 2016-04-26 mengmeng
// for lebedev

double **matrix(int m, int n)
{ int i;
  double **ptr;
  /* allocate the array as one block and initialize all entries to zero */
  /* m by n matrix*/
  ptr = (double **)calloc(m, sizeof(double *));
  ptr[0] = (double *)calloc(m*n, sizeof(double));
  for(i = 1; i < m; i++)
	  ptr[i] = ptr[i-1] + n;
  return (ptr);
}

void *unitvector(double a, double b, double *xi) // a is theta, b is phi in rad
{   
	xi[0]=cos(a)*sin(b);
	xi[1]=sin(a)*sin(b);
	xi[2]=cos(b);
	
}

void Inverse(double **A, double **B)  // inverse of a 4x4 matrix
{
	double det;
	
	det = A[0][0]*A[1][1]*A[2][2]*A[3][3] + A[0][0]*A[1][2]*A[2][3]*A[3][1] + A[0][0]*A[1][3]*A[2][1]*A[3][2]
	     +A[0][1]*A[1][0]*A[2][3]*A[3][2] + A[0][1]*A[1][2]*A[2][0]*A[3][3] + A[0][1]*A[1][3]*A[2][2]*A[3][0]
		 +A[0][2]*A[1][0]*A[2][1]*A[3][3] + A[0][2]*A[1][1]*A[2][3]*A[3][0] + A[0][2]*A[1][3]*A[2][0]*A[3][1]
		 +A[0][3]*A[1][0]*A[2][2]*A[3][1] + A[0][3]*A[1][1]*A[2][0]*A[3][2] + A[0][3]*A[1][2]*A[2][1]*A[3][0]
		 -A[0][0]*A[1][1]*A[2][3]*A[3][2] - A[0][0]*A[1][2]*A[2][1]*A[3][3] - A[0][0]*A[1][3]*A[2][2]*A[3][1]
		 -A[0][1]*A[1][0]*A[2][2]*A[3][3] - A[0][1]*A[1][2]*A[2][3]*A[3][0] - A[0][1]*A[1][3]*A[2][0]*A[3][2]
		 -A[0][2]*A[1][0]*A[2][3]*A[3][1] - A[0][2]*A[1][1]*A[2][0]*A[3][3] - A[0][2]*A[1][3]*A[2][1]*A[3][0]
		 -A[0][3]*A[1][0]*A[2][1]*A[3][2] - A[0][3]*A[1][1]*A[2][2]*A[3][0] - A[0][3]*A[1][2]*A[2][0]*A[3][1];
		 
    if(det==0.0){
		printf("detA is zero!\n");
		exit(1);
	}	

    B[0][0]=(A[1][1]*A[2][2]*A[3][3] + A[1][2]*A[2][3]*A[3][1] + A[1][3]*A[2][1]*A[3][2] - A[1][1]*A[2][3]*A[3][2]
	         - A[1][2]*A[2][1]*A[3][3] - A[1][3]*A[2][2]*A[3][1])/det;
    B[0][1]=(A[0][1]*A[2][3]*A[3][2] + A[0][2]*A[2][1]*A[3][3] + A[0][3]*A[2][2]*A[3][1] - A[0][1]*A[2][2]*A[3][3]
	         - A[0][2]*A[2][3]*A[3][1] - A[0][3]*A[2][1]*A[3][2])/det;
	B[0][2]=(A[0][1]*A[1][2]*A[3][3] + A[0][2]*A[1][3]*A[3][1] + A[0][3]*A[1][1]*A[3][2] - A[0][1]*A[1][3]*A[3][2]
	         - A[0][2]*A[1][1]*A[3][3] - A[0][3]*A[1][2]*A[3][1])/det;
    B[0][3]=(A[0][1]*A[1][3]*A[2][2] + A[0][2]*A[1][1]*A[2][3] + A[0][3]*A[1][2]*A[2][1] - A[0][1]*A[1][2]*A[2][3]
	         - A[0][2]*A[1][3]*A[2][1] - A[0][3]*A[1][1]*A[2][2])/det;			 

    B[1][0]=(A[1][0]*A[2][3]*A[3][2] + A[1][2]*A[2][0]*A[3][3] + A[1][3]*A[2][2]*A[3][0] - A[1][0]*A[2][2]*A[3][3] 
	         - A[1][2]*A[2][3]*A[3][0] - A[1][3]*A[2][0]*A[3][2])/det;
	B[1][1]=(A[0][0]*A[2][2]*A[3][3] + A[0][2]*A[2][3]*A[3][0] + A[0][3]*A[2][0]*A[3][2] - A[0][0]*A[2][3]*A[3][2] 
	         - A[0][2]*A[2][0]*A[3][3] - A[0][3]*A[2][2]*A[3][0])/det;
	B[1][2]=(A[0][0]*A[1][3]*A[3][2] + A[0][2]*A[1][0]*A[3][3] + A[0][3]*A[1][2]*A[3][0] - A[0][0]*A[1][2]*A[3][3] 
	         - A[0][2]*A[1][3]*A[3][0] - A[0][3]*A[1][0]*A[3][2])/det;
	B[1][3]=(A[0][0]*A[1][2]*A[2][3] + A[0][2]*A[1][3]*A[2][0] + A[0][3]*A[1][0]*A[2][2] - A[0][0]*A[1][3]*A[2][2] 
	         - A[0][2]*A[1][0]*A[2][3] - A[0][3]*A[1][2]*A[2][0])/det;

	B[2][0]=(A[1][0]*A[2][1]*A[3][3] + A[1][1]*A[2][3]*A[3][0] + A[1][3]*A[2][0]*A[3][1] - A[1][0]*A[2][3]*A[3][1] 
	         - A[1][1]*A[2][0]*A[3][3] - A[1][3]*A[2][1]*A[3][0])/det;
    B[2][1]=(A[0][0]*A[2][3]*A[3][1] + A[0][1]*A[2][0]*A[3][3] + A[0][3]*A[2][1]*A[3][0] - A[0][0]*A[2][1]*A[3][3] 
	         - A[0][1]*A[2][3]*A[3][0] - A[0][3]*A[2][0]*A[3][1])/det;
    B[2][2]=(A[0][0]*A[1][1]*A[3][3] + A[0][1]*A[1][3]*A[3][0] + A[0][3]*A[1][0]*A[3][1] - A[0][0]*A[1][3]*A[3][1] 
	         - A[0][1]*A[1][0]*A[3][3] - A[0][3]*A[1][1]*A[3][0])/det;
    B[2][3]=(A[0][0]*A[1][3]*A[2][1] + A[0][1]*A[1][0]*A[2][3] + A[0][3]*A[1][1]*A[2][0] - A[0][0]*A[1][1]*A[2][3] 
	         - A[0][1]*A[1][3]*A[2][0] - A[0][3]*A[1][0]*A[2][1])/det;

   	B[3][0]=(A[1][0]*A[2][2]*A[3][1] + A[1][1]*A[2][0]*A[3][2] + A[1][2]*A[2][1]*A[3][0] - A[1][0]*A[2][1]*A[3][2] 
	         - A[1][1]*A[2][2]*A[3][0] - A[1][2]*A[2][0]*A[3][1])/det;
   	B[3][1]=(A[0][0]*A[2][1]*A[3][2] + A[0][1]*A[2][2]*A[3][0] + A[0][2]*A[2][0]*A[3][1] - A[0][0]*A[2][2]*A[3][1] 
	         - A[0][1]*A[2][0]*A[3][2] - A[0][2]*A[2][1]*A[3][0])/det;
   	B[3][2]=(A[0][0]*A[1][2]*A[3][1] + A[0][1]*A[1][0]*A[3][2] + A[0][2]*A[1][1]*A[3][0] - A[0][0]*A[1][1]*A[3][2] 
	         - A[0][1]*A[1][2]*A[3][0] - A[0][2]*A[1][0]*A[3][1])/det;	
   	B[3][3]=(A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][0]*A[1][2]*A[2][1] 
	         - A[0][1]*A[1][0]*A[2][2] - A[0][2]*A[1][1]*A[2][0])/det;			 	
}

void partTGL(size_t n, double i, double j, double k, double l, double *weight, double *theta, double *phi, double *C, double *x, double *result){
	int jj, p, q, r, s, Srow=3, Scol=3, Lrow=3, Lcol=3, indexc; 
	double **Av, **ivAv, *xi, rho, *ft, temp;
	
    //allocate Av and xi, ft
	Av = matrix(4,4);
	ivAv=matrix(4,4);
	xi = (double*) calloc(3, sizeof(double));
	ft = (double*) calloc(n, sizeof(double));
	
	result[0] = 0.0;
	/* calculations */
	for(jj=0; jj<n; jj++)
	{
		unitvector(theta[jj],phi[jj],xi);
		rho = x[0]*x[0]*xi[0]*xi[0]+x[1]*x[1]*xi[1]*xi[1]+x[2]*x[2]*xi[2]*xi[2];
		rho = pow(rho,1.5);   //rho in eq.(4a)
		// Av(p,r) 4 by 4 matrix
		Av[0][3]=Av[3][0]=xi[0];
		Av[1][3]=Av[3][1]=xi[1];
		Av[2][3]=Av[3][2]=xi[2];
		temp=0.0;
		for(p=0; p<3; p++)  
			for(r=0; r<3; r++){
				for(q=0; q<3; q++)
					for(s=0; s<3; s++){
						indexc = (r*Lcol+s)*Scol*Srow+p*Scol+q;
						temp = temp+C[indexc]*xi[q]*xi[s];
					}
			 Av[p][r] = temp;
			 temp=0.0;
			}
	    // inverse Av
		Inverse(Av, ivAv);
		ft[jj]=xi[(int)(j-1)]*xi[(int)(l-1)]*ivAv[(int)(i-1)][(int)(k-1)]/rho;		
	}
	
	for(jj=0; jj<n; jj++)
		result[0] = result[0]+weight[jj]*ft[jj];
	
	free(Av);
	free(ivAv);
	free(xi);
	free(ft);
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	 double i, j, k, l;
	 size_t n;
	 double *weight, *theta, *phi, *C, *x, *result;
	 
	 /* check for proper number of arguments */
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:partTGL:nrhs","Nine inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
	
	/* get input values */
	i      = mxGetScalar(prhs[0]);
	j      = mxGetScalar(prhs[1]);
	k      = mxGetScalar(prhs[2]);
    l      = mxGetScalar(prhs[3]); 
    weight = mxGetPr(prhs[4]);
    theta  = mxGetPr(prhs[5]);
    phi    = mxGetPr(prhs[6]);
    C      = mxGetPr(prhs[7]);
    x      = mxGetPr(prhs[8]);
	
	/* initialize the output result to zero */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	result = mxGetPr(plhs[0]);
	
	/* get numbers of elements */
	n = mxGetN(prhs[4]);

	/* call the computational routine */
	partTGL(n,i,j,k,l,weight,theta,phi,C,x,result);
	 
}