#include <stdlib.h>
#include <math.h>
#include <Rinternals.h>

/* Compile into shared library: gcc -shared -O2 -fPIC -o kernconv.so kernconv.c */
/*                          or: gcc -G -O2 -o kernconv.so kernconv.c      */

double i_pnorm(double x)
{ 
    const double pi = 3.141592653589793; double fac;
    
    fac = 1.0 / sqrt(2.0*pi); x = fabs(x);
    return  fac*exp(-0.5*x*x); 
}

double i_gammln(double xx) // Calculate log(gamma(xx)) a la NR 
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,
			  -86.50532032941677,
			  24.01409824083091,
			  - 1.231739572450155,
			  0.1208650973866179e-2,
			  -0.5395239384953e-5};
    int j;
    
    y=x=xx; tmp=x+5.5; tmp -= (x+0.5)*log(tmp); ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y; return -tmp+log(2.5066282746310005*ser/x);
}

int gammln(double *x, double *l) // Calculate log(gamma(xx)) a la NR 
{
    //double xx; xx=x[0]; l[0] = i_gammln(xx); return(0);
    l[0] = i_gammln(x[0]); return(0);
}

double i_gamm(double xx) // Calculate exp(log(gamma(xx))) a la NR 
{
    return  exp(i_gammln(xx)); 
}

int gamm(double *x, double *l) // Calculate exp(log(gamma(xx))) a la NR 
{
    l[0] = i_gamm(x[0]); return(0);
}

double i_skernel(int d, double *u, int p, int q)  // internal spherical kernel
{ 
    const double pi = 3.141592653589793; 
    int i; double volume, c, v = 1.0, r = 0.0, pp, qq; 
    pp=p; qq=q;
    
    if(p>0)                       // support [-1,1]
    {
	for (i=0; i<d; i++)
	{ 
	    r += u[i]*u[i]; if (r>1.0) return 0.0;
	}
	r = 1 - pow(sqrt(r),pp);
	// printf ("r: %f\n", r);
	volume = pow(pi,0.5*d)/i_gamm(0.5*d+1.0) ;
	if(p==2 && q==2){ c = 0.125*(2.0+d)*(4.0+d)/volume; }           // biweight 
	if(p==2 && q==1){ c = 0.5*(2.0+d)/volume; }                     // epanechnikov
	if(p==2 && q==3){ c = (2.0+d)*(4.0+d)*(6.0+d)/(48.0*volume); }  // triweight
	if(p==1 && q==1){ c = (1.0+d)/volume; }                         // triangle
	if(q==0){         c = 1.0/volume; }                             // uniform
	v = c * pow(r,qq);
    }
    else                         // gaussian
    {                              
	for (i=0; i<d; i++)
	{ 
	    r = i_pnorm(u[i]); v *= r;  
	}
    }
    return (v);
}

int skernel(double *d, double *u, double *p, double *q, double *k)  // external spherical kernel
{ 
    int dd, pp, qq; dd=*d; pp=*p; qq=*q;
    k[0] = i_skernel(dd, u, pp, qq); return(0);
}

int x_pow(double *x, double *p, double *r)  // external power
{ 
    r[0] = pow(x[0],p[0]); return(0);
}

double i_pkernel(int d, double *u, int p, int q)  // internal product kernel
{ 
  int i; double v = 1.0, r;

  if(p>0)                       // support [-1,1]
  {
      for (i=0; i<d; i++)
      { 
	  r = fabs(u[i]); if (r>1.0) return 0.0;
	  if(p==2 && q==2){ r = 1.0-r*r; r = 0.9375*r*r; }    // biweight
	  if(p==2 && q==1){ r = 1.0-r*r; r = 0.75*r; }        // epanechnikov
	  if(p==2 && q==3){ r = 1.0-r*r; r = 1.09375*r*r*r; } // triweight
	  if(p==1 && q==1){ r = 1.0-r; }                      // triangle
	  if(q==0){         r = 0.5; }                        // uniform
	  v *= r;  /* v = v * r */
      }
  }
  else                         // gaussian
  {                              
      for (i=0; i<d; i++)
      { 
	  r = i_pnorm(u[i]); v *= r;  
      }
  }
  return v;
}

int pkernel(double *d, double *u, double *p, double *q, double *k)  // external product kernel
{ 
    int dd, pp, qq; dd=*d; pp=*p; qq=*q;
    k[0] = i_pkernel(dd, u, pp, qq); return(0);
}

double i_kernel(int d, double *u, int p, int q, int product)  // internal kernel
{ 
    if(product==1){ return i_pkernel(d, u, p, q); }else{ return i_skernel(d, u, p, q); }
}

int convol(double *dim, double *x, double *y, double *xg, double *r)

/*    
   Input : dim 4 x 1    n|m|d|c|p|q|product
           x   n x d    has to be sorted after the 1st column!
           y   n x c
           xg  m x d    

   Output: r   m x c 
                                                  */
{
  int pp, qq, prod; long i, j, k, istart=0, nn, mm, dd, cc;
  double *nom, *kern, w, thresh=1.0;

  nn=*(dim+0); mm=*(dim+1); dd=*(dim+2); cc=*(dim+3); 
  pp=*(dim+4); qq=*(dim+5); prod=*(dim+6);
  //printf("cc: %i\n",cc);

  nom   = (double *) malloc(sizeof(double)*cc);
  kern  = (double *) malloc(sizeof(double)*dd);

  if(pp==0){ thresh=5.0; }

  for (j=0; j<mm; j++)  /* loop over grid */
    { 
      for (k=0; k<cc; k++) 
	nom[k] = 0;

      for (i=istart; i<nn; i++)   /* loop over obs  */
	{ 
	  kern[0] = *(x+i) - *(xg+j);
      
	  if (kern[0]<-thresh) 
	    istart=i+1;
	  else
	    { 
	      if (kern[0]>thresh) break;

	      for (k=1; k<dd; k++)
		kern[k] = *(x+i+k*nn) - *(xg+j+k*mm);
	      
	      w = i_kernel(dd, kern, pp, qq, prod);

	      for (k=0; k<cc; k++)        
		nom[k]=nom[k]+w* *(y+i+k*nn);
	    }
	}
      
      for (k=0; k<cc; k++)  
	*(r+j+k*mm) = nom[k];
    }

  free (nom);
  free (kern);

  // printf ("convol (C)\n");
  return 0;
}      


