/* pbtk1cpt.c for R deSolve package
   ___________________________________________________

   Model File:  pbtk1cpt.model

   Date:  Sat Jan 12 13:06:10 2019

   Created by:  "mod v5.6.6"
    -- a model preprocessor by Don Maszle
   ___________________________________________________

   Copyright (c) 1993-2017 Free Software Foundation, Inc.

   Model calculations for compartmental model:

   4 States:
     Agutlument = 0.0,
     Acompartment = 0.0,
     Ametabolized = 0.0,
     AUC = 0.0,

   1 Output:
    "Ccompartment",

   0 Inputs:

   3 Parameters:
     vdist = 0,
     ke = 0,
     kgutabs = 1,
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* Model variables: States */
#define ID_Agutlument 0x00000
#define ID_Acompartment 0x00001
#define ID_Ametabolized 0x00002
#define ID_AUC 0x00003

/* Model variables: Outputs */
#define ID_Ccompartment 0x00000

/* Parameters */
static double parms[3];

#define vdist parms[0]
#define ke parms[1]
#define kgutabs parms[2]

/* Forcing (Input) functions */
static double forc[0];


/* Function definitions for delay differential equations */

int Nout=1;
int nr[1]={0};
double ytau[1] = {0.0};

static double yini[4] = {0.0, 0.0, 0.0, 0.0}; /*Array of initial state variables*/

void lagvalue(double T, int *nr, int N, double *ytau) {
  static void(*fun)(double, int*, int, double*) = NULL;
  if (fun == NULL)
    fun = (void(*)(double, int*, int, double*))R_GetCCallable("deSolve", "lagvalue");
  return fun(T, nr, N, ytau);
}

double CalcDelay(int hvar, double dTime, double delay) {
  double T = dTime-delay;
  if (dTime > delay){
    nr[0] = hvar;
    lagvalue( T, nr, Nout, ytau );
}
  else{
    ytau[0] = yini[hvar];
}
  return(ytau[0]);
}

/*----- Initializers */
void initmod (void (* odeparms)(int *, double *))
{
  int N=3;
  odeparms(&N, parms);
}

void initforc (void (* odeforcs)(int *, double *))
{
  int N=0;
  odeforcs(&N, forc);
}


/* Calling R code will ensure that input y has same
   dimension as yini */
void initState (double *y)
{
  int i;

  for (i = 0; i < sizeof(yini) / sizeof(yini[0]); i++)
  {
    yini[i] = y[i];
  }
}

void getParms (double *inParms, double *out, int *nout) {
/*----- Model scaling */

  int i;

  for (i = 0; i < *nout; i++) {
    parms[i] = inParms[i];
  }


  for (i = 0; i < *nout; i++) {
    out[i] = parms[i];
  }
  }
/*----- Dynamics section */

void derivs (int *neq, double *pdTime, double *y, double *ydot, double *yout, int *ip)
{

  yout[ID_Ccompartment] = y[ID_Acompartment] / vdist ;

  ydot[ID_Agutlument] = - kgutabs * y[ID_Agutlument] ;

  ydot[ID_Acompartment] = kgutabs * y[ID_Agutlument] - ke * y[ID_Acompartment] ;

  ydot[ID_Ametabolized] = ke * y[ID_Acompartment] ;

  ydot[ID_AUC] = yout[ID_Ccompartment] ;

} /* derivs */


/*----- Jacobian calculations: */
void jac (int *neq, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip)
{

} /* jac */


/*----- Events calculations: */
void event (int *n, double *t, double *y)
{

} /* event */

/*----- Roots calculations: */
void root (int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip)
{

} /* root */

