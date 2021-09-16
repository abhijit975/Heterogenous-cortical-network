#include <math.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long SEED = -1973;

//
//  Paul de Bakker, 6 Aug 2000                             
//
//  FROM NUMERICAL RECIPES IN C
//
//  Long period (>2 E18) random number generator of
//  L'Ecuyer with Bays-Durham suffle and added
//  safeguards. Returns a uniform random deviate between
//  0.0 and 1.0 (exclusive of the endpoint values).
//  Call with idum a negative integer to initialise;
//  thereafter, do not alter idum between successive
//  deviates in a sequence. RNMX should approximate the
//  largest floating value that is less than 1.
//

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
} // ran2


//
// Returns a normailly distributed deviate with zero mean and unit variance,
// using ran2(idum) as the source of uniform random deviates.
//
double gauss_ran2(long *idum, double mean, double stddev)
{
  static int iset = 0;
  static double gset;
  double fac,rsq,v1,v2, deviate;

  if (*idum < 0) iset = 0; // Reinitialize
  if (iset == 0) { 
    // 
    // We don't have an extra deviate handy, so pick two uniform numbers
    // in the square extending from -1 to +1 in each direction, and
    // see if they are in the unit circle
    //
    do {
      v1 = 2.0 * ran2(idum) - 1.0; 
      v2 = 2.0 * ran2(idum) - 1.0;
      rsq = v1 * v1 + v2 * v2; 
    } while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
    fac = sqrt(-2.0 * log(rsq)/rsq);

    //
    // Now make the Box-Muller transformation to get two normal deviates.
    // Return one and save the other for next time.
    //
    gset = v1*fac;
    iset = 1;      // Set flag.
    deviate = v2*fac;
  }
  else {       // We have an extra deviate handy,
    iset = 0;      // so unset the flag,
    deviate = gset; // and return it.
  }

  return deviate * stddev + mean;
} // gauss_ran2


