#include "mp.h"

extern real ScalingFactor;

/* Surface density */
real Sigma(r)
     real r;
{
  real cavity = 1.0;
  real density, Hin, DR, newden, fac;
  if (r < CAVITYRADIUS) cavity = 1.0/CAVITYRATIO; 
  /* This is *not* a steady state */
  /* profile, if a cavity is defined. It first needs */
  /* to relax towards steady state, on a viscous time scale */
  density = cavity*ScalingFactor*SIGMA0*pow(r,-SIGMASLOPE);

  Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);
  if(RMAX <= OUTTRUNC){
  SIGMA0 = ASPECTRATIO*pow(RMAX,-0.5)*pow(RMAX, -1.5);
  SIGMA0 /=QOUT*PI*pow(RMAX, -SIGMASLOPE)*(1.0-pow(RMIN/(RMAX+Hin),0.5));
  } else {
	SIGMA0 = ASPECTRATIO*pow(OUTTRUNC,-0.5)*pow(OUTTRUNC, -1.5);
        SIGMA0 /=QOUT*PI*pow(OUTTRUNC, -SIGMASLOPE)*(1.0-pow(RMIN/(OUTTRUNC+Hin),0.5));
  }
  density = SIGMA0*pow(r,-SIGMASLOPE)*(1.0-pow(RMIN/(r+Hin),0.5));
  if(r > OUTTRUNC){
        DR = 2.0*ASPECTRATIO*OUTTRUNC; /*(RMAX - OUTTRUNC)/pow(log(FOUT), 0.5);*/
	fac = exp(pow((r-OUTTRUNC)/DR, 2.0));
        if(fac > FOUT){
	fac = FOUT;
	}
	density /= fac;
/* /\*Figure out appropriate Gaussian width to truncate*\/ */
/*   DELTAOUT = 0.5*pow(RMAX-OUTTRUNC,2.0)/log(FOUT); */
/*   DELTAIN = 0.5*pow(RMIN-INTRUNC,2.0)/log(FIN); */
/*   DELTAOUT = sqrt(DELTAOUT); */
/*   DELTAIN = sqrt(DELTAIN); */
/*   if(r >= OUTTRUNC){ 
     density /= FOUT;
*/
/*     /\*gauss outer trunct*\/ */
/*     /\* density *= pow(OUTTRUNC/r,3.0/2.0); *\/ */
/*     /\*power law outer trunc s.t. constant Q outside outtrunct. but only for background density of r^-0.5*\/ */
   }
/*   if(r <= INTRUNC){ */
/*     density *= exp(-pow((r-INTRUNC)/DELTAIN,2.0)/2.0); */
/*   } */

  if(density < DENFLOOR*SIGMA0) density = DENFLOOR*SIGMA0;
  return density;
}

void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else
    /*energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);*/
    energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
  return energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Cooling time */
real CoolingTime(r)
     real r;
{
  real eps = 0.01, rp=10.0; 
  real omega, omegap;
  real ct0;
  ct0 = COOLINGTIME0*pow(r,1.5);

//  omega = pow(r, -1.5);
//  omegap= pow(rp, -1.5);
//  ct0 = COOLINGTIME0*2.0*PI/sqrt( pow(omega - omegap, 2.0) + pow(eps*omegap, 2.0));
  return ct0;
}

void FillCoolingTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    CoolingTimeMed[i] = CoolingTime(Rmed[i]);
}

/* Heating source term */
real Qplusinit(r)
     real r;
{
  real qp0, viscosity;
  real Hin, csq, dlogcsq, hole, omega_sq, Rsoft;
  real romega_prime, romegak_prime;

  Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);   
  csq = pow(ASPECTRATIO,2.0)*pow(r,2.0*FLARINGINDEX -1.0); 
  dlogcsq = (2.0*FLARINGINDEX -1.0)/r;
  Rsoft = r + Hin;
  hole = 1.0 - sqrt(RMIN/(Hin + r));
  viscosity = FViscosity(r);

//  romegak_prime = (3.0/2.0)*pow(r, -1.5);

  omega_sq = csq*(2.0*FLARINGINDEX - 1.0 - SIGMASLOPE + 0.5*r*sqrt(RMIN)*pow(Rsoft, -1.5)/hole) + 1.0/r;
  omega_sq/= r*r; 

  romega_prime = dlogcsq*(r*r*omega_sq - 1.0/r) + 0.5*sqrt(RMIN)*csq*( pow(Rsoft,-3.0)/(hole*hole) )*( pow(Rsoft, 1.5) - 1.5*r*sqrt(Rsoft)*hole - 0.5*sqrt(RMIN)*r ) - 1.0/(r*r);
  romega_prime-= 2.0*r*omega_sq;
  romega_prime/= 2.0*r*sqrt(omega_sq);

  qp0 = (1.0/2.0)*viscosity*Sigma(r)*pow(romega_prime, 2.0);
  return qp0;
}

void FillQplus() {
  int i;
  for (i = 0; i < NRAD; i++)
    QplusMed[i] = Qplusinit(Rmed[i]);
}
