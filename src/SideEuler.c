#include "mp.h"

extern boolean OpenInner, NonReflecting, OuterSourceMass, Evanescent;
extern boolean SelfGravity, SGZeroMode, Adiabatic;
extern Pair DiskOnPrimaryAcceleration;
extern int dimfxy;
real Hp0, Hg0, Ht0;

real GasTotalMass (array)
     PolarGrid *array;
{
  int i, j, ns;
  real *density, total = 0.0, fulltotal=0.0;
  ns = array->Nsec;
  density = array->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      total += Surf[i]*density[j+i*ns];
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}

real GasMomentum (Density, Vtheta)
     PolarGrid *Density, *Vtheta;
{
  int i,j,ns;
  real *density, *vtheta, total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vtheta = Vtheta->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 1; j < ns; j++) {
      total += Surf[i]*(density[j+i*ns]+density[j-1+i*ns])*Rmed[i]*(vtheta[j+i*ns]+OmegaFrame*Rmed[i]);
    }
    total += Surf[i]*(density[i*ns]+density[i*ns+ns-1])*Rmed[i]*(vtheta[i*ns]+OmegaFrame*Rmed[i]);
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return 0.5*fulltotal;
}

real GasTotalEnergy (Density, Vrad, Vtheta, Energy)
     PolarGrid *Density, *Vrad, *Vtheta, *Energy;
{
  int i, j, l, ns;
  real *density, *vrad, *vtheta, *energy, *pot;
  real vr_cent, vt_cent;
  real total = 0.0, fulltotal=0.0;
  ns = Density->Nsec;
  density = Density->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  energy = Energy->Field;
  pot = Potential->Field;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &stat);
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      /* centered-in-cell radial velocity */
      vr_cent = (Rmed[i]-Rinf[i])*vrad[l+ns] + (Rsup[i]-Rmed[i])*vrad[l];
      vr_cent /= (Rsup[i]-Rinf[i]);
      /* centered-in-cell azimuthal velocity */
      if (j < ns-1)
	vt_cent = 0.5*(vtheta[l]+vtheta[l+1]) + Rmed[i]*OmegaFrame;
      else
	vt_cent = 0.5*(vtheta[l]+vtheta[i*ns]) + Rmed[i]*OmegaFrame;
      total += 0.5*Surf[i]*density[l]*(vr_cent*vr_cent + vt_cent*vt_cent) + \
	Surf[i]*energy[l] -						\
	Surf[i]*density[l]*pot[l];
      /* Gas total energy is the sum of its kinematic energy, internal energy */
      /* and gravitational potential energy, including self-gravity */
    }
  }
  if (FakeSequential) {
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    fulltotal = total;
  }
  return fulltotal;
}


void CheckMomentumConservation (Density, Vtheta, sys)
     PolarGrid *Density, *Vtheta;
     PlanetarySystem *sys;
{
  FILE *fichmom;
  char name[80];
  int k;
  real totalmomentum, plmom;
  real xplanet, yplanet, vxplanet, vyplanet;
  real rpl, thetapl, vazimpl, masspl;
  real gasmom, planetsmom;
  gasmom = GasMomentum (Density, Vtheta);
  planetsmom = 0.;
  
  for ( k = 0; k < sys->nb; k++ ) {
    xplanet     = sys->x[k];
    yplanet     = sys->y[k];
    rpl         = sqrt( xplanet*xplanet + yplanet*yplanet );
    thetapl     = atan2 (yplanet, xplanet);
    vxplanet    = sys->vx[k];
    vyplanet    = sys->vy[k];
    vazimpl     = -vxplanet*sin(thetapl) + vyplanet*cos(thetapl);
    masspl      = sys->mass[k];
    plmom       = masspl*rpl*vazimpl;
    planetsmom += plmom;
  }
  totalmomentum = gasmom + planetsmom;
  if ( PhysicalTime < 1e-10 ) {
    Hp0 = plmom;
    Hg0 = gasmom;
    Ht0 = totalmomentum;
    printf("time = %lg, Hp0 = %lg, Hg0 = %lg et Ht0 = %lg\n", PhysicalTime, Hp0, Hg0, Ht0);
  }
  if (!CPU_Master) return;
  sprintf (name, "%s%s.dat", OUTPUTDIR, "Momentum");
  fichmom = fopen(name, "a");
  if (fichmom == NULL) {
    fprintf (stderr, "Can't write 'Momentum.dat' file. Aborting.\n");
    prs_exit (1);
  }
  plmom = fabs (plmom - Hp0);
  gasmom = fabs (gasmom - Hg0);
  totalmomentum = fabs (totalmomentum - Ht0);
  fprintf (fichmom, "%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n", PhysicalTime, plmom, gasmom, totalmomentum, totalmomentum / Ht0);
  fclose (fichmom);
}


void DivisePolarGrid (Num, Denom, Res)
     PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      abs[l] = Rmed[i] * cos(2.0*PI*(real)j/(real)ns);
      ord[l] = Rmed[i] * sin(2.0*PI*(real)j/(real)ns);
    }
  }
}
  
Pair ComputeAccel (force, Rho, x, y, rsmoothing, mass)
     Force *force;
     PolarGrid *Rho;
     real x, y, rsmoothing, mass;
{
  Pair acceleration;
  ComputeForce (force, Rho, x, y, rsmoothing, mass, dimfxy);
  if (ExcludeHill) {
    acceleration.x = force->fx_ex_inner+force->fx_ex_outer;
    acceleration.y = force->fy_ex_inner+force->fy_ex_outer;
  } else {
    acceleration.x = force->fx_inner+force->fx_outer;
    acceleration.y = force->fy_inner+force->fy_outer;
  }
  return acceleration;
}

void OpenBoundary (Vrad, Rho, Energy)
     PolarGrid *Vrad, *Rho, *Energy;
{
  int i,j,l,ns,nr;
  real *rho, *vr, *energy;
  real mean_vr, mean_rho;
  /*if (CPU_Rank != 0) return;*/
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  energy = Energy->Field;
  if (CPU_Rank == 0){
    i = 1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns; 
      rho[l-ns] = rho[l] ;		// copy first ring into ghost ring
      energy[l-ns] = energy[l];
      if ((vr[l+ns] > 0.0) || (rho[l] < SigmaMed[0]))
	vr[l] = 0.0; /*we just allow outflow [inwards] */
      else
	vr[l] = vr[l+ns];
    }
  }

  if (CPU_Rank == CPU_Number-1){
    i = nr-1;
#pragma omp parallel for private(l)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rho[l] = rho[l-ns];		// copy first ring into ghost ring
      energy[l] = energy[l-ns];
      if ((vr[l-ns] < 0.0) ) //|| (rho[l] < SigmaMed[nr-2]))
	vr[l] = 0.0; // we just allow outflow [outwards]
      else
	vr[l] = vr[l-ns];
    }
  }
}

void NonReflectingBoundary (Vrad, Vtheta, Rho, Energy) //godon's nrbc
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i, j, l, ns, nr, jp, lp, i_angle;
  real *rho, *vr, *vt, *energy;
  real urad0, utheta0, density0, en0, cs0;
  real dW, dU, dV;
  real en, density, urad, utheta, cs; 
  real Hin;

  energy = Energy->Field;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  vt = Vtheta->Field;
  
  Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);
  
  if(CPU_Rank == 0) OpenBoundary (Vrad, Rho, Energy);
  
  if (CPU_Rank == CPU_Highest) {
    
    //eqm states of the active cells
    i = nr-2;
    density0 = SigmaMed[i];
    en0 = EnergyMed[i];
    cs0 = sqrt(ADIABATICINDEX*(ADIABATICINDEX - 1.0)*en0/density0);
    
    urad0 = -3.0*FViscosity(Rmed[i])/Rmed[i]*(-SIGMASLOPE+0.5);
    utheta0 = G*1.0/Rmed[i]						\
      + pow(ASPECTRATIO,2.0)*pow(Rmed[i],2.0*FLARINGINDEX -1.0)*
      (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 +
       0.5*Rmed[i]*sqrt(RMIN)*pow(Rmed[i]+Hin, -1.5)/(1.0 -
						      sqrt(RMIN/(Rmed[i]+Hin))));
    if(SelfGravity) {
      utheta0 += -Rmed[i]*GLOBAL_AxiSGAccr0[IMIN+i];
    }
    utheta0 = sqrt(utheta0);
    
    //eqm states of the ghosts cells
    i = nr-1;
    density = SigmaMed[i];
    en = EnergyMed[i];
    cs = sqrt(ADIABATICINDEX*(ADIABATICINDEX - 1.0)*en/density);
    
    urad = -3.0*FViscosity(Rmed[i])/Rmed[i]*(-SIGMASLOPE+0.5);
    utheta = G*1.0/Rmed[i]						\
      + pow(ASPECTRATIO,2.0)*pow(Rmed[i],2.0*FLARINGINDEX -1.0)*
      (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 +
       0.5*Rmed[i]*sqrt(RMIN)*pow(Rmed[i]+Hin, -1.5)/(1.0 -
						      sqrt(RMIN/(Rmed[i]+Hin))));
    if(SelfGravity) {
      utheta += -Rmed[i]*GLOBAL_AxiSGAccr0[IMIN+i];
    }
    utheta = sqrt(utheta);

    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lp= j+(i-1)*ns;
  
      dU = (rho[lp] - density0)/density0; 
      dW = (energy[lp] - en0)/en0/ADIABATICINDEX;
      dV = (vr[lp] - urad0)/cs0;
      

      /* vt[l] = utheta;  */
      /* vr[l] = urad; */
      /* rho[l] = density; */
      /* energy[l] = en; */
	
      // we will assume radial flow is always subsonic
      if(vr[lp] < 0.0){
      	vt[l] = utheta;
      	vr[l] = urad + cs*0.5*(dW + dV);
      	rho[l] = density + 0.5*density*(dW + dV);
      	energy[l] = en + 0.5*ADIABATICINDEX*en*(dW + dV);
      }

      if(vr[lp] > 0.0){
      	vt[l] = utheta + cs*(vt[lp] - utheta0)/cs0;
      	vr[l] = urad + cs*0.5*(dW + dV);
      	rho[l] = density + density*( dU + 0.5*(dV - dW) );
      	energy[l] = en + 0.5*ADIABATICINDEX*en*(dW + dV);
      }

    }
  }
}



void EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real vrad0, vtheta0, viscosity, dens0, energ0, Hin;
  real DRMIN, DRMAX, damping, Tin, Tout, lambda;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  /* if(CPU_Rank == CPU_Number -1){ /\*non-reflecting outer bc. remember to set dampradout>1*\/ */
/*     NonReflectingBoundary (Vrad, Vtheta, Rho, Energy); */
/*   } */
/*   OpenBoundary (Vrad, Rho, Energy);/\*use open inner boundary*\/   */
 
  /* Orbital period at inner and outer boundary */
  Tin = 2.0*PI*pow(GlobalRmed[0],3./2)*DAMPRATEIN;
  Tout = 2.0*PI*pow(GlobalRmed[GLOBALNRAD-1],3./2)*DAMPRATEOUT;
  /* DRMIN AND DRMAX are global Radii boundaries of killing wave zones */
  DRMIN = GlobalRmed[0]*DAMPRADIN;
  DRMAX = GlobalRmed[GLOBALNRAD-1]*DAMPRADOUT;
  Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);
  lambda = 0.0;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    if ( (Rmed[i] < DRMIN) || (Rmed[i] > DRMAX) ) {
      /* Damping operates only inside the wave killing zones */
      if (Rmed[i] < DRMIN) {
	damping = (Rmed[i]-DRMIN)/(GlobalRmed[0]-DRMIN);
	lambda = damping*damping*10.0*step/Tin;
      }
      if (Rmed[i] > DRMAX) {
	damping = (Rmed[i]-DRMAX)/(GlobalRmed[GLOBALNRAD-1]-DRMAX);
	lambda = damping*damping*10.0*step/Tout;
      }
      viscosity = FViscosity (Rmed[i]);
      if (!SelfGravity) {
	/*vtheta0 = sqrt ( G*1.0/Rmed[i] *				\
			 ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			   pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) );*/
         vtheta0 = G*1.0/Rmed[i]                                          \
        + pow(ASPECTRATIO,2.0)*pow(Rmed[i],2.0*FLARINGINDEX -1.0)*
        (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[i]*sqrt(RMIN)*pow(Rmed[i]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[i]+Hin))));
      vtheta0 = sqrt(vtheta0);
      }
      if (SelfGravity) {
	/*vtheta0 = sqrt (  G*1.0/Rmed[i] *				\
			  ( 1.0 - (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*	\
			    pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) - \
			  Rmed[i]*GLOBAL_AxiSGAccr[i+IMIN] );*/
      
      vtheta0 = G*1.0/Rmed[i]                                          \
        + pow(ASPECTRATIO,2.0)*pow(Rmed[i],2.0*FLARINGINDEX -1.0)*
        (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[i]*sqrt(RMIN)*pow(Rmed[i]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[i]+Hin)))) - \
        Rmed[i]*GLOBAL_AxiSGAccr0[i+IMIN];
      vtheta0 = sqrt(vtheta0);      
      }
      /* this could be refined if CentrifugalBalance is used... */
      vtheta0 -= Rmed[i]*OmegaFrame;
      vrad0 = -3.0*viscosity/Rmed[i]*(-SIGMASLOPE+.5);
      dens0 = SigmaMed[i];
      energ0 = EnergyMed[i];
      
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
	vrad[l]   = (vrad[l]+lambda*vrad0)/(1.0+lambda);
	vtheta[l] = (vtheta[l]+lambda*vtheta0)/(1.0+lambda);
	dens[l]   = (dens[l]+lambda*dens0)/(1.0+lambda);
	if (Adiabatic)
	  energ[l]  = (energ[l]+lambda*energ0)/(1.0+lambda);
      }
    }
  }
}

void AzimuthalAverage (Vrad, Vtheta, Rho, Energy)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real vrad0, vtheta0, dens0, energ0;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  for (i = Zero_or_active; i < Max_or_active; i++) {

       vrad0   =  0.0;
       vtheta0 =  0.0;
       dens0   =  0.0;
       energ0  =  0.0;

         for(j = 0; j < ns; j++){
            l = i*ns + j;
            vrad0  += vrad[l];
            vtheta0+= vtheta[l];
            dens0  += dens[l];
            if(Adiabatic) energ0+=energ[l];      
         }
             vrad0 /= (real)ns;
             vtheta0/= (real)ns;
             dens0  /= (real)ns;
           if(Adiabatic) energ0/=(real)ns;
                       
    for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vrad[l]   = vrad0;
        vtheta[l] = vtheta0;
        dens[l]   = dens0;
        if (Adiabatic)
          energ[l]  = energ0;
      }
    }
}




void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i, j, l, nr, ns;
  real *rho, average_rho = 0.0, *vr, penul_vr;
  if (CPU_Rank != CPU_Highest) return;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    average_rho += rho[l];
  }
  average_rho /= (real)ns;
  average_rho = SigmaMed[nr-1]-average_rho;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l] += average_rho;
  }
  i = nr-1;
  penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),-SIGMASLOPE);
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    vr[l] = penul_vr;
  }
}

void ApplySubKeplerianBoundary (Vtheta)
     PolarGrid *Vtheta;
{
  int i, j, l, nr, ns;
  real VKepIn, VKepOut, Hin;
  real *vt;
  vt = Vtheta->Field;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
#pragma omp single
  {
    if ( !SelfGravity ) {
      Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);

//      VKepIn = G*1.0/Rmed[0]                                            \
//        + pow(ASPECTRATIO,2.0)*pow(Rmed[0],2.0*FLARINGINDEX -1.0)*
//        (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[0]*sqrt(RMIN)*pow(Rmed[0]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[0]+Hin))));
//      VKepIn = sqrt(VKepIn);

      VKepOut = G*1.0/Rmed[nr-1]                                          \
        + pow(ASPECTRATIO,2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX -1.0)*
        (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[nr-1]*sqrt(RMIN)*pow(Rmed[nr-1]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[nr-1]+Hin))));
      VKepOut = sqrt(VKepOut);
    }

    else {
      if ( !SGZeroMode )
	mpi_make1Dprofile (SG_Accr, GLOBAL_AxiSGAccr);
      else
	GLOBAL_AxiSGAccr = SG_Accr;

       Hin = ASPECTRATIO*pow(RMIN,FLARINGINDEX+1.0);

//       VKepIn = G*1.0/Rmed[0]						\
//	+ pow(ASPECTRATIO,2.0)*pow(Rmed[0],2.0*FLARINGINDEX -1.0)*
//	(-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[0]*sqrt(RMIN)*pow(Rmed[0]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[0]+Hin)))) - \
//	Rmed[0]*GLOBAL_AxiSGAccr[0];
//       VKepIn = sqrt(VKepIn);
    
        VKepOut = G*1.0/Rmed[nr-1]						\
	 + pow(ASPECTRATIO,2.0)*pow(Rmed[nr-1],2.0*FLARINGINDEX -1.0)*
	 (-SIGMASLOPE + 2.0*FLARINGINDEX - 1.0 + 0.5*Rmed[nr-1]*sqrt(RMIN)*pow(Rmed[nr-1]+Hin, -1.5)/(1.0 - sqrt(RMIN/(Rmed[nr-1]+Hin)))) - \
	 Rmed[nr-1]*GLOBAL_AxiSGAccr[IMIN+nr-1];
       VKepOut = sqrt(VKepOut);
    }
    /* ----- */
    /* i = 0 */
    /* ----- */
    if ( CPU_Rank == 0 ) {
      i = 0;
      for (j = 0; j < ns; j++) {
	l = i*ns + j;
//	vt[l] = VKepIn-Rmed[i]*OmegaFrame;
        vt[l] = vt[l+ns];
      }
    }
    /* ---------- */
    /* i = nr - 1 */
    /* ---------- */
     if(OpenInner == YES){
       if ( CPU_Rank == CPU_Highest ) {
        i = nr - 1;
         for (j = 0; j < ns; j++) {
   	l = i*ns + j;
   	vt[l] = VKepOut-Rmed[i]*OmegaFrame;
//        vt[l] = vt[l-ns];
         }
       }
	}
  }
}

void ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, step)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real step;
{
  if (OpenInner == YES) OpenBoundary (Vrad, Rho, Energy);
  if (NonReflecting == YES) {
    if (Adiabatic)
      ComputeSoundSpeed (Rho, Energy);
    NonReflectingBoundary (Vrad, Vtheta, Rho, Energy);
  }
  if (Evanescent == YES) EvanescentBoundary (Vrad, Vtheta, Rho, Energy, step);
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
     PolarGrid *vtheta;
     real domega;
{
  int i, j, l, nr, ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}
