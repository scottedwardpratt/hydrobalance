#include "msu_hydrobalance/hbcharge.h"
#include "msu_commonutils/misc.h"
#include "msu_sampler/hyper.h"
#include "msu_hydrobalance/hydro2uds.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;

CHydroBalance *CHBCharge::hb=NULL;
CB3D *CHBCharge::b3d=NULL;

void CHBCharge::Propagate(double newtau){
	double t0,tf,neweta;
	if(active==true){
		t0=tau*cosh(eta);
		neweta=rapidity-asinh((tau/newtau)*sinh(rapidity-eta));
		tf=newtau*cosh(neweta);
		//printf("rapidity=%g, t0=%g, newtau=%g, tf=%g, neweta=%g, eta=%g\n",rapidity,t0,newtau,tf,neweta,eta);
		x+=vx*(tf-t0);
		y+=vy*(tf-t0);
		if(x!=x)
			exit(1);
		tau=newtau;
		eta=neweta;
	}
}

void CHBCharge::SetV(double uxmatter,double uymatter){
	double vz,vperp,phi,vmag=1.0;
	FourVector v,umatter;
	vz=vmag*(1.0-2.0*hb->randy->ran());
	phi=2.0*PI*hb->randy->ran();
	vperp=sqrt(vmag*vmag-vz*vz);
	vx=vperp*cos(phi);
	vy=vperp*sin(phi);
	v[0]=1.0;
	v[1]=vx;
	v[2]=vy;
	v[3]=vz;
	
	//v[1]=v[2]=v[3]=0.0;
	// Now boost to lab frame
	umatter[3]=sinh(eta);
	umatter[1]=uxmatter; umatter[2]=uymatter;
	umatter[0]=sqrt(1.0+uxmatter*uxmatter+uymatter*uymatter+umatter[3]*umatter[3]);
	Misc::Boost(umatter,v);
	rapidity=atanh(v[3]/v[0]);
	vx=v[1]/v[0];
	vy=v[2]/v[0];
	//printf("CHECK: %g\n",1.0-vx*vx-vy*vy-v[3]*v[3]/(v[0]*v[0]));
}	

void CHBCharge::Print(){
	printf("Charge Info:\n");
	if(active)
		printf("active\n");
	else
		printf("dead\n");
	printf("q = (%d,%d,%d)\n",q[0],q[1],q[2]);
	printf("x = %g,y=%g, eta=%g, tau=%g\nvx=%g, vy=%g, rapidity=%g, v=%g\n",
	x,y,eta,tau,vx,vy,rapidity,sqrt(vx*vx+vy*vy+tanh(rapidity)*tanh(rapidity)));
	/*
	if(hb!=NULL){
		hb->GetIXIY(x,y,ix,iy);
		printf("ix = %d, iy=%d, T=(%g,%g)\n",ix,iy,hb->mesh->T[ix][iy],hb->newmesh->T[ix][iy]);
	}
	*/
	printf("------------------------------------------------------------\n");
}

