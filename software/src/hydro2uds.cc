#include "msu_hydrobalance/hydro2uds.h"
#include "msu_commonutils/randy.h"
#include "msu_hydrobalance/hbeos.h"
#include "msu_hydrobalance/hbcharge.h"
using namespace std;

CHydroBalance *CHydroMesh::hb=NULL;

CHydroBalance::CHydroBalance(){
	printf("WTF\n");
};

CHydroBalance::CHydroBalance(string parfilename,int ranseed){
	parmap.ReadParsFromFile(parfilename);
	Tf=parmap.getD("FREEZEOUT_TEMP",0.155);
	SIGMA0=parmap.getD("SIGMA0",0.5);
	DiffusionRatio=parmap.getD("DIFFUSION_RATIO",1.0);
	eos=new CHBEoS(&parmap);
	eos->ReadEoS_PST();
	eos->BuildMap();
	eos->GetEoSFromT_PST(Tf);
	eos->GetChiOverS_Claudia();
	eos->FillOutdDdT();
	//eos->PrintChi();
	CHydroMesh::DELTAU=parmap.getD("MESH_DELTAU",0.02);
	CHydroMesh::TAU0=parmap.getD("MESH_TAU0",0.6);
	CHydroMesh::XMIN=parmap.getD("MESH_XMIN",-13.0);
	CHydroMesh::XMAX=parmap.getD("MESH_XMAX",13.0);
	CHydroMesh::YMIN=parmap.getD("MESH_YMIN",-13.0);
	CHydroMesh::YMAX=parmap.getD("MESH_YMAX",13.0);
	CHydroMesh::NX=parmap.getI("MESH_NX",261);
	CHydroMesh::NY=parmap.getI("MESH_NY",261);
	CHydroMesh::DX=(CHydroMesh::XMAX-CHydroMesh::XMIN)/double(CHydroMesh::NX-1);
	CHydroMesh::DY=(CHydroMesh::YMAX-CHydroMesh::YMIN)/double(CHydroMesh::NY-1);
	CHydroMesh::GetDimensions(NX,NY,DX,DY,DELTAU,TAU0,XMIN,XMAX,YMIN,YMAX);
	WRITE_TRAJ=parmap.getB("HB_WRITE_TRAJ",false);
	
	NSAMPLE_HYDRO2UDS=parmap.getD("NSAMPLE_HYDRO2UDS",2);
	randy=new Crandy(ranseed);
	mesh=newmesh=oldmesh=NULL;
	Ncollisions=0;
	CHydroMesh::hb=this;
	CHBCharge::hb=this;
	tau0check=true;
	tau0readcheck=true;
	itauread=0;
	ransum=0.0;
	ranthresh=randy->ran_exp();
	idmax=0;
	MakeMeshes();
	// Variables below used for testing
	eos->T=Tf;
	eos->GetChiOverS_Claudia();
	chif.setZero();
	chifinv.setZero();
	chitot.setZero();
	chitothyper.setZero();
	chif(0,0)=chif(1,1)=eos->chill;
	chif(0,1)=chif(1,0)=eos->chiud;
	chif(0,2)=chif(1,2)=chif(2,1)=chif(2,0)=eos->chils;
	chif(2,2)=eos->chiss;
	chifinv=chif.inverse();
	printf("chif from EoS\n");
	cout << chif << endl;
	source.resize(30);
	for(int itau=0;itau<30;itau++)
		source[itau].setZero();
	ntraj=0;
}

void CHydroBalance::MakeMeshes(){
	mesh=new CHydroMesh();
	newmesh=new CHydroMesh();
	oldmesh=new CHydroMesh();
}

CHydroBalance::~CHydroBalance(){
}

void CHydroBalance::GetUxyBar(int ix,int iy,double &uxbar,double &uybar){
	uxbar=0.125*(mesh->UX[ix][iy]+mesh->UX[ix+1][iy+1]
		+mesh->UX[ix+1][iy]+mesh->UX[ix][iy+1]
			+newmesh->UX[ix][iy]+newmesh->UX[ix+1][iy+1]
				+newmesh->UX[ix+1][iy]+newmesh->UX[ix][iy+1]);
	uybar=0.125*(mesh->UY[ix][iy]+mesh->UY[ix+1][iy+1]
		+mesh->UY[ix+1][iy]+mesh->UY[ix][iy+1]
			+newmesh->UY[ix][iy]+newmesh->UY[ix+1][iy+1]
				+newmesh->UY[ix+1][iy]+newmesh->UY[ix][iy+1]);
}

void CHydroBalance::GetXYBar(int ix,int iy,double &xbar,double &ybar){
	int delix,deliy;
	double x,y;
	xbar=ybar=0.0;
	for(delix=0;delix<2;delix++){
		for(deliy=0;deliy<2;deliy++){
			mesh->GetXY(ix+delix,iy+deliy,x,y);
			xbar+=0.25*x;
			ybar+=0.25*y;
		}	
	}
}

void CHydroBalance::GetTBar(int ix,int iy,double &Tbar){
	Tbar=0.125*(mesh->T[ix][iy]+mesh->T[ix+1][iy+1]
		+mesh->T[ix+1][iy]+mesh->T[ix][iy+1]);
	Tbar+=0.125*(newmesh->T[ix][iy]+newmesh->T[ix+1][iy+1]
		+newmesh->T[ix+1][iy]+newmesh->T[ix][iy+1]);
}

void CHydroBalance::GetPiTildeBar(int ix,int iy,double &pitildexxbar,double &pitildexybar,
double &pitildeyybar){
	pitildexxbar=0.125*(mesh->pitildexx[ix][iy]+mesh->pitildexx[ix+1][iy+1]
		+mesh->pitildexx[ix+1][iy]+mesh->pitildexx[ix][iy+1]
		+newmesh->pitildexx[ix][iy]+newmesh->pitildexx[ix+1][iy+1]
			+newmesh->pitildexx[ix+1][iy]+newmesh->pitildexx[ix][iy+1]);
	pitildexybar=0.125*(mesh->pitildexy[ix][iy]+mesh->pitildexy[ix+1][iy+1]
		+mesh->pitildexy[ix+1][iy]+mesh->pitildexy[ix][iy+1]
		+newmesh->pitildexy[ix][iy]+newmesh->pitildexy[ix+1][iy+1]
			+newmesh->pitildexy[ix+1][iy]+newmesh->pitildexy[ix][iy+1]);
	pitildeyybar=0.125*(mesh->pitildeyy[ix][iy]+mesh->pitildeyy[ix+1][iy+1]
		+mesh->pitildeyy[ix+1][iy]+mesh->pitildeyy[ix][iy+1]
		+newmesh->pitildeyy[ix][iy]+newmesh->pitildeyy[ix+1][iy+1]
			+newmesh->pitildeyy[ix+1][iy]+newmesh->pitildeyy[ix][iy+1]);
}

void CHydroBalance::MakeCharges(){
	int ix,iy,sign,a,b,itau;
	double DQll,DQud,DQls,DQss,g1,g2;
	Eigen::Matrix3d DQ(3,3);
	CHBCharge *charge1,*charge2;
	for(ix=1;ix<mesh->NX-1;ix++){
		for(iy=1;iy<=mesh->NY-1;iy++){
			if(mesh->T[ix][iy]>Tf || newmesh->T[ix][iy]>Tf){
				if(tau0check){
					CalcDQ0(ix,iy,DQll,DQud,DQls,DQss);
				}
				else{
					CalcDQ(ix,iy,DQll,DQud,DQls,DQss);
				}
				//printf("(%3d,%3d): %10.3f %10.3f %10.3f %10.3f T=%g\n",ix,iy,DQll,DQud,DQls,DQss,newmesh->T[ix][iy]);
				DQ(0,0)=DQ(1,1)=DQll;
				DQ(0,1)=DQ(1,0)=DQud;
				DQ(0,2)=DQ(1,2)=DQ(2,1)=DQ(2,0)=DQls;
				DQ(2,2)=DQss;
				for(a=0;a<3;a++){
					for(b=0;b<3;b++){
						ransum+=NSAMPLE_HYDRO2UDS*fabs(DQ(a,b));
						chitot(a,b)+=DQ(a,b);
						while(ransum>ranthresh){
							ranthresh+=randy->ran_exp();
							charge1=new CHBCharge;
							charge2=new CHBCharge;
							mesh->GetXY(ix,iy,charge1->x,charge1->y);
							charge1->x=charge2->x=XMIN+(ix+randy->ran())*DX;
							charge1->y=charge2->y=YMIN+(iy+randy->ran())*DY;
							charge1->tau=charge2->tau=mesh->tau;
							charge1->eta=charge2->eta=0.0;
							if(tau0check){
								randy->ran_gauss2(g1,g2);
								charge1->eta=SIGMA0*g1;
								charge2->eta=SIGMA0*g2;
								/*
								if(abs(charge1->q[2])==1 && (abs(charge2->q[2])==1)){
									printf("tau=%g, eta=(%g,%g), q=(%d,%d)\n",mesh->tau,charge1->eta,charge2->eta,charge1->q[2],charge2->q[2]);
								}
								*/
							}
							charge1->active=charge2->active=true;
							charge1->weight=charge2->weight=1.0;
							charge1->q[0]=charge1->q[1]=charge1->q[2]=0;
							charge2->q[0]=charge2->q[1]=charge2->q[2]=0;
							charge1->SetV(mesh->UX[ix][iy],mesh->UY[ix][iy]);
							charge2->SetV(mesh->UX[ix][iy],mesh->UY[ix][iy]);
							sign=1;
							if(randy->ran()<0.5)
								sign=-1;
							charge1->q[a]=sign;
							if(DQ(a,b)>0.0)
								sign=-sign;
							charge2->q[b]=sign;
							//chitot(a,b)+=charge1->q[a]*charge2->q[b]
								//	+charge2->q[a]*charge1->q[b];
							cmap.insert(pairic(idmax,charge1));
							idmax+=1;
							cmap.insert(pairic(idmax,charge2));
							if(WRITE_TRAJ && tau0check && randy->ran()<0.2){
								if(charge1->q[2]!=0 && charge2->q[2]!=0){ // write for ss CF
									charge1->trajinfo=new CTrajInfo(ntraj);
									ntraj+=1;
									charge2->trajinfo=new CTrajInfo(ntraj);
									ntraj+=1;
									charge1->addtraj();
									charge2->addtraj();
								}
							}
							idmax+=1;
							itau=floorl(mesh->tau/0.5);
							if(itau<30)
								source[itau](a,b)-=charge1->q[a]*charge2->q[b];

						}
					}
				}
			}
		}
	}
	tau0check=false;
}

void CHydroBalance::PropagateCharges(){
	mapic::iterator it,oldit,its;
	double dTdt,dTdx,dTdy,u0;
	bool GGTt,GGTx,GGTy;
	Chyper *hyper;
	CHBCharge *charge;
	int ix,iy,ic,id;
	double newtau=newmesh->tau;
	it=cmap.begin();
	its=it;
	ic=0;
	while(it!=cmap.end()){
		charge=it->second;
		if(ic!=0){
			its=it;
			its--;
		}
		ic+=1;
		id=it->first;
		charge->Propagate(newtau);
		mesh->GetIXIY_lower(charge->x,charge->y,ix,iy);
		if(ix<=0 || iy<=0 || ix>=CHydroMesh::NX-1 || iy>=CHydroMesh::NY-1){
			printf("CHydroBalance::PropagateCharges() disaster, ix=%d, iy=%d\n",ix,iy);
			exit(1);
		}
		if(!(charge->active) || (charge->active && newmesh->T[ix][iy]<Tf)){
			emap.insert(pairic(id,charge));
			hyper=&(charge->hyper);
			if(WRITE_TRAJ)
				charge->addtraj();
			GetGradT(ix,iy,dTdt,dTdx,dTdy,GGTt,GGTx,GGTy);
			GetUxyBar(ix,iy,hyper->u[1],hyper->u[2]);
			GetXYBar(ix,iy,hyper->r[1],hyper->r[2]);
			hyper->tau=newmesh->tau;
			GetPiTildeBar(ix,iy,hyper->pitilde[1][1],hyper->pitilde[1][2],hyper->pitilde[2][2]);
			hyper->dOmega[0]=-dTdt;
			hyper->dOmega[1]=-dTdx;
			hyper->dOmega[2]=-dTdy;
			u0=sqrt(1.0+hyper->u[1]*hyper->u[1]+hyper->u[2]*hyper->u[2]);
			hyper->u[0]=u0;
			hyper->udotdOmega=(u0*hyper->dOmega[0]-hyper->u[1]*hyper->dOmega[1]-hyper->u[2]*hyper->dOmega[2]);
			charge->active=false;
			oldit=it;
			++it;
			cmap.erase(oldit);	
		}
		else{
			++it;
		}
	}
	it=cmap.begin();
	while(it!=cmap.end()){
		charge=it->second;
		if(!(charge->active)){
			printf("no dead particles should survive\n");
			exit(1);
		}
		++it;
	}
}

void CHydroBalance::ScatterCharges(){
	mapic::iterator it;
	double u0,ux,uy,D;
	CHBCharge *charge;
	int ix,iy;
	it=cmap.begin();
	while(it!=cmap.end()){
		charge=it->second;
		mesh->GetIXIY(charge->x,charge->y,ix,iy);
		D=CHBEoS::GetD(newmesh->T[ix][iy]);
		D*=DiffusionRatio;
		ux=newmesh->UX[ix][iy];
		uy=newmesh->UY[ix][iy];
		u0=sqrt(1.0+ux*ux+uy*uy);
		ransum+=DELTAU/(6.0*D*u0);
		while(ransum>ranthresh){
			ranthresh+=randy->ran_exp();
			charge->SetV(ux,uy);
			if(WRITE_TRAJ)
					charge->addtraj();
			Ncollisions+=1;
		}
		++it;
	}
}

void CHydroBalance::CalcDQ(int ix,int iy,double &DQll,
double &DQud,double &DQls,double &DQss){
	double d4x,s0,sx,sy;
	CHBEoS eos222[2][2][2];
	int j0,jx,jy;
	double ux[2][2][2],uy[2][2][2],u0[2][2][2];
	CHydroMesh *mptr;
	
	for(j0=0;j0<2;j0++){
		if(j0==0)
			mptr=mesh;
		else
			mptr=newmesh;
		for(jx=0;jx<2;jx++){
			for(jy=0;jy<2;jy++){
				ux[j0][jx][jy]=mptr->UX[ix+jx][iy+jy];
				uy[j0][jx][jy]=mptr->UY[ix+jx][iy+jy];
				u0[j0][jx][jy]=sqrt(1.0
					+ux[j0][jx][jy]*ux[j0][jx][jy]+uy[j0][jx][jy]*uy[j0][jx][jy]);
				eos222[j0][jx][jy].GetEoSFromT_PST(mptr->T[ix+jx][iy+jy]);
				if(mptr->T[ix+jx][iy+jy]>199.99 && mptr->T[ix+jx][iy+jy]<200.1){
					printf("T=%g, s=%g\n",mptr->T[ix+jx][iy+jy],eos222[j0][ix+jx][iy+jy].s);
				}
				eos222[j0][jx][jy].GetChiOverS_Claudia();
			}
		}
	}
	DQll=DQud=DQls=DQss=0.0;
	for(j0=0;j0<2;j0++){
		if(j0==0){
			s0=-0.25;
			d4x=DELTAU*DX*DY*mesh->tau;
		}
		else{
			s0=0.25;
			d4x=DELTAU*DX*DY*newmesh->tau;
		}
		for(jx=0;jx<2;jx++){
			if(jx==0)
				sx=-0.25;
			else
				sx=0.25;
			for(jy=0;jy<2;jy++){
				if(jy==0)
					sy=-0.25;
				else
					sy=0.25;
				
				DQll+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chill/DELTAU;
				DQud+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chiud/DELTAU;
				DQls+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chils/DELTAU;
				DQss+=d4x*s0*u0[j0][jx][jy]*eos222[j0][jx][jy].chiss/DELTAU;
				
				DQll+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chill/DX;
				DQud+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chiud/DX;
				DQls+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chils/DX;
				DQss+=d4x*sx*ux[j0][jx][jy]*eos222[j0][jx][jy].chiss/DX;
				
				DQll+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chill/DY;
				DQud+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chiud/DY;
				DQls+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chils/DY;
				DQss+=d4x*sy*uy[j0][jx][jy]*eos222[j0][jx][jy].chiss/DY;
								
			}
		}
	}
}

void CHydroBalance::CalcDQ0(int ix,int iy,double &DQll,double &DQud,
double &DQls,double &DQss){
	double u0,ux,uy,d3x=DX*DY*newmesh->tau,T;
	GetUxyBar(ix,iy,ux,uy);
	u0=sqrt(1.0+ux*ux+uy*uy);
	GetTBar(ix,iy,T);
	if(T>Tf){
		eos->GetEoSFromT_PST(T);
		eos->GetChiOverS_Claudia();
		DQll=d3x*u0*eos->chill;
		DQud=d3x*u0*eos->chiud;
		DQls=d3x*u0*eos->chils;
		DQss=d3x*u0*eos->chiss;
	}
	else{
		DQll=DQud=DQls=DQss=0.0;
	}
}

void CHydroBalance::SwapMeshes(){
	CHydroMesh *swap;
	swap=oldmesh;
	oldmesh=mesh;
	mesh=newmesh;
	newmesh=swap;
}

void CHydroBalance::Reset(){
	Ncollisions=0;
	itauread=0;
	idmax=0;
	tau0readcheck=tau0check=true;
	oldmesh->tau=mesh->tau=newmesh->tau=CHydroMesh::TAU0;
}


