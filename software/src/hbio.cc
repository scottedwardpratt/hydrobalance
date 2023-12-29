#include "msu_hydrobalance/hydro2uds.h"
#include "msu_hydrobalance/hbcharge.h"
#include "msu_sampler/hyper.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

bool CHydroBalance::ReadOSCAR(CHydroMesh *hydromesh){
	double r,x,y,rmax=0.0,highestT,biggestU,ur;
	bool keepgoing=true;
	int olditau,ix,iy,iline,flag,alpha;
	char dummy[300];
	double **pi,**pitilde;
	FourVector u;
	double e,p,t,vx,vy,pi00,pi01,pi02,pi11,pi12,pi22,pi33,Pi;
	pi=new double*[4];
	pitilde=new double*[4];
	for(alpha=0;alpha<4;alpha++){
		pi[alpha]=new double[4];
		pitilde[alpha]=new double[4];
	}
	highestT=biggestU=0.0;
	if(tau0readcheck){
		oscar_filename="../hydrodata/"+qualifier+"/"+parmap.getS("HYDRODATA_FILENAME","OSCAR2008H.dat");
		printf("filename=%s\n",oscar_filename.c_str());
		fptr_oscar=fopen(oscar_filename.c_str(),"r");
		for(iline=0;iline<14;iline++){
			fgets(dummy,300,fptr_oscar);
			//printf("%s",dummy);
		}
	}
	
	for(ix=0;ix<NX;ix++)
		for(iy=0;iy<NY;iy++)
			hydromesh->T[ix][iy]=0.0;
	olditau=itauread;
	hydromesh->tau=TAU0+itauread*DELTAU;
	if(tau0readcheck)
		fscanf(fptr_oscar,"%d",&itauread);
	if(feof(fptr_oscar)){
		keepgoing=false;
		fclose(fptr_oscar);
	}
	else{
		while(itauread==olditau && !feof(fptr_oscar)){
			fscanf(fptr_oscar,"%d %d %lf %lf %lf %d %lf %lf",&ix,&iy,&e,&p,&t,&flag,&vx,&vy);
			//printf("%3d: %3d %3d: %8.6f %8.5f %8.5f\n",itauread,ix,iy,t,e,p);
			fscanf(fptr_oscar,"%lf %lf %lf %lf %lf %lf %lf %lf",
			&pi00,&pi01,&pi02,&pi11,&pi12,&pi22,&pi33,&Pi);
			//printf("pi00=%g, pi01=%g, pi11=%g, pi22=%g,pi33=%g, Tr pi=%g\n",
			//pi00,pi01,pi11,pi22,pi33,pi00-pi11-pi22-pi33);
			//fgets(dummy,300,fptr_oscar);
			pi[0][0]=pi00;
			pi[0][1]=pi[1][0]=pi01;
			pi[0][2]=pi[2][0]=pi02;
			pi[0][3]=pi[3][0]=0.0;
			pi[1][1]=pi11;
			pi[1][2]=pi[2][1]=pi12;
			pi[1][3]=pi[3][1]=0.0;
			pi[2][2]=pi22;
			pi[2][3]=pi[3][2]=0.0;
			pi[3][3]=pi33;
			
			u[0]=1.0/sqrt(1.0-vx*vx-vy*vy);
			u[1]=u[0]*vx;
			u[2]=u[0]*vy;
			u[3]=0.0;
			hydromesh->GetXY(ix,iy,x,y);
			r=sqrt(x*x+y*y);
			if(r>rmax && t>=Tf)
				rmax=r;
			Misc::BoostToCM(u,pi,pitilde);
			hydromesh->pitildexx[ix][iy]=pitilde[1][1];
			hydromesh->pitildexy[ix][iy]=pitilde[1][2];
			hydromesh->pitildeyy[ix][iy]=pitilde[2][2];
		
			fgets(dummy,300,fptr_oscar);
			hydromesh->T[ix][iy]=t;
			hydromesh->UX[ix][iy]=u[1];
			hydromesh->UY[ix][iy]=u[2];
			ur=sqrt(u[1]*u[1]+u[2]*u[2]);
			if(ur>biggestU && t>Tf)
				biggestU=ur;
			olditau=itauread;
			fscanf(fptr_oscar,"%d",&itauread);
			if(t>highestT)
				highestT=t;
		}
		tau0readcheck=false;
		if(highestT<Tf || feof(fptr_oscar)){
			keepgoing=false;
			fclose(fptr_oscar);
		}
	}
	if(fabs(lrint(hydromesh->tau)-hydromesh->tau)<0.001)
		printf("highestT=%g, biggestU=%g\n",highestT,biggestU);
	if(!keepgoing)
		printf("XXXXXXXXXX almost ready to quit\n");
	for(alpha=0;alpha<4;alpha++){
		delete pi[alpha];
		delete pitilde[alpha];
	}
	delete pi;
	delete pitilde;
	
	return keepgoing;
}

void CHydroBalance::WriteCharges(){
	string dirname="udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/"+parmap.getS("CHARGESINFO_FILENAME","uds.dat");
	printf("writing charges to %s\n",filename.c_str());
	mapic::iterator it;
	CHBCharge *charge;
	Chyper *hyper;
	int balanceID;
	unsigned int icharge;
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#  id   u  d  s      weight        tau           eta           x             y\n");
	it=emap.begin();
	for(it=emap.begin();it!=emap.end();++it){
		balanceID=it->first;
		//("writing, balanceID=%d\n",balanceID);
		charge=it->second;
		hyper=&(charge->hyper);

		fprintf(fptr,"%6d %2d %2d %2d %15.9f %15.9f %15.9f %15.9f %15.9f %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e %15.9e\n",
		balanceID,charge->q[0],charge->q[1],charge->q[2],charge->weight,charge->tau,charge->eta,
		charge->x,charge->y,hyper->u[1],hyper->u[2],hyper->dOmega[0],hyper->dOmega[1],
		hyper->dOmega[2],hyper->pitilde[1][1],hyper->pitilde[2][2],hyper->pitilde[1][2]);

		if(WRITE_TRAJ){
			if(charge->trajinfo!=NULL){
				for(icharge=0;icharge<charge->trajinfo->x.size();icharge++){
					printf("writing trajectory, icharge=%d\n",icharge);
					fprintf(charge->trajinfo->fptr,"%8.5f %8.5f %8.5f %8.5f\n",
					charge->trajinfo->x[icharge],charge->trajinfo->y[icharge],charge->trajinfo->eta[icharge],charge->trajinfo->tau[icharge]);
				}
				fclose(charge->trajinfo->fptr);
			}
		}
	}
	fclose(fptr);
	printf("Ncolls/charge=%g\n",2.0*double(Ncollisions)/emap.size());
	
	/*
	// for testing
	Eigen::Matrix3d chitotcharges;
	chitotcharges.setZero();
	int a,b;
	it=emap.begin();
	CHBCharge *charge1,*charge2;
	while(it!=emap.end()){
		charge1=it->second;
		++it;
		charge2=it->second;
		for(a=0;a<3;a++){
			for(b=0;b<3;b++){
				chitotcharges(a,b)+=charge1->q[a]*charge2->q[b]/double(NSAMPLE_HYDRO2UDS);
			}
		}
		++it;
	}
	printf("ChiTot From Charges in Final State\n");
	cout << chitotcharges << endl;*/
}

void CHydroBalance::ClearCharges(){
	mapic::iterator it;
	it=emap.begin();
	CHBCharge *charge1,*charge2;
	while(it!=emap.end()){
		charge1=it->second;
		++it;
		charge2=it->second;
		delete charge1;
		delete charge2;
		++it;
	}
	emap.clear();
}

void CHydroBalance::WriteSource(){
	// Note: DTAU for source mesh was 0.5 fm/c
	int a,b,itau;
	string dirname="udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/udssource.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	for(itau=0;itau<30;itau++){
		fprintf(fptr,"%5.2f ",(0.5+itau)*0.5);
		for(a=0;a<3;a++){
			for(b=a;b<3;b++)
				fprintf(fptr,"%7.0f ",0.5*(source[itau](a,b)+source[itau](b,a))/(0.5*NSAMPLE_HYDRO2UDS));
		}
		fprintf(fptr,"\n");
	}
	fclose(fptr);
}

void CHydroBalance::WriteFinalCF(){
	int ieta,a,b,c;
	int Netabins=parmap.getD("CF_NETABINS",50);
	double Deta=parmap.getD("CF_DETA",0.1);
	double eta,Z=NSAMPLE_HYDRO2UDS*Deta;
	Eigen::Matrix3d *cf=new Eigen::Matrix3d[Netabins];
	for(ieta=0;ieta<Netabins;ieta++)
		cf[ieta].setZero(3,3);
	mapic::iterator it;
	CHBCharge *charge1,*charge2;
	it=emap.begin();
	while(it!=emap.end()){
		charge1=it->second;
		++it;
		charge2=it->second;
		++it;
		eta=fabs(charge1->eta-charge2->eta);
		ieta=floorl(eta/Deta);
		if(ieta<Netabins){
			a=b=-1;
			for(c=0;c<3;c++){
				if(charge1->q[c]!=0)
					a=c;
				if(charge2->q[c]!=0)
					b=c;
			}
			cf[ieta](a,b)+=0.5*charge1->q[a]*charge2->q[b];
			cf[ieta](b,a)+=0.5*charge1->q[a]*charge2->q[b];
		}
	}
	
	string dirname="udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/cf_uds.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Deta       uu            ud          us            ss\n");
	for(ieta=0;ieta<Netabins;ieta++){
		fprintf(fptr,"%4.2f %11.5e %11.5e %11.5e %11.5e\n",
		(0.5+ieta)*Deta,cf[ieta](0,0)/Z,cf[ieta](0,1)/Z,cf[ieta](0,2)/Z,cf[ieta](2,2)/Z);
	}
	fclose(fptr);
	delete [] cf;
}

void CHydroBalance::WriteHyper(){
	string dirname="udsdata/"+qualifier;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	string filename=dirname+"/"+parmap.getS("HYPERDATA_FILENAME","hyper.dat");
	printf("writing hyper info to %s\n",filename.c_str());
	Chyper *hyper;
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"#Tf=%g\n",Tf);
	fprintf(fptr,"#   tau        x            y           Ux              Uy         dOmega0       dOmegaX       dOmegaY\n");
	list<Chyper *>::iterator it;
	for(it=hyperlist.begin();it!=hyperlist.end();++it){
		hyper=*it;
		fprintf(fptr,"%13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e %13.7e\n",
		hyper->tau,hyper->r[1],hyper->r[2],hyper->u[1],hyper->u[2],hyper->dOmega[0],hyper->dOmega[1],hyper->dOmega[2],hyper->pitilde[1][2],hyper->pitilde[2][2],hyper->pitilde[1][2]);
	}
	printf("Wrote %d hyper-elements\n",int(hyperlist.size()));
	fclose(fptr);
}
