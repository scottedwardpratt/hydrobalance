#include "msu_hydrobalance/hydro2uds.h"

using namespace std;
using namespace NMSUPratt;

int CHydroMesh::NX=0;
int CHydroMesh::NY=0;
double CHydroMesh::DX=0.0;
double CHydroMesh::DY=0.0;
double CHydroMesh::XMIN=0.0;
double CHydroMesh::XMAX=0.0;
double CHydroMesh::YMIN=0.0;
double CHydroMesh::YMAX=0.0;
double CHydroMesh::DELTAU=0.0;
double CHydroMesh::TAU0=0.0;

CHydroMesh::CHydroMesh(){
	int ix,iy;
	T=new double*[NX];
	DQ=new double*[NX];
	UX=new double *[NX];
	UY=new double *[NX];
	pitildexx=new double *[NX];
	pitildexy=new double *[NX];
	pitildeyy=new double *[NX];
	for(ix=0;ix<NX;ix++){
		T[ix]=new double[NY];
		DQ[ix]=new double[NY];
		UX[ix]=new double[NY];
		UY[ix]=new double[NY];
		pitildexx[ix]=new double[NY];
		pitildexy[ix]=new double[NY];
		pitildeyy[ix]=new double[NY];
		for(iy=0;iy<NY;iy++){
			T[ix][iy]=DQ[ix][iy]=UX[ix][iy]=UY[ix][iy]
				=pitildexx[ix][iy]=pitildexy[ix][iy]=pitildeyy[ix][iy]=0.0;
		}
	}
	tau=TAU0;
}

void CHydroMesh::GetXY(int ix,int iy,double &x,double &y){
	x=XMIN+ix*DX;
	y=YMIN+iy*DY;
}

void CHydroMesh::GetIXIY(double x,double y,int &ix,int &iy){
	ix=lrint((x-XMIN)/DX);
	iy=lrint((y-YMIN)/DY);
	//if(ix==0 || ix==NX || iy==0 || iy==NY){
		//printf("At edge of mesh\n");
	//}
}

void CHydroMesh::GetIXIY_lower(double x,double y,int &ix,int &iy){
	ix=floorl((x-XMIN)/DX);
	iy=floorl((y-YMIN)/DY);
	//if(ix==0 || ix==NX || iy==0 || iy==NY){
		//printf("At edge of mesh\n");
	//}
}