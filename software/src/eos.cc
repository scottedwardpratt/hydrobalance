#include "msu_hydro2uds/eos.h"
#include "msu_sampler/resonances.h"
#include "msu_commonutils/parametermap.h"
#include "msu_commonutils/constants.h"
#include "msu_sampler/eos.h"

using namespace std;

vector<double> CEoS::epsilon_PST,CEoS::P_PST,CEoS::s_PST,CEoS::T_PST;
vector<double> CEoS::epsilon_claudia,CEoS::P_claudia,CEoS::s_claudia,CEoS::T_claudia;
vector<double> CEoS::twopiTD,CEoS::Tdiff;
vector<double> CEoS::chill_overs_claudia,CEoS::chiud_overs_claudia,CEoS::chils_overs_claudia,CEoS::chiss_overs_claudia;
vector<double> CEoS::chill_HSC,CEoS::chiud_HSC,CEoS::chils_HSC,CEoS::chiss_HSC;
vector<double> CEoS::dDdT;
mapdi CEoS::etmap;

CresList *CEoS::reslist=NULL;

CEoS::CEoS(){	
}

CEoS::CEoS(CparameterMap *parmapset){
	parmap=parmapset;
	reslist=new CresList(parmap);
	ReadDiffusionData();
	ReadChiData_Claudia();
	//FillOutdDdT();
};

void CEoS::ReadDiffusionData(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","../latticedata");
	string filename=dirname+"/diffusion.dat";
	char dummy[100];
	char voldummy[100];
	int ntaudummy;
	double errsysdummy,errsysstatdummy,t,td;
	FILE *fptr=fopen(filename.c_str(),"r");
	fgets(dummy,100,fptr);
	//printf("%s\n",dummy);
	fscanf(fptr,"%s",voldummy);
	while(!feof(fptr)){
		fscanf(fptr,"%lf %d %lf %lf %lf",&t,&ntaudummy,&td,&errsysdummy,&errsysstatdummy);
		Tdiff.push_back(t*0.001);
		twopiTD.push_back(td);
		fscanf(fptr,"%s",voldummy);
	}	
	fclose(fptr);
}

void CEoS::FillOutdDdT(){
	int ndata=epsilon_PST.size();
	dDdT.resize(ndata);
	Eigen::Matrix3d chi0,chi1,chi2,chiinv0,chiinv1,chiinv2,Dmat;
	chi0.setZero(); chi1.setZero(); chi2.setZero(); chiinv0.setZero(); chiinv1.setZero(); chiinv2.setZero();
	Dmat.setZero();
	int ie;
	double T0,T1,T2;
	for(ie=1;ie<ndata-1;ie++){
		T0=T_PST[ie-1];
		T1=T_PST[ie];
		T2=T_PST[ie+1];
		T=T0;
		GetChiOverS_Claudia();  // check whether this gives Chi or ChiOverS
		chi0(0,0)=chi0(1,1)=chill;
		chi0(0,1)=chi0(1,0)=chiud;
		chi0(2,2)=chiss;
		chi0(0,2)=chi0(1,2)=chi0(2,0)=chi0(2,1)=chils;
		chiinv0=chi0.inverse();
		T=T1;
		GetChiOverS_Claudia();
		chi1(0,0)=chi1(1,1)=chill;
		chi1(0,1)=chi1(1,0)=chiud;
		chi1(2,2)=chiss;
		chi1(0,2)=chi1(1,2)=chi1(2,0)=chi1(2,1)=chils;
		chiinv1=chi1.inverse();
		T=T2;
		GetChiOverS_Claudia();
		chi2(0,0)=chi2(1,1)=chill;
		chi2(0,1)=chi2(1,0)=chiud;
		chi2(2,2)=chiss;
		chi2(0,2)=chi2(1,2)=chi2(2,0)=chi2(2,1)=chils;
		chiinv2=chi2.inverse();
		Dmat=(chiinv2-chiinv0)/(T2-T0);
		Dmat=GetD(T1)*chi1*Dmat;
		/*
		printf("---- T=(%g,%g,%g) ----\n",T0,T1,T2);
		chi0.print("chi0");
		chiinv0.print("chiinv0");
		chi1.print("chi1");
		chiinv1.print("chiinv1");
		chi2.print("chi2:");
		chiinv2.print("chiinv2");
		Dmat.print("Dmat\n");
		*/
	}
}

double CEoS::GetD(double T){
	int i0;
	double result,delT;
	if(T>=Tdiff[0] && T<=Tdiff[Tdiff.size()-1]){
		i0=0;
		while(T>Tdiff[i0+1]){
			i0+=1;
		}
		delT=Tdiff[i0+1]-Tdiff[i0];
		result=((T-Tdiff[i0])*twopiTD[i0+1]/delT)+(Tdiff[i0+1]-T)*twopiTD[i0]/delT;
	}
	else if(T<Tdiff[0]){
		delT=Tdiff[1]-Tdiff[0];
		result=twopiTD[0]-(Tdiff[0]-T)*(twopiTD[1]-twopiTD[0])/delT;
	}
	else{
		i0=Tdiff.size()-1;
		delT=Tdiff[i0]-Tdiff[i0-1];
 		result=twopiTD[i0]+(T-Tdiff[i0])*(twopiTD[i0]-twopiTD[i0-1])/delT;
	}
	result=result*HBARC_GEV/(2.0*PI*T);
	//result*=0.001;
	return result;
}

void CEoS::ReadChiData_HSC(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","../latticedata");
	string filename;
	double error,chi0;
	int idata;
	const int ndata=7;
	double chiB[ndata],chiI[ndata],chiQ[ndata],chiSS[ndata],chiLL[ndata],Tarray[ndata];
	FILE *fptr;
	filename=dirname+"/chi-B.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&Tarray[idata],&chiB[idata],&error);
	}
	fclose(fptr);
	filename=dirname+"/chi-I.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&Tarray[idata],&chiI[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-LL.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&Tarray[idata],&chiLL[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-Q.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&Tarray[idata],&chiQ[idata],&error);
	}
	fclose(fptr);
	fclose(fptr);
	filename=dirname+"/chi-SS.dat";
	fptr=fopen(filename.c_str(),"r");
	for(idata=0;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf",&Tarray[idata],&chiSS[idata],&error);
	}
	fclose(fptr);
	
	chill_HSC.resize(ndata);
	chils_HSC.resize(ndata);
	chiss_HSC.resize(ndata);
	chiud_HSC.resize(ndata);
	for(idata=0;idata<ndata;idata++){
		chi0=pow(Tarray[idata]/HBARC,3)*(pow(PI,4)/90.0)*(4.0/(2.0*PI*PI))*2.0;
		chiLL[idata]*=chi0;
		chiB[idata]*=chi0/3.0;
		chiSS[idata]*=chi0;
		chiI[idata]*=0.5*chi0;
		chiQ[idata]*=2.0*chi0/3.0;
		chill_HSC[idata]=chiLL[idata];
		chiss_HSC[idata]=chiSS[idata];
		chiud_HSC[idata]=chiLL[idata]-2.0*chiI[idata];
		chils_HSC[idata]=2.25*(chiB[idata]-(2.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]-(2.0/9.0)*chiud_HSC[idata]);
		//chils[idata]=-4.5*(chiQ[idata]-(5.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]+(4.0/9.0)*chiud[idata]);
		printf("-------  T=%g MeV --------\n",Tarray[idata]);
		printf("%8.5f %8.5f %8.5f\n",chill_HSC[idata],chiud_HSC[idata],chils_HSC[idata]);
		printf("%8.5f %8.5f %8.5f\n",chiud_HSC[idata],chill_HSC[idata],chils_HSC[idata]);
		printf("%8.5f %8.5f %8.5f\n",chils_HSC[idata],chils_HSC[idata],chiss_HSC[idata]);
		/*
		chils_HSC[idata]=2.25*(chiB[idata]-(2.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]-(2.0/9.0)*chiud_HSC[idata]);
		printf("T=%g: chils_HSC[%d]=%g\n",Tarray[idata],idata,chils_HSC[idata]);
		chils_HSC[idata]=-4.5*(chiQ[idata]-(5.0/9.0)*chiLL[idata]-(1.0/9.0)*chiSS[idata]+(4.0/9.0)*chiud[idata]);
		printf("T=%g: chils_HSC[%d]=%g\n",Tarray[idata],idata,chils_HSC[idata]);
		printf("--------\n");
		*/
		
	}
}

void CEoS::ReadChiData_Claudia(){
	string dirname=parmap->getS("LATTICEDATA_DIRNAME","latticedata");
	string filename;
	int idata;
	const int ndata=81;
	char dummy[100];
	FILE *fptr;
	// You will read in chi/s, not chi
	filename=dirname+"/chi.dat";
	fptr=fopen(filename.c_str(),"r");
	T_claudia.resize(ndata);
	chill_overs_claudia.resize(ndata);
	chils_overs_claudia.resize(ndata);
	chiss_overs_claudia.resize(ndata);
	chiud_overs_claudia.resize(ndata);
	fgets(dummy,100,fptr);
	for(idata=4;idata<ndata;idata++){
		fscanf(fptr,"%lf %lf %lf %lf %lf",
		&T_claudia[idata],&chill_overs_claudia[idata],&chiud_overs_claudia[idata],&chils_overs_claudia[idata],&chiss_overs_claudia[idata]);
		T_claudia[idata]*=0.001;
	}
	fclose(fptr);
}

void CEoS::GetChiOverS_Claudia(){
	double delT=5.0,Tmax=0.4,w0,w1,Tm=1000.0*T;
	const int ndata=81;
	int iT0;
	printf("In CEoS::GetChiOverS_Claudia, T=%g, Tm=%g\n",T,Tm);
	if(Tm<20){
		chill=chill_overs_claudia[5]*Tm*s/20.0;
		chiud=chiud_overs_claudia[5]*Tm*s/20.0;
		chils=chils_overs_claudia[5]*Tm*s/20.0;
		chiss=chiss_overs_claudia[5]*Tm*s/20.0;
	}
	else if(Tm>=Tmax){
		chill=chill_overs_claudia[ndata-1]*s;
		chiud=chiud_overs_claudia[ndata-1]*s;
		chils=chils_overs_claudia[ndata-1]*s;
		chiss=chiss_overs_claudia[ndata-1]*s;
	}
	else{
		iT0=lrint(floor(Tm/delT));
		if(iT0>=ndata-1){
			printf("iT0 out of range in CEoS::GetChiOverS_Claudia()\n");
			exit(1);
		}
		w1=(Tm-delT*iT0)/delT;
		w0=1.0-w1;
		chill=chill_overs_claudia[iT0]*s*w0+chill_overs_claudia[iT0+1]*s*w1;
		chiud=chiud_overs_claudia[iT0]*s*w0+chiud_overs_claudia[iT0+1]*s*w1;
		chils=chils_overs_claudia[iT0]*s*w0+chils_overs_claudia[iT0+1]*s*w1;
		chiss=chiss_overs_claudia[iT0]*s*w0+chiss_overs_claudia[iT0+1]*s*w1;
	}
	
}

void CEoS::CalcEoS_PST(){
	Eigen::Matrix3d chi;
	Eigen::Matrix3d sigma;
	string filename
		=parmap->getS("EOS_PSTDATA_FILENAME","../eos/EOS_tables/EOS_PST.dat");
	FILE *fptr=fopen(filename.c_str(),"r");
	FILE *fptrh=fopen("hadroneos.dat","w");
	double eps;
	const int ndata=155500;
	epsilon_PST.resize(ndata);
	P_PST.resize(ndata);
	s_PST.resize(ndata);
	T_PST.resize(ndata);
	
	int ie=0,ne;
	fscanf(fptr,"%lf ",&eps);
	while(!feof(fptr)){
		epsilon_PST[ie]=eps;
		fscanf(fptr,"%lf %lf %lf",&P_PST[ie],&s_PST[ie],&T_PST[ie]);
		//printf("epsilon=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
		ie+=1;
		fscanf(fptr,"%lf ",&eps);
	}
	
	double TH=0.150, TQGP=0.160;
	double TT,Ph,eh,nh,sh,w;
	int nres=reslist->resmap.size();
	vector<double> density;
	vector<double> maxweight;
	density.resize(nres);
	maxweight.resize(nres);
	ne=ie;
	for(ie=0;ie<ne;ie++){
		if(T_PST[ie]<TQGP){
			TT=T_PST[ie];
			MSU_EOS::CalcEoSandTransportCoefficients(TT,reslist,eh,Ph,nh,density,chi,sigma);
			sh=(Ph+eh)/TT;
			if(T<TH){
				P_PST[ie]=Ph;
				s_PST[ie]=sh;
				epsilon_PST[ie]=eh;
			}
			else{
				w=(TQGP-TT)/(TQGP-TH);
				P_PST[ie]=w*Ph+(1.0-w)*P_PST[ie];
				s_PST[ie]=w*sh+(1.0-w)*s_PST[ie];
				epsilon_PST[ie]=w*eh+(1.0-w)*epsilon_PST[ie];
				//printf("Tpst=%g, Ppst=%g, spst=%g, epst=%g\n",
				//T_PST[ie],P_PST[ie],s_PST[ie],epsilon_PST[ie]);
				//printf("T=%g, Ph=%g, sh=%g, eh=%g\n",T,Ph,sh,eh);
			}
		}
		fprintf(fptrh,"%13.7e %13.7e %13.7e %13.7e\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
	}
	//printf("ie=%d =? %d\n",ie,ndata);
	//printf("epsilon=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie-1],P_PST[ie-1],s_PST[ie-1],T_PST[ie-1]);
	fclose(fptr);
	fclose(fptrh);
	
}

void CEoS::ReadEoS_PST(){
	string filename
		=parmap->getS("EOS_PSTDATA_FILENAME","../eos/EOS_tables/EOS_PST.dat");
	FILE *fptr=fopen(filename.c_str(),"r");
	double eps;
	const int ndata=155500;
	epsilon_PST.resize(ndata);
	P_PST.resize(ndata);
	s_PST.resize(ndata);
	T_PST.resize(ndata);
	
	unsigned int ie=0;
	fscanf(fptr,"%lf ",&eps);
	while(!feof(fptr)){
		epsilon_PST[ie]=eps;
		fscanf(fptr,"%lf %lf %lf",&P_PST[ie],&s_PST[ie],&T_PST[ie]);
		//printf("eps=%g, P=%g, s=%g,T=%g\n",epsilon_PST[ie],P_PST[ie],s_PST[ie],T_PST[ie]);
		ie+=1;
		fscanf(fptr,"%lf ",&eps);
	}
	
	int nres=reslist->resmap.size();
	vector<double> density;
	vector<double> maxweight;
	Eigen::Matrix3d chi;
	density.resize(nres);
	maxweight.resize(nres);
	fclose(fptr);
	
}

void CEoS::BuildMap(){
	unsigned int ie;
	//printf("check BuildMap, epsilon_PST.size()=%d\n",int(epsilon_PST.size()));
	for(ie=0;ie<epsilon_PST.size();ie++){
		etmap.insert(pairdi(T_PST[ie],ie));
		//printf("ie=%d, etmap.size()=%d\n",ie,int(etmap.size()));
	}
}

void CEoS::GetEoSFromEpsilon_PST(double epsilonset){
	Eigen::Matrix3d chi;
	Eigen::Matrix3d sigma;
	double depsilon=0.02,w0,TT;
	int ie0,ndata;
	epsilon=epsilonset;
	ie0=lrint(floor(-0.5+epsilon/depsilon));
	if(ie0<0)
		ie0=0;
	if(ie0>=int(epsilon_PST.size()-1)){
		ndata=epsilon_PST.size();
		P=P_PST[ndata]+0.33*(epsilon_PST[ndata]-epsilon);
		TT=T_PST[ndata]*pow(epsilon/epsilon_PST[ndata],0.25);
		s=(P+epsilon)/T;
	}
	else{
		w0=(epsilon-epsilon_PST[ie0+1])/depsilon;
		P=w0*P_PST[ie0]+(1.0-w0)*P_PST[ie0+1];
		TT=w0*T_PST[ie0]+(1.0-w0)*T_PST[ie0+1];
		s=w0*s_PST[ie0]+(1.0-w0)*s_PST[ie0+1];
	}
	if(fabs(T-0.140)<0.01){
		printf("TESTING, s=%g\n",s);
		double Ph,eh,nh;
		int nres=reslist->resmap.size();
		vector<double> density,maxweight;;
		density.resize(nres);
		maxweight.resize(nres);
		MSU_EOS::CalcEoSandTransportCoefficients(TT,reslist,eh,Ph,nh,density,chi,sigma);
		Ph=Ph/1000.0;
		eh=eh/1000.0;
	}
}

void CEoS::Print(){
	printf("---- T=%g, P=%g, epsilon=%g, s=%g  ----\n",T,P,epsilon,s);
	printf("chill/s=%g, chiud/s=%g, chils/s=%g, chiss/s=%g\n",chill,chiud,chils,chiss);
}

void CEoS::PrintChi(){
	Eigen::Matrix3d chimat;
	chimat(0,0)=chimat(1,1)=chill;
	chimat(1,0)=chimat(0,1)=chiud;
	chimat(0,2)=chimat(1,2)=chimat(2,1)=chimat(2,0)=chils;
	chimat(2,2)=chiss;
	printf("--------- Chi(Tf)-------------\n");
	cout<< chimat << endl;
	cout << chimat << endl;
}

void CEoS::GetEoSFromT_PST(double Tset){
	mapdi::iterator iter;
	int ie;
	double w0,delT;
	T=Tset;	
	iter=etmap.lower_bound(T);
	ie=iter->second;
	//printf("first=%g, second=%d\n",iter->first,iter->second);
	//printf("ie=%d\n",ie);
	if(ie!=0)
		ie-=1;
	delT=T_PST[ie+1]-T_PST[ie];
	w0=(T_PST[ie+1]-T)/delT;
	epsilon=w0*epsilon_PST[ie]+(1.0-w0)*epsilon_PST[ie+1];
	P=w0*P_PST[ie]+(1.0-w0)*P_PST[ie+1];
	s=w0*s_PST[ie]+(1.0-w0)*s_PST[ie+1];
	
}
