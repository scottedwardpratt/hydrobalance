#ifndef __CHARGE_H__
#define __CHARGE_H__

#include "msu_commonutils/commondefs.h"
#include "msu_sampler/hyper.h"
#include "msu_hydrobalance/hydro2uds.h"

using namespace std;
using namespace NMSUPratt;

namespace NMSUPratt{

	class CTrajInfo{
	public:
		FILE *fptr;
		int balanceID;
		vector<double> x,y,eta,tau;
		CTrajInfo(int IDset);
		void add(double x1,double y1,double eta1,double tau1);
	};

	class CHBCharge{
	public:
		~CHBCharge(){};
		bool active;
		int q[3];
		double x,y,eta,tau,weight,vx,vy,rapidity;
		Chyper hyper;
		void Propagate(double newtau);
		void SetV(double ux,double uy);
		void Print();
		static CHydroBalance *hb;
		static CB3D *b3d;
		CTrajInfo *trajinfo;
		CHBCharge(){
			trajinfo=NULL;
		};
		void addtraj(){
			if(trajinfo!=NULL){
				trajinfo->add(x,y,eta,tau);
			}
		}
	};

}

#endif