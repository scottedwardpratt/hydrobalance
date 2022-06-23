#include "msu_hydro2uds/charge.h"
#include "msu_commonutils/misc.h"
#include "msu_sampler/hyper.h"
#include "msu_hydro2uds/hydro2uds.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"

using namespace std;

CTrajInfo::CTrajInfo(int IDset){
	balanceID=IDset;
	string filename="trajectories/traj"+to_string(balanceID)+".txt";
	fptr=fopen(filename.c_str(),"w");
}

void CTrajInfo::add(double x1,double y1,double eta1,double tau1){
	x.push_back(x1);
	y.push_back(y1);
	eta.push_back(eta1);
	tau.push_back(tau1);
}

