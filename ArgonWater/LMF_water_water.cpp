// this file is used to calculate the lmfBA and lmfBB potential between 
// for the symmetric LMF theory. The mimic system have the LJ interaction between
// water truncated

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

double Fx(double RR, double r, double theta);

int main()
{
int nbinr=736; // this is the bin number of rdf
vector<double> r(nbinr,0.0),rho(nbinr,0.0); // rho stores the full rdf_BA
ifstream ifile;
ifile.open("rdf_OA.xvg");
for (int i=0; i<nbinr; i++){
	ifile>>r[i]>>rho[i];
}
ifile.close();

vector<double> rho0(nbinr,0.0); // rho stores the RC rdf_BA
ifile.open("rdf_OA_RC.xvg");
for (int i=0; i<nbinr; i++){
	ifile>>r[i]>>rho0[i];
}
ifile.close();

double factor=0.0;
for(int i=nbinr-50; i<nbinr; i++){
  factor+=rho[i]/50;
}
cout<<"factor: "<<factor<<"\n";
for(int i=0; i<nbinr; i++){
  rho[i]=rho[i]/factor;
}
  
double rhobulk=1000/pow(2.96635,3);
for (int i=0; i<nbinr; i++){
	rho[i]=rho[i]*rhobulk;
}

factor=0.0;
for(int i=nbinr-50; i<nbinr; i++){
  factor+=rho0[i]/50;
}
cout<<"factor: "<<factor<<"\n";
for(int i=0; i<nbinr; i++){
  rho0[i]=rho0[i]/factor;
}
  

for (int i=0; i<nbinr; i++){
	rho0[i]=rho0[i]*rhobulk;
}

double stepr=r[1]-r[0];
double deltheta=0.001*M_PI;
int NTHETA=int(M_PI/deltheta);
vector<double> FrBA(nbinr,0.0);
#pragma omp parallel for simd
for (int i=0; i<nbinr; i++){
	for (int j=0; j<nbinr; j++){
		for (int k=0; k<NTHETA; k++){
			double theta=(k+0.5)*deltheta;
			FrBA[i]+=Fx(r[j],r[i],theta)*2*M_PI*pow(r[j],2)*stepr*(rho[j]-rhobulk)*sin(theta)*deltheta;
		}
	}
}

double inner_cut=0;
for (int i=0; i<nbinr; i++){
  if (r[i]<inner_cut){
    FrBA[i]=0.0;
  }
}

// transform AA to BA to BB
vector<double> FrBB_AA(nbinr,0.0);
#pragma omp parallel for simd
for (int i=0; i<nbinr; i++){
	for (int j=0; j<nbinr; j++){
		for (int k=0; k<NTHETA; k++){
			double theta=(k+0.5)*deltheta;
			double dis=sqrt(pow(r[j],2)+pow(r[i],2)-2*r[j]*r[i]*cos(theta));
			int IBIN=int(dis/stepr);
			if (IBIN<nbinr){
				FrBB_AA[i]+=-FrBA[IBIN]*(r[j]*cos(theta)-r[i])/dis*2*M_PI*pow(r[j],2)*stepr*(rho0[j]-rhobulk)*sin(theta)*deltheta;
			}
		}
	}
}



vector<double> urBB_AA(nbinr,0.0);
for (int i=1; i<(nbinr-1); i++){
	for (int j=i; j<(nbinr-1); j++){
		urBB_AA[i]+=(FrBB_AA[j]+FrBB_AA[j+1])/2*stepr;
	}
}
urBB_AA[0]=urBB_AA[1]+FrBB_AA[1]*stepr;


ofstream ofile3("ulmf_OO.txt");
for (int i=0; i<nbinr; i++){
	ofile3<<r[i]<<"  "<<urBB_AA[i]<<"\n";
}
ofile3.close();
// write out the result






return 0;

}


double Fx(double RR, double r, double theta){
	double eps=0.650194;
	double sigma=0.316557;
	double rwca=sigma*pow(2.0,1.0/6.0);
	double rcut=1.25;
	
	double dis=sqrt(pow(RR,2)+pow(r,2)-2*RR*r*cos(theta));
	
	double F=0.0;
	if (dis<=rwca){
		F=0.0;
	}else if(dis<=rcut){
		double Fr=24*eps/dis*pow(sigma/dis,6)*(2*pow(sigma/dis,6)-1);
		F=-Fr*(RR*cos(theta)-r)/dis;
	}else{
		F=0.0;
	}
	
	return F;
}


