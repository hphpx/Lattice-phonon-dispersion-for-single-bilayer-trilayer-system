#ifndef  DynTrilayer_H
#define  DynTrilayer_H

#define  Pi  3.14159265359
#include<vector>
using namespace std;

struct position
{
	double x,y;
};

class DynTrilayer
{
	public:
	//	DynTrilayer()
		void lattice(double);
	//	void calM(int dn);
		void showlattice(vector<position> &);
	    //void filllattice(double d, double shift);
	    //void fillaxis(double d, double shift)
	    
	    
		void setpara(double,double);
		double R,a,theta,ta;
		
		double d,kappa,ktheta,kmax,kmin,ukx,uky,kh;
		
		double x,t;
		void latticeRec(double); 
		void setparaRec(double,double);
	//	double M1,M2,M3;
		vector<position> pos;//,pos1,pos2,pos3;
	//	vector<double> xpos,ypos,xpos2,ypos2,xpos3,ypos3;
};

#endif

