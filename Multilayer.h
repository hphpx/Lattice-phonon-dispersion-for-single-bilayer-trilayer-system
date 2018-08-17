#ifndef  Multilayer_H
#define  Multilayer_H

#include"DynTrilayer.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <list>

using Eigen::MatrixXd;
using namespace Eigen;

struct eigensys
{
	double evr[6];
	double evi[6];
	double val;
};

struct eigenwhole
{
	 double val;
	 VectorXcd ev;
};

class Multilayer: public DynTrilayer
{
	public:
		void gensecondTri();
		void gensecondRho();
		void genfirst();
		void gensecondRec();
		void genthird();
		
		
		complex<double> calMone(double (Multilayer::*f)(double,double,double), double kx, double ky, double dz, vector<position> &vec);
		//complex<double> calMcomplexXY(double kx, double ky, double dz, vector<position> &vec);
		//complex<double> calMcomplexYY(double kx, double ky, double dz, vector<position> &vec);
		
		complex<double> calMmergeXX(double kx, double ky, double dz, vector<position> &vec);
		complex<double> calMmergeXY(double kx, double ky, double dz, vector<position> &vec);
		complex<double> calMmergeYY(double kx, double ky, double dz, vector<position> &vec);
		
		MatrixXcd mat,matfixed,matk,matsub;
	
		void solveM();
		void solveM_pl();
		void solveM_pl_line();
		void filleigens();
	    void checkdiff();
		void initial();
		void outputdata();
		void outputdata_pl();
		void outputdata_pl_line();
		vector<double> eigenresult;
		vector<double> checken;
		vector<vector<double> > allresult;
		vector<position> posq,posqc,pos2,pos2q,pos1,pos3;
	    double yukawaXX(double,double,double);
		double yukawaXY(double,double,double);
		double yukawaYY(double,double,double);
		void showlattice_quater(vector<position>&);
		void showlattice_full(vector<position>&);
		ComplexEigenSolver<MatrixXcd> ces;
		VectorXcd b;
		int dim;
		int lastnum;
		int smart;
		vector<list<eigensys> > allres;
		list<eigensys> res;
	    vector<eigenwhole> res1;
		vector<vector<eigenwhole> > allres1;
		
		string sd,spre;
		
		
		
};



#endif
