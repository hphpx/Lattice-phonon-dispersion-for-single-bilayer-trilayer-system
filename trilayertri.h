#ifndef trilayertri_H
#define trilayertri_H

#include"Multilayer.h"

class trilayertri: public Multilayer
{
  public:
	void buildfixedmatrix();
	void buildfixedmatrix_single_tri();
	void buildfixedmatrixSP();
	void buildkmatrix(double,double);
	void buildkmatrix_single_tri(double,double);
	void buildkmatrixSP(double,double);
	void buildfixedmatrix_no_shift();
	void buildkmatrix_no_shift(double,double);

	void solveALL();
	void smart_solveALL();
	void solveALL_SP();
	void solveALL_single_tri();
	void solveALL_no_shift();
	void clearall();
	void subinitial();
	
};

#endif
