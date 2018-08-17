#include"trilayertri.h"
#include<iostream>

using namespace std;
int main()
{
	trilayertri sol;
	sol.dim=6;
	sol.kmin=0.0001;
	sol.kh=0.1;
	sol.initial();
	sol.d=0.2;
	sol.kappa=0.4;
	sol.setpara(50,30);
	sol.lattice(sol.R+sol.a*2);
	sol.genfirst();
//	sol.showlattice(sol.pos);
    sol.gensecondTri();
    //sol.genthird();
   // sol.showlattice_quater(sol.posq);
    //sol.showlattice_quater(sol.pos2q);
//	sol.showlattice_full(sol.pos3);

   
   sol.sd="0-2";
   sol.spre="sparse1k04";
   int i;
   double degree[5]={0,30,60,90};
   double kmaxii[5]={8,6,8,6};
   for(i=0;i<=0;i++)
   {
   	sol.ktheta=degree[i];
   	sol.kmax=kmaxii[i];
   	sol.subinitial();
   	cout<<sol.allres1.size()<<" ***** "<<endl;
   	sol.solveALL(); //cout<<"dddd\n";
    sol.outputdata_pl_line(); // cout<<"112345\n";
    sol.clearall();
   }
   
   
	return 0;
}
