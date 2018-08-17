#include"trilayertri.h"
#include"Multilayer.cpp"
void trilayertri::buildfixedmatrix()
{
	matsub.resize(2,2);
	matfixed=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(0,0,0,posq);
	matsub(1,1)=calMmergeYY(0,0,0,posq);
	matsub(0,1)=0;//calMmergeXY(0,0,0,posq); 
	matsub(1,0)=matsub(0,1);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d,pos2);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d*2,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d*2,pos2);
	
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
}

void trilayertri::buildkmatrix(double kx,double ky)
{
	matk=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(kx,ky,0,posq);
	matsub(1,1)=calMmergeYY(kx,ky,0,posq);
	matsub(0,1)=calMmergeXY(kx,ky,0,posq);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,0,2,2)+=matsub;
	matk.block(2,2,2,2)+=matsub;
	matk.block(4,4,2,2)+=matsub;
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,kx,ky,d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,kx,ky,d,pos2);
	matsub(0,1)=calMone(&Multilayer::yukawaXY,kx,ky,d,pos2);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,2,2,2)+=matsub;
	matk.block(4,0,2,2)+=matsub;
	
    matsub(1,1)=conj(matsub(1,1)); 
    matsub(0,0)=conj(matsub(0,0));
    matsub(0,1)=conj(matsub(0,1));
    matsub(1,0)=conj(matsub(1,0));
	
	matk.block(0,4,2,2)+=matsub;
	matk.block(2,0,2,2)+=matsub;
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,kx,ky,2*d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,kx,ky,2*d,pos2);
	matsub(0,1)=calMone(&Multilayer::yukawaXY,kx,ky,2*d,pos2);
	matsub(1,0)=matsub(0,1);
	
	matk.block(2,4,2,2)+=matsub;
    
	matsub(1,1)=conj(matsub(1,1)); 
    matsub(0,0)=conj(matsub(0,0));
    matsub(0,1)=conj(matsub(0,1));
    matsub(1,0)=conj(matsub(1,0));
	
	matk.block(4,2,2,2)+=matsub;
	
	
}

//////////
void trilayertri::buildfixedmatrix_no_shift()
{
	matsub.resize(2,2);
	matfixed=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(0,0,0,posq);
	matsub(1,1)=calMmergeYY(0,0,0,posq);
	matsub(0,1)=0;//calMmergeXY(0,0,0,posq); 
	matsub(1,0)=matsub(0,1);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
	//matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d,pos2);
	//matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d,pos2);
	
	matsub(0,0)=calMmergeXX(0,0,d,pos2q);
	matsub(1,1)=calMmergeYY(0,0,d,pos2q);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
//	matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d*2,pos2);
//	matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d*2,pos2);

    matsub(0,0)=calMmergeXX(0,0,d*2,posqc);
	matsub(1,1)=calMmergeYY(0,0,d*2,posqc);
	
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
}

void trilayertri::buildkmatrix_no_shift(double kx,double ky)
{
	matk=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(kx,ky,0,posq);
	matsub(1,1)=calMmergeYY(kx,ky,0,posq);
	matsub(0,1)=calMmergeXY(kx,ky,0,posq);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,0,2,2)+=matsub;
	matk.block(2,2,2,2)+=matsub;
	matk.block(4,4,2,2)+=matsub;
	
	matsub(0,0)=calMmergeXX(kx,ky,d,pos2q);
	matsub(1,1)=calMmergeYY(kx,ky,d,pos2q);
	matsub(0,1)=calMmergeXY(kx,ky,d,pos2q);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,2,2,2)+=matsub;
	matk.block(4,0,2,2)+=matsub;
	
   // matsub(1,1)=conj(matsub(1,1)); 
   // matsub(0,0)=conj(matsub(0,0));
   // matsub(0,1)=conj(matsub(0,1));
   // matsub(1,0)=conj(matsub(1,0));
	
	matk.block(0,4,2,2)+=matsub;
	matk.block(2,0,2,2)+=matsub;
	
//	matsub(0,0)=calMone(&Multilayer::yukawaXX,kx,ky,2*d,pos2);
//	matsub(1,1)=calMone(&Multilayer::yukawaYY,kx,ky,2*d,pos2);
//	matsub(0,1)=calMone(&Multilayer::yukawaXY,kx,ky,2*d,pos2);
   	matsub(0,0)=calMmergeXX(kx,ky,d*2,posqc);
	matsub(1,1)=calMmergeYY(kx,ky,d*2,posqc);
	matsub(0,1)=calMmergeXY(kx,ky,d*2,posqc);
	matsub(1,0)=matsub(0,1);
	
	matk.block(2,4,2,2)+=matsub;
    
//	matsub(1,1)=conj(matsub(1,1)); 
   // matsub(0,0)=conj(matsub(0,0));
 //   matsub(0,1)=conj(matsub(0,1));
   // matsub(1,0)=conj(matsub(1,0));
	
	matk.block(4,2,2,2)+=matsub;
	
	
}

/////////////////////////////////////////////////////////////////////////

void trilayertri::buildfixedmatrixSP()
{
	matsub.resize(2,2);
	matfixed=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(0,0,0,posq);
	matsub(1,1)=calMmergeYY(0,0,0,posq);
	matsub(0,1)=0;//calMmergeXY(0,0,0,posq); 
	matsub(1,0)=matsub(0,1);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d,pos2);
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(2,2,2,2)+=matsub;
	
	matfixed.block(0,0,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	//matsub(0,0)=calMone(&Multilayer::yukawaXX,0,0,d*2,pos2);
	//matsub(1,1)=calMone(&Multilayer::yukawaYY,0,0,d*2,pos2);
	
	matsub(0,0)=calMmergeXX(0,0,d*2,posqc);
	matsub(1,1)=calMmergeYY(0,0,d*2,posqc);
	
	matfixed.block(2,2,2,2)+=matsub;
	matfixed.block(4,4,2,2)+=matsub;
	
	
}

void trilayertri::buildkmatrixSP(double kx,double ky)
{
	matk=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(kx,ky,0,posq);
	matsub(1,1)=calMmergeYY(kx,ky,0,posq);
	matsub(0,1)=calMmergeXY(kx,ky,0,posq);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,0,2,2)+=matsub;
	matk.block(2,2,2,2)+=matsub;
	matk.block(4,4,2,2)+=matsub;
	
	matsub(0,0)=calMone(&Multilayer::yukawaXX,kx,ky,d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,kx,ky,d,pos2);
	matsub(0,1)=calMone(&Multilayer::yukawaXY,kx,ky,d,pos2);
	matsub(1,0)=matsub(0,1);
	
	matk.block(0,2,2,2)+=matsub;
	matk.block(0,4,2,2)+=matsub;

	
    matsub(1,1)=conj(matsub(1,1)); 
    matsub(0,0)=conj(matsub(0,0));
    matsub(0,1)=conj(matsub(0,1));
    matsub(1,0)=conj(matsub(1,0));
	
	matk.block(4,0,2,2)+=matsub;
	matk.block(2,0,2,2)+=matsub;
	
/*	matsub(0,0)=calMone(&Multilayer::yukawaXX,kx,ky,2*d,pos2);
	matsub(1,1)=calMone(&Multilayer::yukawaYY,kx,ky,2*d,pos2);
	matsub(0,1)=calMone(&Multilayer::yukawaXY,kx,ky,2*d,pos2);
	matsub(1,0)=matsub(0,1);*/
	
	matsub(0,0)=calMmergeXX(kx,ky,2*d,posqc);
	matsub(1,1)=calMmergeYY(kx,ky,2*d,posqc);
	matsub(0,1)=calMmergeXY(kx,ky,2*d,posqc);
	matsub(1,0)=matsub(0,1);
	
	matk.block(2,4,2,2)+=matsub;	
	matk.block(4,2,2,2)+=matsub;
	
	
}

//////////////////////
void trilayertri::buildfixedmatrix_single_tri()
{
	matsub.resize(2,2);
	matfixed=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(0,0,0,posq);
	matsub(1,1)=calMmergeYY(0,0,0,posq);
	matsub(0,1)=0;//calMmergeXY(0,0,0,posq); 
	matsub(1,0)=matsub(0,1);
	
	matfixed=matsub;

	
	
}

void trilayertri::buildkmatrix_single_tri(double kx,double ky)
{
	matk=MatrixXcd::Zero(dim,dim);
	
	matsub(0,0)=calMmergeXX(kx,ky,0,posq);
	matsub(1,1)=calMmergeYY(kx,ky,0,posq);
	matsub(0,1)=calMmergeXY(kx,ky,0,posq);
	matsub(1,0)=matsub(0,1);
	
	matk=matsub;
	
	
}
///////////////////////
void trilayertri::solveALL()
{
	double len; 
	lastnum=0;
	int g;
	buildfixedmatrix(); // cout<<"havvv\n";
	for(len=kmin;len<kmax;len+=kh)
	{
	//	fillmatrixbilayer(len*ukx,len*uky);
	 // cout<<"888888888\n";
	
	      buildkmatrix(len*ukx,len*uky);
	      mat=matfixed-matk;
	     // cout<<matk<<endl;
	    // cout<<matfixed<<endl;
		//	fillmatrixbilayerRec(len*ukx,len*uky);  
		//	 cout<<len<<endl; //cin>>i;
			// cout<<mat<<endl; cout<<ktheta<<"  "; cout<<"------------------------\n";
			// cin>>g;
		solveM_pl_line();  
		lastnum++; 
	}
}

void trilayertri::smart_solveALL()
{
	double len=kmin; 
	lastnum=0;
	int g;
	buildfixedmatrix(); // cout<<"havvv\n";
	for(;len<kmax;)
	{
	//	fillmatrixbilayer(len*ukx,len*uky);
	 // cout<<"888888888\n";
	 if(lastnum!=0)
	 {
	 
	   checkdiff();
	   if(smart==1) kh=0.0002;
	   else kh=0.01;
	 
       len+=kh;
     }
	 
	      buildkmatrix(len*ukx,len*uky);
	      mat=matfixed-matk;
	     // cout<<matk<<endl;
	    // cout<<matfixed<<endl;
		//	fillmatrixbilayerRec(len*ukx,len*uky);  
		//	 cout<<len<<endl; //cin>>i;
			// cout<<mat<<endl; cout<<ktheta<<"  "; cout<<"------------------------\n";
			// cin>>g;
		solveM_pl_line();  
		lastnum++; 
		cout<<lastnum<<" "<<len<<"  "<<smart<<endl; 
		cin>>g;
	}
}

void trilayertri::solveALL_single_tri()
{
	double len; 
	lastnum=0;
	int g;
	buildfixedmatrix_single_tri();
	for(len=kmin;len<kmax;len+=kh)
	{
	//	fillmatrixbilayer(len*ukx,len*uky);
	
	      buildkmatrix_single_tri(len*ukx,len*uky);
	      mat=matfixed-matk;
		//	fillmatrixbilayerRec(len*ukx,len*uky);  
		//	 cout<<len<<endl; //cin>>i;
	//	cout<<matfixed<<endl;
	//	cout<<matk<<endl;
	//		 cout<<mat<<endl; 
			 //cout<<ktheta<<"  ";
		//	  cout<<"------------------------\n";
		//	 cin>>g;
		solveM_pl_line();  
		lastnum++; 
	}
}

void trilayertri::solveALL_no_shift()
{
	double len; 
	lastnum=0;
	int g;
	buildfixedmatrix_no_shift();
	for(len=kmin;len<kmax;len+=kh)
	{
	//	fillmatrixbilayer(len*ukx,len*uky);
	
	      buildkmatrix_no_shift(len*ukx,len*uky);
	      mat=matfixed-matk;
		//	fillmatrixbilayerRec(len*ukx,len*uky);  
		//	 cout<<len<<endl; //cin>>i;
			// cout<<mat<<endl; cout<<ktheta<<"  "; cout<<"------------------------\n";
		//	 cin>>g;
		solveM_pl_line();  
		lastnum++; 
	}
}

void trilayertri::solveALL_SP()
{
	double len; 
	lastnum=0;
	int g;
	buildfixedmatrixSP();
	for(len=kmin;len<kmax;len+=kh)
	{
	//	fillmatrixbilayer(len*ukx,len*uky);
	
	      buildkmatrixSP(len*ukx,len*uky);
	      mat=matfixed-matk;
		//	fillmatrixbilayerRec(len*ukx,len*uky);  
		//	 cout<<len<<endl; //cin>>i;
			// cout<<mat<<endl; cout<<ktheta<<"  "; cout<<"------------------------\n";
		//	 cin>>g;
		solveM_pl_line();  
		lastnum++; 
	}
}



void Multilayer::initial()
{
	ktheta=0;
	//kmax=4;
	//kmin=0.03;
//	kmin=0.0001;
//	d=0;
//	kh=0.0002;
//	kappa=0.4;
	//dim=2;
	ukx=cos(ktheta*Pi/180);
	uky=sin(ktheta*Pi/180);
	cout<<ukx<<" u "<<uky<<endl;
	
	mat.resize(dim,dim);
	
	
	
}

void trilayertri::clearall()
{
	allres1.clear();
}

void trilayertri::subinitial()
{
	
	ukx=cos(ktheta*Pi/180);
	uky=sin(ktheta*Pi/180);
	cout<<ukx<<" u "<<uky<<endl;
	
}
