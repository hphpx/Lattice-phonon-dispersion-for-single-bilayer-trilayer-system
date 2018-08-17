#include"Multilayer.h"
#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>


bool compare_eigen(eigensys &a, eigensys &b)
{
	return a.val<b.val;
}

bool compare_eigen1(eigenwhole a, eigenwhole b)
{
	return a.val<b.val;
}

void Multilayer::genfirst()
{  
   int i;
   double x,y,tR=R*R;
	for(i=0;i<pos.size();i++)
	  {
	  	x=pos[i].x; y=pos[i].y;
	  	if(x!=0 || y!=0) pos1.push_back(pos[i]);
	  	
	  	if(x*x+y*y<=tR)
	  	  if(x+y>0 && x*y>=0)
		   {
		    posq.push_back(pos[i]); posqc.push_back(pos[i]);
		   }
		 // else { cout<<x*y<<" - "<<x+y<<endl;
		  //}
	  }
	  
	  position temp;
	  temp.x=0; temp.y=0;
	  posqc.push_back(temp);
}




void Multilayer::genthird()
{  
   int i;
   double x,y,z,temp,shift,tR=R*R;
   z=4*d*d;
   position p;
   shift=a/sqrt(3);
   temp=tR-z;
	for(i=0;i<pos.size();i++)
	  {
	  	x=pos[i].x; y=pos[i].y-shift;
	  	if(x*x+y*y<=temp)
	  	 // if(x+y>0 && x*y>=0)
		   {  p.x=x; p.y=y;
		     pos3.push_back(p);
		   }
	  }
}

void Multilayer::gensecondRec()
{  
   int i;
   double y,z,shiftx,shifty,temp,tR=R*R;
   position p;
   z=d*d;
   shiftx=x/2;  shifty=t*x/2; //cout<<shift<<" 7777777"<<endl;
   cout<<shiftx<<" -8- "<<shifty<<" "<<t<<endl;
   double x;
   temp=tR-z;
	for(i=0;i<pos.size();i++)
	  {
	  	x=pos[i].x+shiftx; y=pos[i].y+shifty;
	  	//if(pos[i].x==0 && pos[i].y==0) cout<<x<<" *** "<<y<<endl;
	  	if(x>0 && y>0 && x*x+y*y<=temp)
	  	  {
	  	  	p.x=x; p.y=y; pos2q.push_back(p);  //cout<<x<<" * "<<y<<endl;
	  	  }
	  }
}


void Multilayer::gensecondTri()
{  
   int i;
   double x,y,z,shift,temp,tR=R*R;
   position p;
   z=d*d;
   
   shift=a/sqrt(3);
    cout<<shift<<" 7777777"<<endl;
  
   temp=tR-z;
	for(i=0;i<pos.size();i++)
	  {
	  	x=pos[i].x; y=pos[i].y+shift;
	  	if(x*x+y*y<=temp)
	  	  {
	  	  	p.x=x; p.y=y; pos2.push_back(p);
	  	  }
	  }
}

void Multilayer::gensecondRho()
{  
   int i;
   double x,y,z,shift,temp,tR=R*R;
   position p;
   z=d*d;
  
   shift=a*cos(theta*Pi/180); cout<<shift<<" 7777777"<<endl;
  
   temp=tR-z;
	for(i=0;i<pos.size();i++)
	  {
	  	x=pos[i].x; y=pos[i].y+shift;
	  	if(x+y>0 && x*y>=0 && x*x+y*y<=temp)
	  	  {
	  	  	p.x=x; p.y=y; pos2q.push_back(p);
	  	  }
	  }
}

complex<double> Multilayer::calMone(double (Multilayer::*f)(double,double,double), double kx, double ky, double dz, vector<position> &vec)
{
   	double real, img, temp,temp1,sumr=0,sumi=0;
   	int i; 
	if(kx==0 && ky==0)
	 {
	 	for(i=0;i<vec.size();i++)
	  {
	
		temp1=(*this.*f)(vec[i].x, vec[i].y, dz);
		sumr+=temp1;
		
	  }
	  sumi=0;
	 }
	else 
	{
	 
	for(i=0;i<vec.size();i++)
      {
	
		temp=vec[i].x*kx+vec[i].y*ky;
		real=cos(temp);
		img=sin(temp);
		temp1=(*this.*f)(vec[i].x, vec[i].y, dz);
		sumr+=(real*temp1);
		sumi+=(img*temp1);
	  }
    }
//	temp1*=0.5;
	complex<double> result(sumr,sumi);
	return result;
}


complex<double> Multilayer::calMmergeXX(double kx, double ky, double dz, vector<position> &vec)
{
	double temp,sum=0;
	int i;
	if(kx==0 && ky==0)
	{
		for(i=0;i<vec.size();i++)
	   {
		if(vec[i].x==0) 
		  temp=2;
		else if(vec[i].y==0)
		  temp=2;
		else temp=4;
		
		temp*=yukawaXX(vec[i].x,vec[i].y,dz);
		
	
		sum+=temp;
	   }
	}
	else
	{
	
	for(i=0;i<vec.size();i++)
	 {
		if(vec[i].x==0) 
		  temp=cos(vec[i].y*ky)*2;
		else if(vec[i].y==0)
		  temp=cos(vec[i].x*kx)*2;
		else temp=cos(vec[i].x*kx)*cos(vec[i].y*ky)*4;
		
		temp*=yukawaXX(vec[i].x,vec[i].y,dz);
		
	
		sum+=temp;
	 }
    }
	
	//	temp*=0.5;
		
		complex<double> result(sum,0);
		return result;
}

complex<double> Multilayer::calMmergeYY(double kx, double ky, double dz, vector<position> &vec)
{
	double temp,sum=0;
	int i;
	if(kx==0 && ky==0)
	{
		for(i=0;i<vec.size();i++)
	   {
		if(vec[i].x==0) 
		  temp=2;
		else if(vec[i].y==0)
		  temp=2;
		else temp=4;
		
		temp*=yukawaXX(vec[i].y,vec[i].x,dz);
		
	
		sum+=temp;
	   }
	}
	else
	{
	
	for(i=0;i<vec.size();i++)
	 {
		if(vec[i].x==0) 
		  temp=cos(vec[i].y*ky)*2;
		else if(vec[i].y==0)
		  temp=cos(vec[i].x*kx)*2;
		else temp=cos(vec[i].x*kx)*cos(vec[i].y*ky)*4;
		
		temp*=yukawaXX(vec[i].y,vec[i].x,dz);
		sum+=temp;
	
		
	 }
    }
	
	//	temp*=0.5;
		
		complex<double> result(sum,0);
		return result;
}




complex<double> Multilayer::calMmergeXY(double kx, double ky, double dz, vector<position> &vec)
{
	double temp,sum=0;
	int i;
	if(kx*ky!=0)
	{
	
	for(i=0;i<vec.size();i++)
	{
	    if(vec[i].x==0 || vec[i].y==0) continue;
	    
	    temp=sin(vec[i].x*kx)*sin(vec[i].y*ky)*4;
		
		temp*=yukawaXY(vec[i].x,vec[i].y,dz);
		
		sum+=temp;
		
	}
    }
    
	
	//	sum*=0.5
		
		complex<double> result(sum,0);
		return result;
}


double Multilayer::yukawaXX(double x, double y, double dz)
{
	double sum=0, temp,r,xs,ys,ks;
	
	xs=x*x; ys=y*y; ks=kappa*kappa;
	r=sqrt(xs+ys+dz*dz);
	sum=ks*xs*xs;
	temp=dz*dz*(ks*xs-kappa*r-1);
	sum+=temp;
	temp=ys*(1+kappa*r);
	sum-=temp;
	temp=xs*(2+ks*ys+2*kappa*r);
	sum+=temp;
	
	sum*=exp(-kappa*r);
	sum/=pow(r,5);
	
	return sum;
}

double Multilayer::yukawaYY(double x, double y, double dz)
{
	return yukawaXX(y,x,dz);
}

double Multilayer::yukawaXY(double x, double y, double dz)
{
	if(x==0 || y==0) return 0;
	double r,sum=0;
	r=sqrt(x*x+y*y+dz*dz);
	
	sum=x*y*(3+3*kappa*r+kappa*kappa*r*r);
	sum*=exp(-kappa*r);
	sum/=pow(r,5);
	
	return sum;
}


void Multilayer::solveM()
{
//	ComplexEigenSolver<MatrixXcd> ces;
    int i;
	ces.compute(mat);
	b.resize(dim);
	 
	b=ces.eigenvalues();
	for(i=0;i<dim;i++)
	 {
	 	eigenresult.push_back(b(i).real());    // cout<<b(i).real()<<endl; 
	 }
//	 cout<<"--------------------------------\n";
	 
	sort(eigenresult.begin(),eigenresult.end());
	
	
	
	allresult.push_back(eigenresult);
	eigenresult.clear();
	
}

void Multilayer::solveM_pl()
{
//	ComplexEigenSolver<MatrixXcd> ces;
    int i,j; 
     //cout<<"* "<<mat<<endl;
	ces.compute(mat); 
	b.resize(dim);
	eigensys temp;
	b=ces.eigenvalues();
	for(i=0;i<dim;i++)
	 {
	 	eigenresult.push_back(b(i).real());    // cout<<b(i).real()<<endl; 
	 	temp.val=b(i).real();
	 	for(j=0;j<dim;j++)
	 	  {
	 	    temp.evr[j]=ces.eigenvectors().col(i)(j).real();
	 	    temp.evi[j]=ces.eigenvectors().col(i)(j).imag();
	 	  }
	 	res.push_back(temp);
	 }
//	 cout<<"--------------------------------\n";
	res.sort(compare_eigen);
	sort(eigenresult.begin(),eigenresult.end());
	
	
	allres.push_back(res);
	res.clear();
	allresult.push_back(eigenresult);
	eigenresult.clear();
	
}

void Multilayer::solveM_pl_line()
{
//	ComplexEigenSolver<MatrixXcd> ces;
    int i,j; 
    // cout<<"* "<<mat<<endl;
	ces.compute(mat); 
	b.resize(dim);
	eigenwhole temp;
//	b=ces.eigenvalues();
//	cout<<dim<<"  -- "<<endl;
	if(lastnum==0)
	{
	//	cout<<mat<<endl; cout<<ktheta<<"  ---------------\n";
	for(i=0;i<dim;i++)
	 {
	 	//eigenresult.push_back(b(i).real());    // cout<<b(i).real()<<endl; 
	 	temp.val=ces.eigenvalues()(i).real();
	  //  cout<<ces.eigenvalues()(i)<<endl;
	 	temp.ev=ces.eigenvectors().col(i);
			// temp.ev[j]=ces.eigenvectors().col(i)(j).real();
	 	    //temp.evi[j]=ces.eigenvectors().col(i)(j).imag();
	 	  
	 	res1.push_back(temp);
	 }
	 //cout<<"--------------------------------\n";
//	res1.sort(compare_eigen1);
    sort(res1.begin(),res1.end(),compare_eigen1);
	//sort(eigenresult.begin(),eigenresult.end());
	
	
//	allres1.push_back(res1);
//	res1.clear();
//	allresult.push_back(eigenresult);
	//eigenresult.clear();
    }
    
    else {
    	
    	 filleigens();
    }
   // cout<<"uuiuiu\n";
    allres1.push_back(res1);
	res1.clear();
}

void Multilayer::filleigens()
{
	int used[6]={0};
	double v1,diff,diff_min;
	eigenwhole temp;
	int minp;
	int i,j;
	for(i=0;i<dim;i++)
	 {
	    b=allres1[lastnum-1][i].ev;
	    diff_min=9999; minp=0;
	    for(j=0;j<dim;j++)
	     {  
	     	if(!used[j])
	     	 {
	     	 	v1=norm(b.dot(ces.eigenvectors().col(j)));
	     	 	diff=fabs(v1-1);
	     	 	if(diff<diff_min)
	     	 	{
	     	 		diff_min=diff; minp=j;
	     	 	}
	     	 }
	     	 
	     	 
	     }
	     
	    temp.val=ces.eigenvalues()(minp).real();
	    temp.ev=ces.eigenvectors().col(minp);
	    res1.push_back(temp);
	    
	    
	    used[minp]=1;
	 }
}

void Multilayer::checkdiff()
{   
    int i;
    checken.clear();
   	for(i=0;i<dim;i++)
   	 {
   	 	checken.push_back(allres1[lastnum-1][i].val);
   	 	
   	 }
   	 
   	 sort(checken.begin(),checken.end());
   	 
   	 double temp, mini=999999.0;
   	 for(i=1;i<dim;i++)
   	  {
   	  	if(checken[i]<0 || checken[i-1]<0) continue;
   	  
   	  	
		temp=sqrt(checken[i]/2)-sqrt(checken[i-1]/2);
		cout<<sqrt(checken[i]/2)<<endl;
   	  	if(mini>temp) mini=temp;
   	  }
   	  cout<<mini<<" ** \n";
   	  if(mini<0.01) smart=1;
   	  else smart=0;
   	  
   	 // cout<<"here\n";
}

void Multilayer::outputdata()
{
	ofstream output;
	cout<<"input data file name:";
	string s;
	cin>>s;
	output.open(s.c_str());
	int i,j;
	double temp;
	for(i=0;i<allresult.size();i++)
	{
	  output<<kmin+kh*(1+i)<<" ";
	  for(j=0;j<allresult[i].size();j++)
	  // output<<allresult[i][j]<<" ";
	  {
	  	if(allresult[i][j]>0) temp=sqrt(allresult[i][j]/2);
	  	else temp=-sqrt(-allresult[i][j]/2);
	  	output<<temp<<" ";
	  } 
	  output<<endl;
	}
	
	output.close();
}

void Multilayer::outputdata_pl()
{
	ofstream output,output1;
	cout<<"input data file name:";
	string s,s1;
	cin>>s;
	s1=s+"-plor.txt";
	s+=".txt";
	output.open(s.c_str());
	output1.open(s1.c_str());
	int i,j;
	double temp;
	for(i=0;i<allresult.size();i++)
	{
	  output<<kmin+kh*(1+i)<<" ";
	  for(j=0;j<allresult[i].size();j++)
	  // output<<allresult[i][j]<<" ";
	  {
	  	if(allresult[i][j]>0) temp=sqrt(allresult[i][j]/2);
	  	else temp=-sqrt(-allresult[i][j]/2);
	  	output<<temp<<" ";
	  } 
	  output<<endl;
	}
	
	output.close();
	list<eigensys>::iterator it;
	output1.setf(ios::fixed);
	for(i=0;i<allres.size();i++)
	 {
	 	output1<<setw(10)<<setprecision(5)<<kmin+kh*(1+i)<<endl;
	 	for(it=allres[i].begin();it!=allres[i].end();it++)
	 	{
	 		output1<<setw(10)<<setprecision(5)<<(*it).val;
	 		for(j=0;j<dim;j++)
	 		 {
	 		 	if(fabs((*it).evr[j])>1e-10) output1<<setw(10)<<setprecision(5)<<(*it).evr[j];
	 		 	else output1<<setw(10)<<setprecision(5)<<0;
	 		 	output1<<",";
			//	output1<<setw(11)<<setprecision(5)<<(*it).evr[j]<<",";
	 		 	output1.setf(ios::left);
	 		 	if(fabs((*it).evi[j])>1e-10) output1<<setw(10)<<setprecision(5)<<(*it).evi[j];
	 		 	else output1<<setw(10)<<setprecision(5)<<0;
	 		 	output1.unsetf(ios::left);
	 		 }
	 		 
	 		 output1<<endl;
	 	}
	 }
	 
	 output1.close();
}

string to_string(int a)
{  
   if(a==0) return "0";
   if(a<0) return "neg";
   string s;
   
   int temp;
	while(a>0)
	{
		temp=a%10;
		s+=(temp+'0');
		a/=10;
	}
	
	reverse(s.begin(),s.end());
	return s;
}


void Multilayer::outputdata_pl_line()
{
	ofstream output,output1;
	cout<<"input data file name:";
	string sp=""; //"C:\\Users\\Hong\\Desktop\\gssg\\deskGS\\qusi2d\\triplelayer\\NewEIGN\\trilayertri\\JINR50\\";
	string s,s1;
//	cin>>s;
    s=spre+"d"+sd+"t"+to_string(ktheta);
	sp+=s;
	s1=sp+"-plor.txt";
	sp+=".txt";
	output.open(sp.c_str());
	output1.open(s1.c_str());
	int i,j,k;
	double temp;
	
	for(i=0;i<allres1.size();i++)
	 {
	 	output<<kmin+kh*i<<" ";
	 //	cout<<kmin+kh*i<<" ";
	 	for(j=0;j<allres1[0].size();j++)
	 	 {
	 	 	if(allres1[i][j].val>0) temp=sqrt(allres1[i][j].val/2);
	    	else temp=-sqrt(-allres1[i][j].val/2);
	  	    output<<temp<<" ";
	 	 }
	 	 output<<endl;
	 }
	
		output.close();
		cout<<"finished1\n";

	for(i=0;i<allres1.size();i++)
	 {
	 	output1<<kmin+kh*i<<endl;
	 	for(j=0;j<allres1[0].size();j++)
	 	 {
	 	 	output1<<allres1[i][j].val<<" ";
	 	 	for(k=0;k<dim;k++)
	 	 	 {
	 	 	 	output1<<allres1[i][j].ev(k).real()<<" "<<allres1[i][j].ev(k).imag()<<" ";
	 	 	 }
	 	 	output1<<endl;
	 	 }
	 }
	


	

	 
	 output1.close();
}



void Multilayer::showlattice_quater(vector<position>  &pos)
{
	ofstream output;
	int i;
	string s;
	cout<<"input file name:";
	cin>>s;
	s+=".txt";
	output.open(s.c_str()); cout<<"hello\n";
	//output.setf(ios::scientific);
	for(i=0;i<pos.size();i++)
	{
	 //cout<<pos[i].x<<" ^^ "<<pos[i].y<<endl;
	 if(pos[i].x==0) 
	   {
	     output<<0<<" "<<pos[i].y<<endl; output<<0<<" "<<-pos[i].y<<endl;
	   }
	 else if(pos[i].y==0) 
	   {
	     output<<pos[i].x<<" "<<0<<endl; output<<-pos[i].x<<" "<<0<<endl;
	   }
	 else{
	 //	cout<<pos[i].x<<" ^--^ "<<pos[i].y<<endl;
	  output<<pos[i].x<<" "<<pos[i].y<<endl;
	  output<<-pos[i].x<<" "<<-pos[i].y<<endl;
	  output<<-pos[i].x<<" "<<pos[i].y<<endl;
	  output<<pos[i].x<<" "<<-pos[i].y<<endl;
         }
     }
/*	for(i=0;i<xpos.size();i++)
	 {
	   output<<xpos[i]<<" "<<0<<endl;
	   output<<-xpos[i]<<" "<<0<<endl;
     }
	for(i=0;i<ypos.size();i++)
	 {
	   output<<0<<" "<<ypos[i]<<endl;
	   output<<0<<" "<<-ypos[i]<<endl;
     }*/

	output.close();

}


void Multilayer::showlattice_full(vector<position>  &pos)
{
	ofstream output;
	int i;
	string s;
	cout<<"input file name:";
	cin>>s;
	s+=".txt";
	output.open(s.c_str()); cout<<"hello\n";
	//output.setf(ios::scientific);
	for(i=0;i<pos.size();i++)
	{
	 
	
	  output<<pos[i].x<<" "<<pos[i].y<<endl;
    }


	output.close();

}
