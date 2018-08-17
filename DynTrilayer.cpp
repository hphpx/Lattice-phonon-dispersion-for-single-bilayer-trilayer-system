#include"DynTrilayer.h"
#include<iostream>
#include<fstream>
#include<string>
#include<iomanip>
#include<cmath>


void DynTrilayer::setpara(double a, double b)
{
    R=a; theta=b; 	
    this->a=sqrt(3*Pi/sin(2*Pi*theta/180));
}

void DynTrilayer::setparaRec(double a, double b)
{
    R=a; t=b; 	
    this->x=sqrt(2*Pi/t);
}

void DynTrilayer::latticeRec(double R)
{
	double Rx=R/x;
	double temp;
	int m,n,end;
	position tempos;
	for(m=1;;m++)
	{
		temp=Rx*Rx-m*m*t*t;
		if(temp<0) break;
		
		end=sqrt(temp);
		tempos.y=m*t*x;
		for(n=1;n<=end;n++)
		{
			tempos.x=n*x;
			pos.push_back(tempos);
			
			tempos.x=-tempos.x; pos.push_back(tempos);
			tempos.y=-tempos.y; pos.push_back(tempos);
			tempos.x=-tempos.x; pos.push_back(tempos);
			
		}
	}
	
	end=Rx/t; tempos.x=0;
	for(m=1;m<=end;m++)
	 {
	 	tempos.y=m*t*x;  pos.push_back(tempos);
	 	tempos.y=-tempos.y; pos.push_back(tempos);
	 }
	 
	end=Rx; tempos.y=0;
	for(n=1;n<=end;n++)
	 {
	 	tempos.x=n*x;  pos.push_back(tempos);
	 	tempos.x=-tempos.x; pos.push_back(tempos);
	 }
	
	tempos.x=0; tempos.y=0;
	pos.push_back(tempos);
}

void DynTrilayer::lattice(double R)
{
//	R+=2;
	
	double ct=cos(2*Pi*theta/180);
	double st=sin(2*Pi*theta/180);
	double c=cos(Pi*theta/180);
	double s=sin(Pi*theta/180);
//	a=sqrt(3*pi/st);
	ta=R/a;
//	cout<<2*pi*theta/180<<endl;
//	cout<<pi*theta/180<<endl;
    cout<<a<<" "<<R<<endl;
	cout<<ta<<endl;
	cout<<ct<<" "<<st<<endl;
	position tempos;
	int m,n,begin,end;
	double delta,fl,fr,fln,frn;
	for(m=0;;m++)
	 {  cout<<m<<endl;
	 	delta=ta*ta-m*m*st*st;
	 	if(delta<0) break;
	 	delta=sqrt(delta);
	 	fl=-m*ct-delta;
	 	fr=-m*ct+delta;
	 	
	 	cout<<fl<<" -- "<<fr<<endl;
	 //	if(m>fl)  begin=int(m); else begin=int(fl);
	    if(fl<0) begin=int(fl);
	    else begin=int(fl)+1;
	 	
		if(fr>0) end=int(fr);
		else end=int(fr)-1;
	 	
	 	for(n=begin;n<=end;n++)
	 	{
	 		tempos.x=(n-m)*a*s;
	 		tempos.y=(n+m)*a*c;
	 		pos.push_back(tempos);
	 	}
	 	
	 	if(m==0) continue;
	 	
	 	fln=m*ct-delta;
	 	frn=m*ct+delta;
	 	
	 	//if(m>fln)  begin=int(m); else begin=int(fln);
	 	//end=int(frn);
	 	cout<<fln<<" ++ "<<frn<<endl;
	 	if(fln<0) begin=int(fln);
	    else begin=int(fln)+1;
	 	if(frn>0) end=int(frn);
		else end=int(frn)-1;
	 	
	 	for(n=begin;n<=end;n++)
	 	{
	 		tempos.x=(n+m)*a*s;
	 		tempos.y=(n-m)*a*c;
	 		pos.push_back(tempos);
	 	}
	 }
	 
//	 R-=2;
	 
/*	 int i;
	 int cut=R/(2*a*s);
	 for(i=0;i<=cut;i++)
	  xpos.push_back(i*2*a*s);
	 cut=R/(2*a*c);
	 for(i=1;i<=cut;i++)
	  ypos.push_back(i*2*a*c);*/
}

/*void DynTrilayer::filllattice(double d, double shift, vector<position> &vec, vector<double> &vecx, vector<double &vecy)
{
	double x,y,td
	position tempos;
	td=d*d; tR=R*R;
	for(i=0;i<pos.size();i++)
	 {
	 	x=pos[i].x;
	 	y=pos[i].y+shift;
	 	r=x*x+y*y+td;
	 	if(r<=tR) 
	 	 {
	 	   tempos.x=x; tempos.y=y; vec.push_back(pos[i]);
	     }
	 }
	 
	  r=sqrt(tR-td);
	 if(shift==0)
	 {
	  cut=r/(2*a*s);
	  for(i=0;i<=cut;i++)
	   vecx.push_back(i*2*a*s);
     }
     
     else{
     	    if(fabs(shift-a*cos(pi*theta/180))<1e-6)
			     
     	
     }
	  
	
	 cut=r/(2*a*c);
	 for(i=1;i<=cut;i++)
	  vecy.push_back(i*2*a*c+shift);
}	 */

void DynTrilayer::showlattice(vector<position> &pos)
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
	 // output<<-pos[i].x<<" "<<-pos1[i].y<<endl;
	 // output<<-pos[i].x<<" "<<pos1[i].y<<endl;
	 // output<<pos[i].x<<" "<<-pos1[i].y<<endl;
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
