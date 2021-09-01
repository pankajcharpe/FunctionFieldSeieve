/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   polynomial_selection.cpp
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:40 AM
 */

#include <cstdlib>


#include <iostream>
#include<string>
#include<NTL/tools.h>
#include<NTL/vector.h>
#include<NTL/matrix.h>
#include<NTL/ZZX.h>
#include<NTL/mat_ZZ.h>
#include<NTL/GF2X.h>
#include <NTL/pair_GF2X_long.h>
#include<NTL/GF2XFactoring.h>
#include<NTL/GF2E.h>
#include<NTL/mat_GF2E.h>
#include<math.h>
#include"ffspol.h"
#include"ffstools.h"
#include<string>


//Resultant by sylvester's matrix method
void resultant(GF2X& r,const ffs_poly& f,const ffs_poly& g)
{
  GF2X p=BuildSparseIrred_GF2X(1000); //set some large bound of degree for sparse irreducible polynomial which will not affect the computation  
  GF2E::init(p);

  Mat<GF2E> M;	
 long d=f.deg+g.deg;
 
 M.SetDims(d,d);
 long n=f.deg;
 long m=g.deg; 
 int i,j;

//Initialize the matrix with zero's
 for(i=0;i<d;i++)
 for(j=0;j<d;j++)
  M[i][j]=conv<GF2E>(0);

for(i=0;i<m;i++)
 for(j=0;j<d;j++)
{
   if(i>=0 && m>i)
  {
   M[i][j]=conv<GF2E>(f.coeffs[n+i-j]);
   
  } 
} 

 for(;i<d;i++)
 for(j=0;j<d;j++)
{
   if((i-j)>=0 && (i-j)<2)
   {
     M[i][j]=conv<GF2E>(g.coeffs[i-j]);
   }
}
r=conv<GF2X>(determinant(M)); 
}

ffs_poly polynomial_selection(const ffs_poly& f,const long n)
{
   ffs_poly g;
   long t=2;
   Vec< Pair< GF2X,long > > factors;
   // vec_GF2X factors;
   factors.SetLength(100);
   long d=ceil(n/f.deg);
   g.deg=1;
  g.coeffs.SetLength(t);
 
   GF2X g0,g1;
   GF2X res;
    	
  while(1)
  {
     g0=random_GF2X(d+1);
     g1=random_GF2X(d+1);

     cout<<"g0="<<g0<<endl;
     cout<<"g0="<<g0<<endl;

     g.coeffs[0]=g0;
     g.coeffs[1]=g1;          
     resultant(res,f,g);  

    cout<<"\n\nResulant="<<res<<endl;

    CanZass(factors,res,(long)0); //defined in GF2XFactoring.h

    cout<<"Factors="<<factors<<endl;
 
    for(int i=0;i<factors.length();i++)
    {
     if(IterIrredTest(factors[i].a) && (deg(factors[i].a)==n)) 
     {
       cout<<"Irreducible factor of degree n:-"<<factors[i].a<<endl;
       return g;
     }
   }
  
 }

}


int main()
{
  ffs_poly f,g;
  long n=607;
  string pol_f="4,0,3,0,0,1";
  cout<<"\n polynomial f:\n";
  read_ffs_poly(f,pol_f);
  print_ffs_poly(f); 
  cout<<"\nPolynomial g:\n";  
  g=polynomial_selection(f,n); 
  print_ffs_poly(g);
  return 0;
}


