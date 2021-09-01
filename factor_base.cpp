/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:34 AM
 */

#include <iostream>
#include<fstream>
#include<limits.h>
#include"poly_q.h"
#include"ffspol.h"
#include<NTL/GF2XFactoring.h>
#include<cstring>
#include<string>

bool flag=0;
#define MAX_FFS_DEGREE 10

/********************FactorBase****************************************************************/
/* Factor Base format: q:n1,n2:r1,r2,r3
   In the (frequent) case where n1,n2=1,0 this can be abridged with:
              q: r1,r2,r3 
 Here, q is a prime ideal , ri are the corresponding roots .          

/***********************************************************************************************/ 
  
struct factor_base_entry
{
  GF2X q;
  GF2X r;              
  long n1;
  long n0;
};

typedef struct factor_base_entry fb_entry;

struct factor_base_list
{
  Vec<fb_entry> list;
  long length;
  long alloc;
};
typedef factor_base_list fb_list;

//Factor base intialization 
void factor_base_list_intialization(fb_list& fb_entry_list)
{
  fb_entry_list.list.SetLength(10);
  fb_entry_list.length=0;
  fb_entry_list.alloc=10;
}

//Free storage of factor base
void factor_base_list_clear(fb_list& fb_entry_list)
{
  fb_entry_list.list.kill();  
}

//Inserting an factor base entry in factor base list
void insert_entry(fb_list& fb_entry_list,GF2X& q,GF2X& r,long n1,long n0)
{
   if(fb_entry_list.length==fb_entry_list.alloc)
   {
      fb_entry_list.alloc+=10;
      fb_entry_list.list.SetLength(fb_entry_list.alloc);
   }
   fb_entry_list.list[fb_entry_list.length].q=q;
   fb_entry_list.list[fb_entry_list.length].r=r;
   fb_entry_list.list[fb_entry_list.length].n1=n1;
   fb_entry_list.list[fb_entry_list.length].n0=n0;  
   fb_entry_list.length++;
}

//Writing an entry to the factor base file
void write_entry_to_file(fb_list& fb_entry_list)
{
   ofstream file;

  if(flag==1)
   file.open("FB_G",ios::app);
  else
  file.open("FB_F",ios::app);

  if(!file.is_open())
  {  
    cout<<"Error while openning the file for writing.\n";
    exit(EXIT_FAILURE);
  }
  if(fb_entry_list.length==0)
     return ;
  long prev_n0=-1,prev_n1=-1;

  GF2X prev_q;
  clear(prev_q);

/**********************************************************/
  //Binary to Hexadecimal

/**********************************************************/

 for(int i=0;i<fb_entry_list.length;i++)
 {
   GF2X q,r;
   q=fb_entry_list.list[i].q;
   r=fb_entry_list.list[i].r;
   long n0=fb_entry_list.list[i].n0;
   long n1=fb_entry_list.list[i].n1;
   
   if(prev_q==q && n1==prev_n1 && n0==prev_n0)
   {
     file<<","<<r;
   }
   else
   {
      if(i>0)
        file<<"\n";
      prev_q=q;
      prev_n0=n0;
      prev_n1=n1;

      if(n1==1 && n0==0)
      {
        file<<q<<":"<<r;
      } 
      else
      {
        file<<q<<":"<<n1<<" "<<n0<<":"<<r;
      }
    }
 }
 file<<"\n";
}

   
//Linear composition :FF(x):=F(phi1*x + phi0)
void ffs_poly_linear_composition_GF2X(ffs_poly& FF,ffs_poly& F,GF2X& phi1,GF2X& phi0)
{
  int deg=F.deg;
  if(IsZero(phi1))
   exit(EXIT_FAILURE);

  ffs_poly phi,phik,temp;
 
  phi.coeffs.SetLength(2);
  phik.coeffs.SetLength(deg+1);
  temp.coeffs.SetLength(deg+1);

  phi.deg=1;
  phi.coeffs[0]=phi0;
  phi.coeffs[1]=phi1;
  phik.deg=0;

  SetCoeff(phik.coeffs[0],0);

  if(FF.coeffs.length()<deg+1)
    exit(EXIT_FAILURE);

  FF.deg=0;
  FF.coeffs[0]=F.coeffs[0];

  for(int i=1;i<=deg;i++)
  {
    ffs_poly_mul(phik,phik,phi);
    ffs_poly_scalar_mul(temp,phik,F.coeffs[i]);
    ffs_poly_add(FF,FF,temp);
  }
  
  phi.coeffs.kill();
  phik.coeffs.kill();
  temp.coeffs.kill();  
}
 
//Valuation of univariate polynomial
int poly_t_valuation(GF2X& pol_t,GF2X& temp)
{
  int valuation=0;
  if(IsZero(pol_t))
    return INT_MAX;
  GF2X XX,YY,aux;
   
  XX=pol_t;
  YY=temp;

  do
  {
    DivRem(XX,aux,XX,YY);
    if(IsZero(aux))
      valuation++;
  }while(IsZero(aux) && deg(XX)>0);
   
  XX.kill();
  YY.kill();
  aux.kill();
  return valuation;
}


//Valuation of bivariate polynomial
int ffs_poly_valuation(ffs_poly& ffspol,GF2X& p)
{
  int value=INT_MAX;
  GF2X pol_t; 
  pol_t=p;
    
  for(int k=0;k<=ffspol.deg;k++)
  {
     int valuation=poly_t_valuation(ffspol.coeffs[k],pol_t);

     if(valuation==0)
       return 0;
 
     value=value>valuation?valuation:value;
  }
  return value;
}

//Calculation for the lifted root
void find_lifted_root_unramified(GF2X& new_r,ffs_poly& ffspol,GF2X& r,GF2X& p,long k_max)
{
  GF2X temp1,temp2,P,R;
  long k=1;

  P=p;
  R=r;
  
  while(k < k_max)
  {
     if(2*k<=k_max)
      mul(P,P,P);  //P^2k
     else
     {
       for(long i=k+1;i<=k_max;i++)
         mul(P,P,p);         
     }
     
    ffs_poly_evaluation(temp1,ffspol,R);
    ffs_poly_evaluation_diff(temp2,ffspol,R);
    rem(temp2,temp2,P);
    InvMod(temp2,temp2,P);
    mul(temp1,temp1,temp2);
    sub(temp1,R,temp1);
    rem(R,temp1,P);
    k=k*2;   
  }
  new_r=R;
  temp1.kill();
  temp2.kill();
  P.kill();
  R.kill();
}

//Find all the affine roots of a polynomial modulo p
int find_all_affine_roots(fb_list& fb_entry_list,ffs_poly& ffspol, GF2X& p,long k_max,long k0,long m,GF2X& phi_1,GF2X& phi_0)
{
 
  if(k0>=k_max)
   return -1;
     
  int projective;
  
  Vec<GF2X> roots;
  roots.SetLength(MAX_FFS_DEGREE);
  
  poly_q_info q_info;
  poly_q_info_init(q_info,p);  
  
  poly_q f;
  poly_q_init(f);
  
  poly_q_set_ffs_poly(f,ffspol,q_info);  
 
  if(f.degree!=ffspol.deg)
  {
    projective=1;
  }
  else
  {
    projective=0;
  }
      
  int no_of_roots=poly_q_roots(roots,f,q_info);
  for(int k=0;k<no_of_roots;k++)
  {
    int multiplicity=poly_q_root_multiplicity(f,roots[k],q_info);
        if(multiplicity==1)
        {
          //k_max=1 and k=0,no lift is needed
          if(k_max==1 && k==0)
            insert_entry(fb_entry_list,p,roots[k],1,0);
          else
          {
             GF2X new_r,phi_r,temp;
             find_lifted_root_unramified(new_r,ffspol,roots[k],p,k_max-k0);
             mul(phi_r,phi_1,new_r);
             add(phi_r,phi_r,phi_0);
             SetCoeff(temp,0);
             
             for(int j=0;j<m;j++)
              mul(temp,temp,p);
              
             for(int j=1;j<=k_max-k0;j++)
             {
               mul(temp,temp,p);
               GF2X temp2;
               rem(temp2,phi_r,temp);
               insert_entry(fb_entry_list,temp,temp2,k0+j,k0+j-1);
             }                       
           }
        } 
       else
       {
          GF2X r;
          r=roots[k];
          ffs_poly FF;
          FF.coeffs.SetLength(ffspol.coeffs.length());
          ffs_poly_linear_composition_GF2X(FF,ffspol,p,r);

          int valuation=ffs_poly_valuation(FF,p);
           GF2X pmp1;
           SetCoeff(pmp1,0);
           for(int j=0;j<m+1;j++)
             mul(pmp1,pmp1,p);
             
           GF2X phi_r;
           mul(phi_r,phi_1,r);
           add(phi_r,phi_r,phi_0);
           rem(phi_r,phi_r,pmp1);
           insert_entry(fb_entry_list,pmp1,phi_r,k0+valuation,k0);
           GF2X new_phi1,new_phi0;
           mul(new_phi1,phi_1,p);
           mul(new_phi0,phi_1,r);
           add(new_phi0,new_phi0,phi_0);
           {
             GF2X pv;
             set(pv);
             LeftShift(pv,pv,0);
             for(int j=0;j<valuation;j++)
              mul(pv,pv,p);
             
             for(int j=0;j<=FF.deg;j++)
              div(FF.coeffs[j],FF.coeffs[j],pv);
              
             pv.kill();        
          }        
          find_all_affine_roots(fb_entry_list,FF,p,k_max,k0+valuation,m+1,new_phi1,new_phi0);
          FF.coeffs.kill();
          FF.deg=-1;
        }
       }
       poly_q_clear(f);
       return !projective;   
}
 
//Find all roots such that f(r)=0 mod p
void find_all_roots_ffsmodp(fb_list& fb_entry_list,ffs_poly& ffspol,GF2X& p,long power_bound)
{
  
  long k_max=power_bound/deg(p);
  if(k_max==0)
   k_max=1;
   
  GF2X phi_1,phi_0;
  clear(phi_0);  
  SetCoeff(phi_1,0);
 
 
  long projective=find_all_affine_roots(fb_entry_list,ffspol,p,k_max,0,0,phi_1,phi_0);
     


  if(!projective)
  {
    long deg_f=ffspol.deg;
    ffs_poly ff;
    ff.coeffs.SetLength(deg_f+1);

    GF2X temp;
    set(temp);
    LeftShift(temp,temp,0);

    for(long i=0;i<=deg_f;i++)
    {
       mul(ff.coeffs[i],ffspol.coeffs[deg_f-i],temp);
     
       if(i<deg_f)
        temp=temp*p;
    }
     
    ff.deg=deg_f;
    int valuation=ffs_poly_valuation(ff,p); 
    
    if(valuation<0)
     exit(EXIT_FAILURE);
     
    insert_entry(fb_entry_list,p,p,valuation,0);  
    temp=p;
    
    for(int i;i<valuation;i++)
      temp=temp*p;
      
    for(int i=0;i<=deg_f;i++)
      div(ff.coeffs[i],ff.coeffs[i],temp);
      
     fb_list new_list;
     factor_base_list_intialization(new_list);
     
     find_all_affine_roots(new_list,ff,p,k_max-1,0,0,phi_1,phi_0);
         

     //convert back to roots
     for(int i=0;i<new_list.length;i++)
     {
       GF2X new_p,new_r;
       mul(new_p,new_list.list[i].q,new_p);
       mul(new_r,new_list.list[i].r,new_p);
       add(new_r,new_r,new_p);
       insert_entry(fb_entry_list,new_p,new_r,new_list.list[i].n1+valuation,new_list.list[i].n0+valuation);
     }
     factor_base_list_clear;
     temp.kill();
     ff.deg=-1;
     ff.coeffs.kill();
   }      
      
}


//calculation of alpha for selected p
double alpha_for_p(fb_list& fb_entry_list,GF2X& p,long power_bound)
{
   Vec<double> powers;
   
   if(powers.length()==0)
   {
     powers.SetLength(20);
     powers[0]=1.0;
     powers[1]=1.0/(double)2;
     
     for(int i=2;i<20;i++)
     {
       powers[i]=powers[i-1]*powers[1];
     }  
   }    
    
   unsigned long d=deg(p);
   double alpha=0.0;
  
  //Initialize with value for rational polynomial(up to the power_bound)
   alpha=alpha+powers[d];
   
   for(unsigned long i=2;i*d<power_bound;i++)
    alpha=alpha+powers[i*d]; 
    
    //Subtract contribution of each entry
      if(fb_entry_list.length!=0)
      {
          double alpha_temp=0.0;
          for(unsigned long i;i<fb_entry_list.length;i++)
           alpha_temp=alpha_temp-(fb_entry_list.list[i].n1-fb_entry_list.list[i].n0)*powers[deg(fb_entry_list.list[i].q)];
           
            alpha=alpha+alpha_temp/(powers[d]+1);   //for (a,b) pair to be coprime
       }     
    
       //normalize alpha
       alpha*=d;  //multiply the valuation by the degree
       return alpha;
}


//Factor base creation
void build_factor_base(ffs_poly ffspol,int factor_base_bound,int power_bound)
{
  cout<<"\nFactor Base is loading:\n";
  cout<<"Factor base bound="<<factor_base_bound<<endl;
    
  double alpha=0.0;
  GF2X p;
  ZZ max_limit;
  long ith_bit;
  long num_bits,i=0;
  ZZ num;
   
     num=2; 
     max_limit=power2_ZZ(factor_base_bound+1);
     //max_limit=118;
     while(num<max_limit)  
    {
       num_bits=NumBits(num);
        while(i<num_bits)
        {
          ith_bit=bit(num,i);
          SetCoeff(p,i,ith_bit);
          i++;
        }
             
       if(IterIrredTest(p))
      {
         
         fb_list fb_entry_list;
         factor_base_list_intialization(fb_entry_list);
         fb_entry_list.length=0;
         
         find_all_roots_ffsmodp(fb_entry_list,ffspol,p,power_bound);
         
         alpha+=alpha_for_p(fb_entry_list,p,power_bound);
         
         write_entry_to_file(fb_entry_list);
         factor_base_list_clear(fb_entry_list);
      }
       i=0;
     
       if(num==2)       
          num+=1;
       else 
       num+=2; //only for odd cases
        
  }
   cout<<"alpha="<<alpha<<endl;
} 
   
  
int main(int argc,char **argv)
{
 //the FFS polynomials
  ffs_poly f,g;

//max degree of the irreducible ideals
  long degree_bound;

//max degree of the powers of ideals
  long power_bound;

//factor base bound on both the sides
  int factor_base_bound[2]={0,0};

//which side 0 or 1
  int side=0;
 

//Reading Function Fields Polynomials

string pol_f="7,1,1,1,1,1" ;//for 313 bits in hexadecimal
string pol_g="c18ee2f3ac15c0f1,b45a8852c772037f";

cout<<"\nReading of polynomial f:\n";
read_ffs_poly(f,pol_f);
print_ffs_poly(f);
cout<<"\n";
build_factor_base(f,18,1); 
 
cout<<"\nReading of polynomial g:\n";
fflush(stdin);
flag=1;
read_ffs_poly(g,pol_g);
print_ffs_poly(g);
cout<<"\n";
build_factor_base(g,18,1);
cout<<"\nTotal time: "<<GetTime()<<"secs"<<"\n\n";
return 0;
}
