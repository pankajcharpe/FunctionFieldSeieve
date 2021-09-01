/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include"ffspol.h"
#include"ffstools.h"
#include<NTL/GF2X.h>
#include<cstring>
#include<string>

void read_ffs_poly(ffs_poly& f, string pol)
{
  string temp;
  char *str=new char[pol.size()+1] ;
  int j=0,ii=0;
  int length=pol.length();
  unsigned int n=0;

   while(ii<length)
   {
     if(pol[ii]==',')
      n++;
      ii++;
   }  
   n=n+1;
   f.coeffs.SetLength(n);
   f.deg=n-1;

    string::iterator i;
     i=pol.begin();
  

  while(1)
  {
      while(*i!=',' && i!=pol.end())
     {
       temp+=*i;
       i++;
     }   
     strcpy (str, temp.c_str());

     ZZ bin_num=hex_to_binary(str);
     f.coeffs[j]=set_GF2X(bin_num);
     j++;
      if(i==pol.end())
     break;
     i++;
     temp.clear();
    
  }
   
}


void print_ffs_poly(const ffs_poly& f)
{
  cout<<"\n"<<f.coeffs<<endl;
}


//ffs polynomial evaluation using efficient horner's method 
//y=f(x)

void ffs_poly_evaluation(GF2X& y,ffs_poly& ffspol, GF2X& x)
{
   
  long deg=ffspol.deg;
  
  if(deg==-1)     //degree(0)=-1
  {
    clear(y);
    return;
  }
   
  if(deg==0)  //case of constant polynomial
  {
     y=ffspol.coeffs[0];
     return;
  }
  GF2X temp;
 
  mul(temp,x,ffspol.coeffs[deg]);
  for(int i=deg-1;i>0;i--)
  {
    add(temp,temp,ffspol.coeffs[i]);
    mul(temp,temp,x);
  }
  add(y,temp,ffspol.coeffs[0]);
  temp.kill();    
}


//Evaluation of y=f'(x)
void ffs_poly_evaluation_diff(GF2X& y, ffs_poly& ffspol, GF2X& x)
{
   long deg=ffspol.deg;
  
  if(deg==-1)     //degree(0)=-1
  {
    clear(y);
    return;
  }
   
  if(deg==0)  //case of constant polynomial
  {
     y=ffspol.coeffs[1];
     return;
  }
   long s;
  GF2X temp,aux;
  s=deg;

  aux=ffspol.coeffs[deg]*s;
  mul(temp,x,aux);
  
  for(int i=deg-1;i>1;--i)
  {
     s=i;
     aux=ffspol.coeffs[i]*s;
     add(temp,temp,aux);
     mul(temp,temp,x);
  }
  add(y,temp,ffspol.coeffs[1]);
  temp.kill();
  aux.kill();    
}

//Multiplication of ffs polynomials
void ffs_poly_mul(ffs_poly& ffspol,ffs_poly& x,ffs_poly& y)
{
    //If any one polynomial is zero
     if(x.deg==-1 || y.deg==-1)
     {
        ffspol.deg=-1;
        return ;
     } 
      
     ffs_poly temp;
     temp.coeffs.SetLength(x.deg+y.deg+1);    
     
   //  for(int i=0;i<=(x.deg+y.deg);i++)
         // clear(temp.coeffs[i]);
     
     GF2X pol_temp;

     
     for(int i=0;i<=x.deg;i++)
      for(int j=0;j<=y.deg;j++)
       {
         mul(pol_temp,x.coeffs[i],y.coeffs[j]);
         add(temp.coeffs[i+j],temp.coeffs[i+j],pol_temp);
       }
       clear(pol_temp);
       
       ffspol.coeffs.SetLength(x.deg+y.deg+1);
       
         for(int i=0;i<=x.deg+y.deg;i++)
           ffspol.coeffs[i]=temp.coeffs[i];
           
           ffspol.deg=x.deg+y.deg;

           temp.coeffs.kill();
           temp.deg=-1;     
}


//Addition of ffs polynomial
void ffs_poly_add(ffs_poly& ffspol,ffs_poly& x,ffs_poly& y)
{
   int index;
   bool flag;
   
    if(x.deg==y.deg)
      flag=1;
    else
      flag=0;
      
    ffs_poly *xx,*yy;
           
    if(x.deg > y.deg)
    {
      xx=&y;
      yy=&x;
    }
    else
    {
      xx=&x;
      yy=&y;    
    }  
    
    
     ffspol.coeffs.SetLength(yy->deg+1);

    for(index=0;index<=xx->deg;index++)
      add(ffspol.coeffs[index],xx->coeffs[index],yy->coeffs[index]);
      
    for(;index<=yy->deg;index++)
      ffspol.coeffs[index]=yy->coeffs[index];
      
      ffspol.deg=yy->deg;
      
      if(flag==1)
        ffspol.deg=ffspol.coeffs.length()-1;
        
}


//Scalar multiplication of ffs polynomial
void ffs_poly_scalar_mul(ffs_poly& ffspol,ffs_poly& x,GF2X& y)
{
    if(IsZero(y))
    {
      ffspol.deg=-1;
      return;
    }
    
    ffspol.coeffs.SetLength(x.deg+1);
    
    for(int index=0;index<=x.deg;index++)
      mul(ffspol.coeffs[index],x.coeffs[index],y);
      
      ffspol.deg=x.deg;
 }   

//Print univariate polynomial in variable t
void print_GF2X(GF2X& f)
{
 
    if(deg(f)==-1)
     return ;
   
   if(deg(f)==0)
   {
     cout<<f;
     return;
   }
   
   unsigned long d=deg(f);
   
   for(int i=d;i>=0;i--)
   {
     GF2 term=f[i];
     
     if(i==0 && IsOne(term) )
       cout<<i;
     if(i==1 && IsOne(term))
     {
        cout<<"t";
        if(IsOne(f[0]))
         cout<<"+";
     }
     
     if(IsOne(term))
     { 
       cout<<"t^"<<i;
       cout<<"+";
     } 
   }         
        
  }















