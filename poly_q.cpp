/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



#include"poly_q.h"

/*********************************************************************************************************/

//         Information of a polynomial q

/*********************************************************************************************************/

 void poly_q_info_init(poly_q_info& q_info ,GF2X& pol_q)
{
   q_info.poly_q=pol_q;
   q_info.deg_q=deg(pol_q);
   q_info.order=1;
   
   for(long i=0;i<q_info.deg_q;i++)
      q_info.order*=2;
}

void poly_q_info_add(GF2X& res,GF2X& pol_p,GF2X& pol_q)
{
    add(res,pol_p,pol_q);
}

void poly_q_info_sub(GF2X& res,GF2X& pol_p,GF2X& pol_q)
{
    sub(res,pol_p,pol_q);
}  
    
void poly_q_info_opp(GF2X& res,GF2X& pol_p)
{
    //negate(res,pol_p);
}
  
void poly_q_info_mul(GF2X& res,GF2X& pol_p,GF2X& pol_q,poly_q_info& q_info)
{
    MulMod(res,pol_p,pol_q,q_info.poly_q);
} 

void poly_q_info_square(GF2X& res,GF2X& pol_p,poly_q_info& q_info)
{
    MulMod(res,pol_p,pol_p,q_info.poly_q);
}

void poly_q_info_inv(GF2X& res,GF2X& pol_p,poly_q_info& q_info)
{
    InvMod(res,pol_p,q_info.poly_q);
}

void poly_q_info_reduce(GF2X& res,GF2X& pol_p,poly_q_info& q_info)    
{
    rem(res,pol_p,q_info.poly_q);
}  

void poly_q_info_multi_precision(GF2X& res,GF2X& pol_t,poly_q_info& q_info)
{
  GF2X pol_p,temp;
 
  pol_p=q_info.poly_q;

  rem(temp,pol_t,pol_p);
  res=temp;
  
  temp.kill();
  pol_p.kill();
}
  
void poly_q_info_set_ti(GF2X& res,long i,poly_q_info& q_info)
{
   if(i>=q_info.deg_q)
   {
      GF2X ti,pol_q;
      set(ti);
      LeftShift(ti,ti,i);
      rem(ti,ti,pol_q);
      
      res=ti;
      ti.kill();
      pol_q.kill();
   }
   else
   {
      GF2X_set_ti(res,i);
   }
}

#define RAND_BITS 31
inline void random_poly(GF2X& f,poly_q_info& q_info)
{
  uint64_t r;
  r = 0;
    unsigned int b = 0;
    while ((b*RAND_BITS) < (unsigned int)q_info.deg_q) 
    {
        r |= ((uint64_t)rand()) << (b*RAND_BITS);
        b++;
    }
    r &= (((uint64_t)1) << (unsigned int)(q_info.deg_q))-1;
    f=uint64_t_to_GF2X(r);
}


void print_q(poly_q& q)
{
  cout<<"\n-------------poly_q------------------------\n";
  uint64_t pol;
  cout<<"Degree of poly_q="<<q.degree<<endl;
  cout<<"poly_q is=";  
   for(int i=0;i<=q.degree;i++)
  {
    pol=GF2X_to_uint64_t(q.coeff[i]);
    cout<<" "<<pol;
  }
  cout<<"\n-------------------------------------------\n";
}

void print_q_info(poly_q_info& q_info)
{
  cout<<"\n-------------poly_q_info------------------------\n";
  cout<<"poly_q="<<q_info.poly_q<<endl;
  cout<<"degree of q_info="<<q_info.deg_q<<endl;
  cout<<"order="<<q_info.order;
  cout<<"\n-------------------------------------------\n";
}
          
/*********************************************************************************************************/

//         Manipulating functions for poynomial q

/*********************************************************************************************************/

//Initialize the polynomial
void poly_q_init(poly_q& pol_q)
{
  pol_q.degree=-1;
  pol_q.alloc=0;
  pol_q.coeff.SetLength(0);
}

//Free the storage space of coefficients
void poly_q_clear(poly_q& pol_q)
{
  pol_q.coeff.kill();
}

//Reallocate space for n coefficients
inline void poly_q_realloc(poly_q& pol_q,long n)
{
   pol_q.coeff.SetLength(n);
   pol_q.alloc=n;
}

//Update the degree of the polynomial
void poly_q_update_degree(poly_q& pol_q)       
{
  while(pol_q.degree>=0 && IsZero(pol_q.coeff[pol_q.degree]))
  {
    pol_q.degree--;
  }
    
}

//set the polynomial to zero
void poly_q_set_zero(poly_q& pol_q)
{
    pol_q.degree=-1;
}

//check whether the function is zero
int poly_q_is_zero(poly_q& pol_q)
{
   if(pol_q.degree==-1)
     return 1;
   else
      return 0;
}

//set polynomial res to pol_p
void poly_q_set(poly_q& res,poly_q& pol_p,poly_q_info& q_info)
{
     poly_q_realloc(res,pol_p.degree+1);
     
     for(int i=0;i<=pol_p.degree;i++)
       res.coeff[i]=pol_p.coeff[i];
       
        res.degree=pol_p.degree;
}

//Set random polynomial of degree atmost d
void poly_q_set_random(poly_q& pol_q,long d,poly_q_info& q_info)
{
   poly_q_realloc(pol_q,d+1);
      
   for(int i=0;i<=d;i++)
   {
      random_poly(pol_q.coeff[i],q_info);
   }
   pol_q.degree=d;
   poly_q_update_degree(pol_q);
}

//set polynoomial to ffs polynomial
void poly_q_set_ffs_poly(poly_q& pol_q,ffs_poly& ffspol,poly_q_info& q_info)
{
   poly_q_realloc(pol_q,ffspol.deg+1);
   
   for(int i=0;i<=ffspol.deg;i++)
   {
     poly_q_info_multi_precision(pol_q.coeff[i],ffspol.coeffs[i],q_info);
   }
     
   pol_q.degree=ffspol.deg;
   
   poly_q_update_degree(pol_q);
}      

//Set ith coefficient of the polynomial          
void poly_q_set_ti(poly_q& pol_q,long n,poly_q_info& q_info)
{
   pol_q.coeff.SetLength(n+1);
   pol_q.alloc=n+1;
   
   for(int i=0;i<n;i++)
    clear(pol_q.coeff[i]);
    
   poly_q_info_set_ti(pol_q.coeff[n],0,q_info);
    
   pol_q.degree=n;
}

//Addition of two polynomials
void poly_q_add(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info)
{
  int i;
  
  //In case of updation of degree
  int flag=(x.degree==y.degree);
  
  poly_q xx,yy;
  
  if(x.degree >y.degree)
  {
    xx=y;
    yy=x;
  }
  else
  {
    xx=x;
    yy=y;
  }
  
  z.coeff.SetLength(yy.degree+1);
  z.alloc=yy.degree+1;
  
  for(i=0;i<=xx.degree;i++)
   add(z.coeff[i],xx.coeff[i],yy.coeff[i]);  
   
  while(i<=yy.degree)
  {
     z.coeff[i]=yy.coeff[i];
     i++;
  }     
  
  z.degree=yy.degree;
  
  if(flag==1)
   poly_q_update_degree(z);
 
 }
      
    
//Subtraction of two polynomials
void poly_q_sub(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info)             
{
   int opposite,i,flag;
   flag=(x.degree==y.degree);
   poly_q xx,yy;
   
    if(x.degree >y.degree)
  {
    xx=y;
    yy=x;
    opposite=1;
  }
  else
  {
    xx=x;
    yy=y;
    opposite=0;
  }                 
  
  z.coeff.SetLength(yy.degree+1);
  z.alloc=yy.degree+1;
  
   for(i=0;i<=xx.degree;i++)
   {
     sub(z.coeff[i],x.coeff[i],y.coeff[i]);  
   }
   while(i<=yy.degree)
  {
     if(opposite)
      z.coeff[i]=yy.coeff[i];
     else
      z.coeff[i]=-yy.coeff[i];  
     i++;
    
  }     
  
  z.degree=yy.degree;       
  if(flag==1)
   poly_q_update_degree(z);
 
}

//Multiplication by ti
void poly_q_mul_ti(poly_q& res,poly_q& p,long i,poly_q_info& q_info)
{
   if(poly_q_is_zero(p))
   {
     poly_q_set_zero(res);
     return;
   }
   
   poly_q_realloc(res,p.degree+i+1); 
   
   for(int j=p.degree;j>=0;j--)
   {
      res.coeff[i+j]=p.coeff[j];
   }  
   
   for(int j=0;j<i;j++)
    res.coeff[j]=0;
   
   res.degree=p.degree+i;
}

//Polynomial multiplication
void poly_q_mul(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info)
{
  
  if(poly_q_is_zero(x)||poly_q_is_zero(y))
  {
    poly_q_set_zero(z);
    return;
  }
  
  long d1=x.degree;
  long d2=y.degree;
  Vec<GF2X> res;
  
  res.SetLength(d1+d2+1);
  
  for(int i=0;i<(d1+d2+1);i++)
    res[i]=0;
    
    for(int i=0;i<=d1;i++)
     for(int j=0;j<=d2;j++)
      {
         GF2X temp;         
         mul(temp,x.coeff[i],y.coeff[j]);
         add(res[i+j],res[i+j],temp);                 
      }
      
      poly_q_realloc(z,d1+d2+1);
      
      for(int i=0;i<=d1+d2;i++)   
       rem(z.coeff[i],res[i],q_info.poly_q);
    
       res.kill();
       
       z.degree=x.degree+y.degree;       
}


//squaring of a polynomial
void poly_q_square(poly_q& res,poly_q& p,poly_q_info& q_info)     
{
       
       poly_q_mul(res,p,p,q_info);
}

//Scalar multiplication of a polynomial with s
void poly_q_scalar_mul(poly_q& res,poly_q& pol_p,GF2X& s,poly_q_info& q_info) 
{
  poly_q_realloc(res,pol_p.degree+1);
  
  for(int k=0;k<=pol_p.degree;k++)
  {
    MulMod(res.coeff[k],pol_p.coeff[k],s,q_info.poly_q);
  }
  
  res.degree=pol_p.degree;   

}        
         
//Scalar division of a polynomial with s
void poly_q_scalar_div(poly_q& res,poly_q& pol_p,GF2X& s,poly_q_info& q_info)
{
    if(IsZero(s))
    {
      poly_q_set_zero(res);
      return ;
    }
    
    GF2X inverse_s;
    
    InvMod(inverse_s,s,q_info.poly_q);
    
    poly_q_realloc(res,pol_p.degree+1);
    
    for(int k=0;k<pol_p.degree+1;k++)
      MulMod(res.coeff[k],pol_p.coeff[k],inverse_s,q_info.poly_q);
      
    res.degree=pol_p.degree;
}

//Get the ith coefficient of a polynomial
void poly_q_get_coeff(GF2X& coeff,poly_q& pol_p,long i,poly_q_info& q_info)
{
   coeff=pol_p.coeff[i];
}

//Set the ith coefficient of a polynomial
void poly_q_set_coeff(poly_q& res,GF2X& temp,long i,poly_q_info& q_info)
{
   long ii=i;
   
   if(ii>res.degree && IsZero(temp))
     return ;
     
   if(ii>res.degree)
   {
      poly_q_realloc(res,ii+1);
      for(int i=0;i<(ii-res.degree);i++)
         res.coeff[i]=0;
   }
    
   res.coeff[i]=temp;
   
   if(ii>res.degree)
      res.degree=ii;
   else if(IsZero(temp) && ii==res.degree)
      poly_q_update_degree(res);
}

//degree of the polynomial
long poly_q_degree(poly_q& pol_p)
{
   return pol_p.degree;
}

//Find the remainder and the quotient of the polynomail
int poly_q_divrem(poly_q& que,poly_q& rem,poly_q& x,poly_q& y,poly_q_info& q_info)
{  

  //Inverse of the leading coefficient of y
   GF2X inverse;                       
   
   long deg_x,deg_y;
   
   deg_y=y.degree;
   
   if(deg_y==-1)
    return 0;
    
   if(deg_y==0)
   {
     poly_q_get_coeff(inverse,y,0,q_info);
     poly_q_scalar_div(que,x,inverse,q_info);
     poly_q_set_zero(rem);
     return 1;
   }
  
   deg_x=x.degree;
   
   if(deg_x<deg_y)
   {
     poly_q_set_zero(que);
     poly_q_set(rem,x,q_info);
     return 1;
   }       
   
   poly_q temp_que,temp_rem;
   poly_q qq ,rr;
   
   if(&que==&x||&que==&y)
   {
      poly_q_init(temp_que);
      qq=temp_que;
   }
   else
      qq=que;
                
  if(&rem==&x||&rem==&y)
  {
      poly_q_init(temp_rem);                   
      rr=temp_rem;
  }
  else
      rr=rem;
    
  GF2X r,q;
  poly_q temp;
  poly_q_init(temp);         

  // x=y*qq+rr  
  poly_q_set(rr,x,q_info);
  poly_q_set_zero(qq);  
  poly_q_get_coeff(inverse,y,deg_y,q_info);
  InvMod(inverse,inverse,q_info.poly_q);
  long deg_rem=deg_x;  
  while(deg_rem>=deg_y)
  {
    poly_q_get_coeff(r,rr,deg_rem,q_info);
    
  
    
    MulMod(q,r,inverse,q_info.poly_q);

  
    poly_q_set_coeff(qq,q,deg_rem-deg_y,q_info);

    
    poly_q_scalar_mul(temp,y,q,q_info);
    
    poly_q_mul_ti(temp,temp,deg_rem-deg_y,q_info);  
   
    poly_q_sub(rr,rr,temp,q_info);
    
    deg_rem=rr.degree;
    
  }
  poly_q_clear(temp); 
  
  if(&que==&x||&que==&y)
  {
   poly_q_set(que,qq,q_info);
   poly_q_clear(temp_que);   
  }
  if(&rem==&x||&rem==&y)
  {
   poly_q_set(rem,rr,q_info);
   poly_q_clear(temp_rem); 
  }   
  
  
  
  return 1;
}

//Return 1 if division is possible and store quotient 
int poly_q_div(poly_q& que,poly_q& x,poly_q& y,poly_q_info& q_info)
{
  poly_q rem;
  poly_q_init(rem);
  int ret_val=poly_q_divrem(que,rem,x,y,q_info);
  poly_q_clear(rem);
  return ret_val;
}

//Return 1 if division is possible and store remainder    
int poly_q_rem(poly_q& rem,poly_q& x,poly_q& y,poly_q_info& q_info)
{
 
  poly_q que;
  poly_q_init(que);
   
  int ret_val=poly_q_divrem(que,rem,x,y,q_info);
  
  poly_q_clear(que);
 
  return ret_val;
}

//Gcd of two polynomials
void poly_q_gcd(poly_q& gcd,poly_q& x,poly_q& y,poly_q_info& q_info)
{
  if(poly_q_is_zero(x))
  {
     poly_q_set(gcd,y,q_info);
     return ;
  }
  if(poly_q_is_zero(y))
  {
     poly_q_set(gcd,x,q_info);
     return ;
  }     
  if((x.degree==0)||(y.degree==0))
  {
     poly_q_set_ti(gcd,0,q_info);
     return ;
  }
  
  poly_q p,q;
  poly_q_init(p);       
  poly_q_init(q);
  
  if(x.degree>=y.degree)
  {
    poly_q_set(p,x,q_info);
    poly_q_set(q,y,q_info);
  }
  else
  {
    poly_q_set(p,y,q_info);
    poly_q_set(q,x,q_info);
  }
  
  poly_q pp=p,qq=q,temp;
  
  while(!poly_q_is_zero(qq))
  {     
    poly_q_rem(pp,pp,qq,q_info);
        
    temp=pp;
    pp=qq;
    qq=temp;
   
  }
   
  GF2X coeff;
  
  poly_q_get_coeff(coeff,pp,pp.degree,q_info);
  poly_q_scalar_div(gcd,pp,coeff,q_info);
  poly_q_clear(p);
  poly_q_clear(q);
} 

//power modulo
void poly_q_powerMod(poly_q& res,poly_q& pol_p,uint64_t power,poly_q& pol_f,poly_q_info& q_info)
{


   if(pol_p.degree>pol_f.degree)
     exit(EXIT_FAILURE);
     
   if(power==0)
   {
     poly_q_set_ti(res,0,q_info);
     return ;
   }
   if(power==1)
   {
     poly_q_set(res,pol_p,q_info);
     return ;     
   }  
   

     
  // select msb:
    uint64_t mask = ((uint64_t)1)<<63;
    while ((power & mask) == 0)
        mask >>= 1;
    mask >>= 1;
   
   poly_q temp;
   poly_q_init(temp);
   poly_q_set(temp,pol_p,q_info);
   
   while(mask!=0)
   {

      poly_q_square(temp,temp,q_info);

 
      poly_q_rem(temp,temp,pol_f,q_info);

      if(power & mask)
      {
        poly_q_mul(temp,temp,pol_p,q_info); 
        poly_q_rem(temp,temp,pol_f,q_info);
 
      }
      mask>>=1;
      
   }
   
   poly_q_set(res,temp,q_info);
   poly_q_clear(temp);
}

//Testing for irreducibility
int poly_q_is_irreducible(poly_q& pol_p,poly_q_info& q_info)
{
  if(poly_q_is_zero(pol_p))
    exit(EXIT_FAILURE);
    
   if(pol_p.degree==0)
    return 0;
    
  if(pol_p.degree==1)
    return 1;
  
  poly_q f;
  poly_q_init(f);
  
  GF2X lead_coeff;
  poly_q_get_coeff(lead_coeff,pol_p,pol_p.degree,q_info);
  poly_q_scalar_div(f,pol_p,lead_coeff,q_info);
  
  //Compute poroduct{i<=d/2}{x^q^i-x}mod f
  poly_q xqi,x,acc,temp;
  
  poly_q_init(xqi);
  poly_q_init(x);
  poly_q_init(acc);
  poly_q_init(temp);
  
  poly_q_set_ti(x,1,q_info);
  poly_q_powerXmod(xqi,q_info.order,f,q_info);
  poly_q_sub(acc,xqi,x,q_info);
  
  int i=2;
  
  while(2*i<=pol_p.degree)
  {
     poly_q_powerMod(xqi,xqi,q_info.order,f,q_info);
     poly_q_sub(temp,xqi,x,q_info);
     poly_q_mul(acc,acc,temp,q_info);
     poly_q_rem(acc,acc,f,q_info);
     i++;
  }
  
  poly_q_gcd(acc,acc,f,q_info);
  
  int return_val;
  
  if(acc.degree>0)
   return_val=0;
  else
   return_val=1;
   
   poly_q_clear(temp);    
   poly_q_clear(acc);    
   poly_q_clear(xqi);                
   poly_q_clear(x); 
   poly_q_clear(f);
   
   return return_val;
}


//power X modulo
void poly_q_powerXmod(poly_q& res,uint64_t power,poly_q& f,poly_q_info& q_info)
{
   poly_q temp;
   poly_q_init(temp);
   poly_q_set_ti(temp,1,q_info);
   poly_q_powerMod(res,temp,power,f,q_info);
   poly_q_clear(temp);
}
  
//Splitting of a polynomial in a linear terms
void poly_q_split_linear(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info,long term)
{
  if(f.degree==1)
  {
    poly_q_get_coeff(roots[0+term],f,0,q_info);
    //negate
    roots[0+term]=-roots[0+term];
    return;
  }
  poly_q a,trace;
  poly_q_init(a);
  poly_q_init(trace);
  
  do{
     poly_q_set_random(a,f.degree-1,q_info);
     poly_q_set_zero(trace);
     for(long i=0;i<q_info.deg_q-1;i++)
     {
       poly_q_add(trace,trace,a,q_info);
       poly_q_square(a,a,q_info);
       poly_q_rem(a,a,f,q_info);
     }
     poly_q_add(trace,trace,a,q_info);
     poly_q_gcd(trace,trace,f,q_info);
    }while((trace.degree==0)||(trace.degree==f.degree));
    
   if(trace.degree==1)
   {
     poly_q_get_coeff(roots[0+term],trace,0,q_info);
     //negate
     roots[0+term]=-roots[0+term];
   }
   else
     poly_q_split_linear(roots,trace,q_info,term);
     
     term+=trace.degree;
     poly_q_div(trace,f,trace,q_info);
     
     if(trace.degree==1)
     {
        poly_q_get_coeff(roots[0+term],trace,0,q_info);
        //negate
        roots[0+term]=-roots[0+term];
     }    
     else
       poly_q_split_linear(roots,trace,q_info,term);
       
     poly_q_clear(a);
     poly_q_clear(trace);
  
}     
    
    
//This function returns the number of roots,if there is a multiple root,it is returned only once
int poly_q_roots(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info)
{ 

  if(poly_q_is_zero(f))
   exit(EXIT_FAILURE);

  if(f.degree==0)
   return 0;

  poly_q F;
  poly_q_init(F);

  GF2X lead_coeff;
  poly_q_get_coeff(lead_coeff,f,f.degree,q_info);
    
  poly_q_scalar_div(F,f,lead_coeff,q_info);
    
  if(F.degree==1)
  {
    GF2X c;
    poly_q_get_coeff(c,F,0,q_info);
    //negate
    roots[0]=-c;
    poly_q_clear(F);
    return 1;
  }
  
 //extract factors of f that contains only linear factors
  poly_q xq,x;
  poly_q_init(xq);
  poly_q_init(x);
  int ret_val; 
  poly_q_powerXmod(xq,q_info.order,F,q_info); 
  
  poly_q_set_ti(x,1,q_info); 
  poly_q_sub(xq,xq,x,q_info); 
  poly_q_gcd(xq,xq,F,q_info);
 
  if(xq.degree<=0)
   ret_val=0;
  else
  {
    ret_val=xq.degree;
    poly_q_split_linear(roots,xq,q_info,0);
  }
  poly_q_clear(xq);
  poly_q_clear(x);
  poly_q_clear(F);
  return ret_val;
}
 
//split
int poly_q_split(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info)
{
  if(poly_q_is_zero(f))
   exit(EXIT_FAILURE);

  if(f.degree==0)
   return 0;

  poly_q F;
  poly_q_init(F);

  GF2X lead_coeff;
  poly_q_get_coeff(lead_coeff,f,f.degree,q_info);
  poly_q_scalar_div(F,f,lead_coeff,q_info);
  
  if(F.degree==1)
  {
    GF2X c;
    poly_q_get_coeff(c,F,0,q_info);
    //negate
    roots[0]=-c;
    poly_q_clear(F);
    return 1;
  }
  
 //extract factors of f that contains only linear factors
  poly_q xq,x;
  poly_q_init(xq);
  poly_q_init(x);
  int ret_val;
  
  poly_q_powerXmod(xq,q_info.order,F,q_info);
  poly_q_set_ti(x,1,q_info);
  poly_q_sub(xq,xq,x,q_info);
  poly_q_gcd(xq,xq,F,q_info);
  
  if(xq.degree < f.degree )
   ret_val=0;
  else
  {
    ret_val=1;
    poly_q_split_linear(roots,xq,q_info,0);
  }
  poly_q_clear(xq);
  poly_q_clear(x);
  poly_q_clear(F);
  return ret_val;
}

       
//polynomial evaluation
void poly_q_evaluation(GF2X& y,poly_q& f,GF2X& x,poly_q_info& q_info)
{
  int deg_f=f.degree;
  if(deg_f==-1)
  {
    clear(y);
    return ;
  }
  if(deg_f==0)
  {
    poly_q_get_coeff(y,f,0,q_info);
    return;
  }
  
  GF2X temp;
  MulMod(temp,x,f.coeff[deg_f],q_info.poly_q);
  
  for(long k=deg_f-1;k>0;k--)
  {
    add(temp,temp,f.coeff[k]);
    MulMod(temp,temp,x,q_info.poly_q);
  }
  
  add(y,temp,f.coeff[0]);
}


//multiplicity of the root
int poly_q_root_multiplicity(poly_q& pol_p,GF2X& root,poly_q_info& q_info)          
{
   poly_q x_minus_root,f;
   poly_q_init(x_minus_root);
   poly_q_init(f);
   
   poly_q_set_ti(x_minus_root,1,q_info);
   GF2X mroot;
   //negate
   mroot=-root;
   poly_q_set_coeff(x_minus_root,mroot,0,q_info);
   
   int res=0;
   GF2X z;
   poly_q_set(f,pol_p,q_info);
   
   do{
       res++;
       //f=f/(x-root)
       poly_q_div(f,f,x_minus_root,q_info);
       
       //z=f(root)
       poly_q_evaluation(z,f,root,q_info);
     }while(IsZero(z)&& (f.degree>0));
     
     poly_q_clear(x_minus_root);
     poly_q_clear(f);
     return res;
}

     
   
  
  
  
  
  
  
    

       