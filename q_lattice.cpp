/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */




#include <iostream>
#include<NTL/tools.h>
#include<NTL/vector.h>
#include<NTL/matrix.h>
#include<NTL/ZZX.h>
#include<NTL/GF2X.h>
#include"q_lattice.h"
#include"ffsnorm.h"


using namespace std;
using namespace NTL;
#define MAX_PREC_N 32
#define SCALE 0 

void ffs_poly_norm(GF2X& norm, ffs_poly &ffspol, GF2X& a, GF2X& b);

// The skewness is a (difference of) degree -> an unsigned int.
// The function returns 1 if it succeeded, i.e. the result fits in the
// bounds for the a_i.
int skewness_Gaussian(q_lattice& q_lat, unsigned int skewness)
{
    if (!q_lat.want_long_q) 
    {
        Vec<GF2X> a;
        Vec<GF2X> b;

        a.SetLength(2);
        b.SetLength(2); 

        a[0]=q_lat.q;
        clear(b[0]);
        a[1]=q_lat.r;
        
        GF2X_set_ti(b[1], skewness);  // this is one if there is no skewness.

        do {
              GF2X qq;

              div(qq,a[0],a[1]);
      
              if ((deg(qq) + deg(b[1])) > deg(a[0]))
                 break;

              GF2X_submul(a[0], a[1], qq);
              GF2X_submul(b[0], b[1], qq);

              div(qq, a[1], a[0]);
              if ((deg(qq) + deg(b[0])) > deg(a[1]))
                  break;
              GF2X_submul(a[1], a[0], qq);
              GF2X_submul(b[1], b[0], qq);
          } while (deg(a[0]) > deg(b[0]));

        // Compensate for the skewness
        RightShift(b[0], b[0], (long)skewness);
        RightShift(b[1], b[1], (long)skewness);

        // cast and check that the result fits!
        q_lat.a0=a[0];
        q_lat.a1=a[1];
        q_lat.b0=b[0];
        q_lat.b1=b[1];
        return 1;
   
}
else
{
  cout<<"skewness Gaussian Computstion is not defined for long q\n";
  exit(EXIT_FAILURE);
}
}

void print_q_lattice_info(q_lattice& q_lat)
{
    if (!q_lat.want_long_q) {
        cout<<"# q-lattice info:\n";
        cout<<"#   q = "<<q_lat.q; 
        cout<<" ; rho = "<<q_lat.r<<endl;
        cout<<"#   a0 = "<<q_lat.a0;  
        cout<<" ; a1 = "<< q_lat.a1;
        cout<<" ; b0 = "<<q_lat.b0;
        cout<<" ; b1 = "<<q_lat.b1<<endl;
    } 
}

int is_valid_special_q(q_lattice& q_lat, ffs_poly& F)
{
    // F(rho) = Norm_F(rho, 1)
    GF2X q, rho, one, norm;
   
    if (!q_lat.want_long_q) {
        rho=q_lat.r;
        q=q_lat.q;
    } else {
        rho=q_lat.long_r;
        q=q_lat.long_q;
    }
    one=LeftShift(one, 0);

    ffs_poly_norm(norm, F,rho,one);

    rem(norm,norm,q);
    int ret = IsZero(norm);
    norm.kill();
    q.kill();
    rho.kill();
    one.kill();
    return ret;
}



// a = i*a0 + j*a1
// b = i*b0 + j*b1

void ij_to_ab(GF2X& a, GF2X& b, GF2X& i, GF2X& j, q_lattice& q_lat)
{
    GF2X tmp;     
    GF2X f;

    f=q_lat.a0 * i;
    tmp=f;
    f=q_lat.a1*j  ; 
    a=f;
    add(a,a,tmp); 

    f=q_lat.b0 * i;
    tmp=f;   
    f=q_lat.b1*j ;  
    b=f;
       
    add(b, b,tmp);

    tmp.kill();
    f.kill();
}

// i = (a*b1 - b*a1) / q
// j = (-a*b0 + b*a0) / q

// The input must be an (a,b) pair in the q-lattice, so that the
// divisions by q are exact (q is the determinant of the base-change
// matrix). If this is not the case, abort.
// In principle, this function will always be called with inputs such
// that i and j fits within their types. So we don't return an error code
// but also abort if this is not the case.

void ab_to_ij(GF2X& i, GF2X& j, GF2X& a, GF2X& b, q_lattice& q_lat)
{
    GF2X tmp, tmp2, tmpq;
    //int fit;
    
    tmpq=q_lat.q;
   
    mul(tmp, a, q_lat.b1);
    mul(tmp2, b, q_lat.a1);
    
    sub(tmp,tmp,tmp2);
    DivRem(tmp,tmp2,tmp,tmpq);

    if(!IsZero(tmp2));
     exit(EXIT_FAILURE);

    i=tmp;
    mul(tmp, a, q_lat.b0);
    mul(tmp2, b, q_lat.a0);
    
    sub(tmp,tmp2,tmp);
    DivRem(tmp,tmp2,tmp,tmpq);

     if(!IsZero(tmp2));
     exit(EXIT_FAILURE);

    j=tmp; 
    
    tmp.kill(); 
    tmp2.kill();
    tmpq.kill();
}


// Compute lambda for an element of the factor base.
// If the result is projective, the set lambda to p.
void fb_lambda_compute(GF2X& lambda,GF2X& p, GF2X& r, q_lattice& q_lat)
{
    GF2X t0, t1;
    GF2X a0, a1, b0, b1;

    a0=q_lat.a0; 
    a1=q_lat.a1;
    b0=q_lat.b0;
    b1=q_lat.b1;
   

    int was_proj = deg(r) == deg(p);

    if (was_proj)
   {
        sub(t0, r, p);
        MulMod(t0, a0, t0, p);
        sub(t0, t0, b0);
        rem(t0, t0, p);
    } else 
      {
        MulMod(t0, b0, r, p);
        sub(t0, a0, t0);
        rem(t0, t0, p);
       }
     
    if (IsZero(t0))
   {
      lambda=p;
      return;
    }
    InvMod(t0, t0, p);
    if (IsZero(t0)) 
    { 
        
        lambda=p;
        return;
    }

    if (was_proj) {
        MulMod(t1, a1, r, p);
        sub(t1, b1, t1);
        rem(t1, t1, p);
    } else {
        MulMod(t1, b1, r, p);
        sub(t1, t1, a1);
        rem(t1, t1, p);
    }
    MulMod(lambda, t0, t1, p);
}

