/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include <iostream>
#include<string.h>
#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include<NTL/GF2XFactoring.h>
#include<time.h>
#include<assert.h>
#include"factorBase.h"
#include"ij_vector.h"
#include<cstring>
#include<string>

#define FP_SIZE 2


void small_factor_base_precomputation(small_factor_base_t& FB,
                               unsigned I, unsigned J, q_lattice& qlat)
{
  for (unsigned i = 0; i < FB.n; ++i) {
    small_ideal_t& gothp = FB.elts[i];
    fb_lambda_compute(gothp.lambda, gothp.q, gothp.r, qlat);
    if (gothp.lambda==gothp.q) {
      gothp.proj = 1;
    } else {
      gothp.proj = 0;
      ijbasis_compute_small(gothp.basis, gothp.adjustment_basis,
          gothp, gothp.lambda, I, J);
    }
  }
}



 void normalized_echelon_multiples( Vec< GF2X>& basis, GF2X& p, long degp, int J)
{
    if(deg(p) != degp)
     exit(EXIT_FAILURE);
    
    if (degp >= J)
        return;

    // Make it diagonal (classical reduced-echelon form).
    basis[0]=p;
    for (int i = 1; i < J-degp; ++i)
    {
        GF2X tip;
        tip=p* (long)i;
        
        basis[i]=tip;

        for (int j = i-1; j >= 0; --j) 
       {
            GF2 c;
            GF2X aux;
            c=coeff(basis[i], long(degp + j));
            aux=basis[j]*c;
            basis[i]=basis[i]-aux;
        }
    }

    // Put 1's in the lower triangle.
    for (int i = 1; i < J-degp; ++i)
        for (int j = 0; j < i; ++j)
            basis[i]=basis[i]+basis[j];
}

 void push_small_ideal(small_factor_base_t& FB, GF2X& p, GF2X& r,
    unsigned degp, int power, unsigned I, unsigned J)
{
  long size=FB.n+1;   
  FB.elts.SetLength(size);

  FB.elts[FB.n].q=p;
  FB.elts[FB.n].r=r;
  FB.elts[FB.n].degp = degp;
  FB.elts[FB.n].degq = deg(p);
  FB.elts[FB.n].power = power;
 

  FB.elts[FB.n].basis.SetLength(long(I+J));            
  FB.elts[FB.n].adjustment_basis.SetLength(long(J)); 
  FB.elts[FB.n].projective_basis.SetLength(long(J));

  
  normalized_echelon_multiples(FB.elts[FB.n].projective_basis, p,
     deg(p), J);

  // other fields are recomputed for each special-q.
 
  FB.n++;
}

 void push_large_ideal(large_factor_base_t& FB, GF2X& p, GF2X& r,
    unsigned degp)
{
  long size=FB.n+1;    
  FB.elts.SetLength(size);

  FB.elts[FB.n].p=p;
  FB.elts[FB.n].r=r;
  FB.elts[FB.n].data = degp;

  FB.n++;
}

 void push_ideal(large_factor_base_t& LFB, small_factor_base_t& SFB,
    GF2X& p, GF2X& r, unsigned degp, int power, unsigned min_degp,
    unsigned I, unsigned J)
{
  if (degp < min_degp)
    push_small_ideal(SFB, p, r, degp, power, I, J);
  else {
         if (power) 
         {
           cout<<"Warning: large power in factor base. Ignoring...\n";
           return;
         }
         push_large_ideal(LFB, p, r, degp);
       }
}



long factor_base_max_degp(large_factor_base_t& FB)
{
  return deg(FB.elts[FB.n-1].p);
}

void factor_base_clear(large_factor_base_t& LFB, small_factor_base_t& SFB)
{
  for (unsigned i = 0; i < SFB.n; ++i) {
    SFB.elts[i].basis.kill();
    SFB.elts[i].adjustment_basis.kill();
    SFB.elts[i].projective_basis.kill();
  }
  LFB.elts.kill();
  SFB.elts.kill();
}

double expected_hit_number(large_factor_base_t& LFB,
    unsigned I, unsigned J)
{  
  double powers[I+J];
  powers[0] = 1.0;
  for (unsigned i = 1; i < I+J; ++i)
    powers[i] = powers[i-1]*FP_SIZE;
  double nb = 0;
  large_ideal_t* ptr = LFB.elts.data(); 
  for (unsigned i = 0; i < LFB.n; ++i,ptr++) {
    long L = deg(ptr[0].p);
    if (L <= I+J)
      nb += powers[I+J-L];
    else
      nb += 1/powers[L-(I+J)];
  }
  return nb;
}




int factor_base_init(large_factor_base_t& LFB, small_factor_base_t& SFB,
    const char *filename, long sorted_min_degp, long max_degp,
    unsigned I, unsigned J)
{
   string line;
   char * cstr;
   ifstream myfile (filename);
 
   char str[10];
  
   ZZ bin_num;

  LFB.n=0;  
  SFB.n=0;

  
  long last_degp = 0;
  long previous_prime_degp = -1;
  int cpt = 0;
  int fbb_ok = 0;
  int count=0;
    if (myfile.is_open())
  {
    
    while ( getline (myfile,line) )
    {
      GF2X p, r;
      long degp;
      int power;
      cpt++; 
      cstr = new char [line.length()+1]; 
      strcpy (cstr, line.c_str());
      char * ptr = strtok (cstr,",: ");
      strcpy(str,ptr);
      bin_num=hex_to_binary(str);
      p=set_GF2X(bin_num);
      degp = deg(p);
      if (degp == max_degp)
         fbb_ok = 1;
      if (max_degp && degp > max_degp) {
          break;
        }
        last_degp = degp;

      ptr = strtok (NULL,",: ");
      power = 0;
    
      // short version is always prime
      degp =deg(p);
      previous_prime_degp = degp;

      while (ptr!=0)
      {
        
        strcpy(str,ptr);
        bin_num=hex_to_binary(str);
        r=set_GF2X(bin_num);
        count++;
        push_ideal(LFB, SFB, p, r, degp, power, sorted_min_degp, I, J); 
        ptr =strtok(NULL,",: ");
 
      }
       cout<<"\n";
       delete[] cstr;    
    }
      
    
    myfile.close();
  }
  else {
         cout << "Unable to open file";
         cout <<"Error in reading factor base";
         return 0;
  }
  
  cout<<"Large Factor Base";
  for(long i=0;i<=LFB.n;i++)
   {
     cout<<LFB.elts[i].p<<endl  ;
     cout<<LFB.elts[i].r<<"\n\n"; 
    
   }
   cout<<"Small Factor Base";
    for(long i=0;i<=SFB.n;i++)
   {
     cout<<SFB.elts[i].q<<endl  ;
     cout<<SFB.elts[i].r<<"\n\n"; 
    
   }
   cout<<"Count="<<LFB.n+SFB.n<<endl;
   cout<<"Actual count= "<<count;
  return 1;
}









