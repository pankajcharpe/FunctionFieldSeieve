/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Sieve.cpp
 * Author: pankaj
 *
 * Created on 20 February, 2016, 3:10 PM
 */


#include <iostream>

#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include<NTL/GF2XFactoring.h>
#include<time.h>
#include<assert.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#include"factorBase.h"
#include"poly_q.h"
#include"q_lattice.h"
#include"sub_lattice.h"
#include"bucketInfo.h"
#include"macros.h"
#include"bucketInfo.h"
#include"latticeSieve.h"
#include"smoothness.h"
#include"ffsnorm.h"
#include"ffspol.h"
#include<cstring>
#include<string>

#define SCALE 0

#define ASSERT_ALWAYS(x)						\
    do {								\
        if (!(x)) {							\
            croak__("code BUG() : condition " #x " failed",		\
                    "Abort");						\
            abort();							\
        }								\
    } while (0)



using namespace std;
using namespace NTL;



int factor_survivor(GF2X& a, GF2X& b,
        MAYBE_UNUSED ijpos_t pos, 
        MAYBE_UNUSED void *buckets,
        MAYBE_UNUSED large_factor_base_t* FB,
        ffs_poly* F, int *B, q_lattice& qlat) 
{
    GF2X Nab;
    
    for (int twice = 0; twice < 2; twice++) {
        ffs_poly_norm(Nab, F[twice], a, b);
        if (qlat.side == twice) {
            GF2X qq;
            
            if (!qlat.want_long_q)
                qq=qlat.q;
            else
                qq=qlat.long_q;
            div(Nab, Nab, qq);
            qq.kill();
        }
#ifdef BUCKET_RESIEVE
        bucket_apply_at_pos(Nab, pos, buckets[twice], FB[twice]);
#endif
        if (!GF2X_is_smooth(Nab, B[twice])) {
            Nab.kill();
            return 0;
        }
    }

    {
        GF2X g;
        
        GCD(g, a, b);
        int dg = deg(g);
        g.kill();
        if (dg > 0) {
            Nab.kill();
            return 0;
        } 
    }
// 
  //  if (!fppol_is_monic(b)) {
    //    GF2 lc;
      //  lc=coeff(b, deg(b));
        //div(b, b, lc);
        //div(a, a, lc);
   // }
//    
    Vec< Pair< GF2X,long > >factors_1;
    Vec< Pair< GF2X,long > >factors_2;    
    Vec< Pair< GF2X,long > > factors;
  
    for (int twice = 0; twice < 2; twice++) {
        ffs_poly_norm(Nab, F[twice], a, b);
      
        if(twice==0){
         CanZass(factors_1,Nab,(long)0);
         factors=factors_1;
         }
        else{
         CanZass(factors_2,Nab,(long)0);
         factors=factors_2;
         }
      
        for (int i = 0; i < factors.length(); ++i) {
            if (deg(factors[i].a) > B[twice]) {
                factors_1.kill();
                factors_2.kill();
                factors.kill();
                Nab.kill();
                return 0;
            }
        }
    }

    
    cout<<a;
    cout<<",";
    cout<<b;
    cout<<":";
    
    for (int twice = 0; twice < 2; twice++) {
        cout<<factors[twice];
        if (!twice)
            cout<<":";
    }
    cout<<"\n";
   factors_1.kill();
   factors_2.kill();
    Nab.kill();
    return 1;
}


int main()
{
    ffs_poly f[2];
    q_lattice q_lat;
    large_factor_base_t LFB[2];
    small_factor_base_t SFB[2];
    int factor_base_bound[2] = {0, 0};
    int I=0, J=0; 
    int large_prime_bound[2]= {0, 0};
    unsigned int threshold[2] = {0, 0};  
    int want_sublattice=0;
    int sq_side = 0;
    int first_sieve =0;
    GF2X q0,q1;
    int rho_given = 0;
    int skewness = 0;	 
    int gf = 0;
    int want_reliable_yield = 0;
    int want_reliable_nrels = 0;
    double reliablerange = 0.03;
    int sqt = 3;
    int bench = 0;
    int bench_end = 0;
    double bench_tot_rels = 0;
    double bench_tot_time = 0;
    int want_longq = 0;

 
   //Parameter initialization for 127 bits in hexadecimal

   char pol_f[]="4,1,1,1,1,1" ;
   char pol_g[]="7b26a5c,3a3a6a9";


   //Read Function Field Polynomials
   read_ffs_poly(f[0],pol_f);
   read_ffs_poly(f[1],pol_g);

   gf=2;
   I=9;
   J=9;
   factor_base_bound[0]=12;
   factor_base_bound[1]=12;
   large_prime_bound[0]=15;
   large_prime_bound[1]=15;
   threshold[0]=18;
   threshold[1]=18;

    q_lat.want_long_q = 0;  
   
    cout<<"Sieve Configuration :"<<endl;
    
    //Special-q ranges from num to max_limit
    ZZ num=ZZ(16384); 
    ZZ max_limit=ZZ(16640);
    long ith_bit;
    long num_bits,i=0;

    num_bits=NumBits(num);
    while(i<num_bits)
    {
       ith_bit=bit(num,i);
       SetCoeff(q0,i,ith_bit);
       i++;
    }
    
   sub_lattice *sub_lat;
   sub_lat=&no_sublat; //Ignoring the 'sublattice' option

   // Most of what we do is at the sublattice level. 
   // So we fix I and J accordingly.
   I -= sub_lat->deg;
   J -= sub_lat->deg;

   //Read factor bases
   {
      char *filename1="Aroots.2.607";
      char *filename2="Rroots.2.127";
      factor_base_init(LFB[0], SFB[0], filename1, I, factor_base_bound[0],I, J);
      factor_base_init(LFB[1], SFB[1], filename2, I, factor_base_bound[1],I, J);   
      
   }

   // Allocate storage space for the buckets.
    buckets_t buckets[2];
    buckets_init(buckets[0], I, J, expected_hit_number(LFB[0], I, J), I, 1+factor_base_max_degp(LFB[0]));
    buckets_init(buckets[1], I, J, expected_hit_number(LFB[1], I, J), I, 1+factor_base_max_degp(LFB[1]));
    ASSERT_ALWAYS(buckets[0].n == buckets[1].n);
    print_bucket_info(buckets[0], buckets[1]);
    fflush(stdout);

    void * replayable_bucket = NULL;
     // Size of a bucket region.
    ijpos_t size = bucket_region_size();

    // Allocate space for a bucket region.
    uint8_t *S;
    S = (uint8_t *) malloc(size*sizeof(uint8_t));
    ASSERT_ALWAYS(S != NULL);

    q_lat.side = sq_side; 

    double tot_time = GetTime();
    double tot_norms = 0;
    double tot_sieve = 0;
    double tot_buck_fill = 0;
    double tot_buck_apply = 0;
    double tot_cofact = 0;

    int tot_no_relations = 0;
    int tot_no_special_q = 0;
    int no_relations_per_sq = 0;

    Vec<GF2X> roots;
    roots.SetLength(f[sq_side].deg);
   
   // number of roots still to work on for current q.
    int no_of_roots = 0; 
  
    if(IterIrredTest(q0))
    {
      poly_q_info q_info;
      poly_q_info_init(q_info,q0);  
  
      poly_q fq;
      poly_q_init(fq);
  
      poly_q_set_ffs_poly(fq,f[sq_side],q_info);    
      
      no_of_roots=poly_q_roots(roots,fq,q_info);
      q_lat.q=q0;
      cout<<"############################################\n";
      cout<<"# Roots for q = ";
      cout<<q0;
      cout<<":"<<roots<<endl;
      num+=2;     
    }

   //Begin of loop over special-q's
     
    do{
       // Select next special-q
      // find next q,exit early if rho was given
       
         if (no_of_roots == 0) 
        { 
             
            if (rho_given)
                break;
            // otherwise, compute next valid q.
       do{
 
          while(num<max_limit)  
          {
             num_bits=NumBits(num);
             clear(q0);
             while(i<num_bits)
            {
              ith_bit=bit(num,i);
              SetCoeff(q0,i,ith_bit);
              i++;
            }
             
            if(IterIrredTest(q0))
            {
              break;   
            }
   
            num=num+1;
            i=0;  
          }

          poly_q_info q_info;
          poly_q_info_init(q_info,q0);  
  
          poly_q fq;
          poly_q_init(fq);
  
          poly_q_set_ffs_poly(fq,f[sq_side],q_info);    
             
          no_of_roots=poly_q_roots(roots,fq,q_info);

         }while(no_of_roots==0);

        if(num>=max_limit)
        {
          if (!bench)
          {                  
             break;
          }
        }
        cout<<"############################################\n";
        cout<<"# Roots for q = ";
        cout<<q0;
        cout<<":"<<roots<<endl;
        q_lat.q=q0;
        q_lat.r=roots[no_of_roots-1];
        no_of_roots--;
      }
      else
      {
        q_lat.r=roots[no_of_roots-1];
        no_of_roots--;
      } // end of selection of next special-q.

      tot_no_special_q++;
      
      double t_tot = GetTime();

        double t_norms = 0;
        double t_sieve = 0;
        double t_buck_fill = 0;
        double t_buck_apply = 0;
        double t_cofact = 0;
        int nrels = 0;
 
        //Reduction of q-lattice
        int no_err=skewness_Gaussian(q_lat, skewness);
        if(no_err==0)
        {
          cout<<"error in skewness guassian\n";
          exit(EXIT_FAILURE);
        }

        printf("############################################\n");
        print_q_lattice_info(q_lat);
        fflush(stdout);

        // If the reduced q-lattice is still too unbalanced, then skip it.
        // the optimal degree is ceiling( (s + deg(q))/2 ).
 
        int optimal_degree = (skewness + deg(q_lat.q) + 1) / 2;
        int sq_size=max(max(deg(q_lat.a0),deg(q_lat.a1)),max(skewness+deg(q_lat.b0),skewness+deg(q_lat.b1)));
        cout<<"#   qlattice vector degree:"<<sq_size<<endl;

         if (sq_size > optimal_degree + sqt) 
         {
                cout<<"# Special-q lattice basis is too unbalanced, let's skip it!\n";
                tot_no_special_q--;
                continue;
         }

         // Precompute all the data for small factor base elements.

         for (int i = 0; i < 2; ++i)
          small_factor_base_precomputation(SFB[i], I, J, q_lat);

        // Loop on all sublattices
        // In the no_sublat case, this loops degenerates into one pass, since
        // nb = 1.
         for (sub_lat->n = 0; sub_lat->n < sub_lat->nb; sub_lat->n++) 
      {
        // Fill the buckets.

        t_buck_fill -= GetTime();
        buckets_fill(buckets[0], LFB[0], sub_lat, I, J, q_lat);
        buckets_fill(buckets[1], LFB[1], sub_lat, I, J, q_lat);
        t_buck_fill += GetTime();

        // j0 is the first valid line in the current bucket region.

        GF2X j0;
        clear(j0);
        ijpos_t pos0=0;
     
        for (unsigned k = 0; k < buckets[0].n;++k, pos0 += size)
        {
          // Skip empty bucket regions.

          if (ijvec_get_start_pos(j0, I, J) >= pos0+size)
                continue;

           // Init the bucket region
            memset(S, 0, size*sizeof(uint8_t));

           // Kill trivial positions.
              // When there are no sublattices:
              //   (i,0) for i != 1
              //   (0,j) for j != 1
              // When using sublattices, just the position (0,0)


               
                if(!k)
              {
                S[0] = 255;  // that's (0,0)
                GF2X i;
                for (i=1; ij_set_next_return(i, i, I); )
                 S[ijvec_get_offset(i, I)] = 255;
              }
              
               GF2X j;
               j=j0;
               for (int rc = 1; rc; rc = ij_monic_set_next_return(j, j, J)) {
                    if (ij_in_fp(j))                                      //need to change
                      continue;
                    long pos = ijvec_get_start_pos(j, I, J) - pos0;
                    if (pos >= size)
                      break;
                    S[pos] = 255;
                  }
             
           
              for (int twice = 0; twice < 2; twice++)
              {
                  int side = (first_sieve)?(1-twice):twice;   // Select the side to be sieved
                
                 // Norm initialization.
                // convention: if a position contains 255, it must stay like
                // this. It means that the other side is hopeless.
                
                t_norms -= GetTime();
                init_norms(S, f[side], I, J, j0, pos0, size,q_lat, q_lat.side == side,sub_lat, side);
                t_norms += GetTime();

                // Line sieve.

                unsigned int sublat_threshold;
                t_sieve -= GetTime();
                sieveSFB(S, &sublat_threshold, SFB[side], I, J,
                        j0, pos0, size, sub_lat);
                t_sieve += GetTime();

                // Apply the updates from the corresponding bucket.

                t_buck_apply -= GetTime();
                bucket_apply(S, buckets[side], k);
                t_buck_apply += GetTime();

                 // since (0,0) is divisible by everyone, its position might
                // have been clobbered.
 
                if(!k)
                 S[0]=255;

                 // mark survivors
                // no need to check if this is a valid position

                 for (unsigned i = 0; i < size; ++i)
                 {
                     if (S[i] > (threshold[side] + sublat_threshold)>>SCALE)
                     {
                        S[i] = 255;
                     }
                     else
                        S[i]=0;
                 }
               }

               t_cofact -= GetTime();

              // survivors cofactorization

              {
                  GF2X a,b;
                  GF2X i, j, g;
                  GF2X hati, hatj;
                  
                   


                int rci, rcj = 1;
                for (j=j0; rcj; rcj = ij_monic_set_next_return(j, j, J))
                {
                  long start = ijvec_get_start_pos(j, I, J) - pos0;
                  if (start >= size)
                    break;
                  rci = 1;
                  for (i=0; rci; rci = ij_set_next_return(i, i, I)) 
                  {
                    ijpos_t pos = start + ijvec_get_offset(i, I);
                   
                     if (S[pos] != 255) 
                     {
                        ij_convert_sublat(hati, hatj, i, j, sub_lat);
                      GCD(g, hati, hatj);
                      if (deg(g) != 0 && deg(hati)>0  && deg(hatj)>0)
                        continue;
                      ij_to_ab(a, b, hati, hatj, q_lat);
                      nrels += factor_survivor(a, b, pos, replayable_bucket,
                              LFB, f, large_prime_bound, q_lat);

                     }
                   }
                  }

                  j0=j;
                  a.kill();
                  b.kill(); 
                }
                t_cofact += GetTime();
              }
            }// End of loop on sublattices.

           cout<<"# Total for this special-q:"<<nrels<<" "<<"relations found in "<<t_tot<<"s\n";
  	   cout<<"Time for main steps:\n";
           cout<<t_norms<<"s   (norms);\n";
           cout<<t_sieve<<"s   (sieve);\n";
           cout<<t_buck_fill<<"+"<<t_buck_apply<<"s (buckets: fill+apply);\n";
           cout<<t_tot/nrels <<"s (cofact).\n";
           fflush(stdout);
          
        tot_no_relations += nrels;
        tot_norms   += t_norms;
        tot_sieve   += t_sieve;
        tot_buck_apply += t_buck_apply;
        tot_buck_fill += t_buck_fill;
        tot_cofact  += t_cofact; 

        if (nrels == 0) {
            no_relations_per_sq++;
        } 
        if (want_longq)
            break;

   }while(1); //End loop over special-q's
   
    free(S);
    factor_base_clear(LFB[0], SFB[0]);
    factor_base_clear(LFB[1], SFB[1]);
    buckets_clear(buckets[0]);
    buckets_clear(buckets[1]);
    roots.kill();

    tot_time = GetTime()-tot_time;
    cout<<"\n###### General statistics ######\n";
    cout<<"#   Total time:"<<tot_time<<"s\n";
    cout<<"# Time of main steps:"<<tot_norms<<"s      (norms);      ";
    cout<<tot_sieve<<"s (sieve);\n";
    cout<<"  #                           " ;
    cout<<tot_buck_fill<<"+"<<tot_buck_apply<<"s (buckets: fill+apply);   ";
    cout<<tot_cofact<<"s (cofact).\n";
    cout<<"#   Computed "<< tot_no_special_q<< "special-q\n";
    cout<<tot_no_relations<<"relations found ("<<((double)tot_no_relations / (double)tot_no_special_q)<<"\n";
     

  //ffs_poly_clear(f[0]);
  //ffs_poly_clear(f[1]);
}//main end
















 

