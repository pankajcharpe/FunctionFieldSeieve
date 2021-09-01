/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */




#include <iostream>
#include<NTL/tools.h>
#include<NTL/vector.h>
#include<NTL/GF2X.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include<assert.h>
#include <stdlib.h>

//#include"i_j_vector.h"
#include"gray.h"
#include"ij_vector.h"
#include"macros.h"
/*
 * Given j is a multiple of p. 
 * Compute the next multiple of p, in lex order. 
 * The output is of degree less than J. 
 */

#define FP_CHAR 2
#define SCALE 0
#define ENUM_LATTICE_UNROLL 5

static int next_projective_j(GF2X& rj, GF2X& j, Vec<GF2X>& basis, int degp, int J)
{
    
    // First use monic_set_next() on the high part of j.

    GF2X jhi, njhi;
    RightShift(jhi, j, (long)degp);
/****************************************************************************************************/
                        //need to change  

  
   int rc = ij_monic_set_next_return(njhi, jhi, J-degp);
    if (!rc)
        return 0;

    njhi=njhi-jhi;
    int d = deg(njhi);
    rj=j+basis[d];

 
    if (d >=deg(jhi)) 
    {
        if(d <= deg(jhi)) // monic case
        {
          cout<<"degree jhi is greater than d\n";
          exit(EXIT_FAILURE);
        }
        for (int k = 2; k < FP_CHAR; ++k) 
       {
            if (d > 0)
                rj=rj+basis[d-1];
            if (d > 1)
                rj=rj+basis[d-2];
        }
    }
    return 1;
}

/**********************************************************************************/
static inline
void sieve_hit(uint8_t *S, uint8_t scaledegp, ijpos_t pos,ijpos_t pos0)
{
 
  if (S[pos] < scaledegp)
    cout<<"faulty pos is "<<pos0+pos<<endl;
  if(S[pos] < scaledegp)
  {
    cout<<"S[pos] < scaledegp\n";
    exit(EXIT_FAILURE);
  }
  S[pos] -= scaledegp;
}


void sieveSFB(uint8_t *S, unsigned int *thr,
    small_factor_base_t& FB, unsigned I, unsigned J,
    GF2X&  j0, ijpos_t pos0, ijpos_t size, sub_lattice *sublat)
{
    *thr = 0;
    for (unsigned int ii = 0; ii < FB.n; ++ii) 
    {
        small_ideal_t *gothp = &(FB.elts[ii]);
        int L = gothp->degq;
        int degp = gothp->degp;
        int scaledegp = degp >> SCALE;

     
        if((unsigned)L >= I)
        {
           cout<<"Larger primes are left to the bucket sieve.";
           exit(EXIT_FAILURE);
        }


        // projective roots are handled differently

        if (gothp->proj)
       {
          // large projective roots are just skipped.
          if ((unsigned)L > J)
            continue;

          // First time round?
          if (!pos0) 
         {
            // Find the first line to fill. If no sublat, this is zero.
            // Otherwise, there is a bit of computation.

            clear(gothp->current);
          }

          GF2X j;
          int rcj = 1;
          j=gothp->current;
          while (rcj) {
            ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
            if (start >= size)
              break;

            // Sieve the whole line

            for(int i=0; i < 1<<I; ++i)
              S[start+i] -= scaledegp;


            rcj = next_projective_j(j, j, gothp->projective_basis, L, J);
          }
          gothp->current=j; // remember the next line to sieve.
          continue;
        }

        // Only the first time round.
        if (!pos0)
        {
            clear(gothp->current);          
        }

        GF2X i, j, jj;
        int rcj = 1, degj, degjj;
        for (j=j0, degj = max(deg(j), 0); rcj; )
        {
          ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
          if (start >= size)
            break;
          i=gothp->current;
          ijpos_t pos = start + ijvec_get_offset(i, I);
          if (pos0 || pos)
            sieve_hit(S, scaledegp, pos,pos0);

          // Unrolled p-ary Gray code of size ENUM_LATTICE_UNROLL.
          static const uint8_t gray[] = { GRAY(ENUM_LATTICE_UNROLL) };
          unsigned             ngray  = GRAY_LENGTH(ENUM_LATTICE_UNROLL);

          // Just in case the dimension of the vector space is lower than
          // ENUM_LATTICE_UNROLL.
          unsigned gray_dim = MIN(I-L, ENUM_LATTICE_UNROLL);
          unsigned k0       = ngray - GRAY_LENGTH(gray_dim);

         GF2X s, t;
          clear(t);
          int rc = I-L > ENUM_LATTICE_UNROLL;

          do {
            // Inner-level Gray code enumeration: just go through the Gray code
            // array, each time adding the indicated basis vector.

            if (k0 == 0) 
           {
             
#             define DOGRAY(n)                      \
                "xorl     " #n "(%[B]), %k[i]\n\t"  \
                "subb     %[degp], (%[S],%[i])\n\t"
#             define DOALLGRAY2 DOGRAY(0)  DOGRAY(4)  DOGRAY(0)
#             define DOALLGRAY3 DOALLGRAY2 DOGRAY(8)  DOALLGRAY2
#             define DOALLGRAY4 DOALLGRAY3 DOGRAY(12) DOALLGRAY3
#             define DOALLGRAY5 DOALLGRAY4 DOGRAY(16) DOALLGRAY4
#             define DOALLGRAY6 DOALLGRAY5 DOGRAY(20) DOALLGRAY5
#             define DOALLGRAY7 DOALLGRAY6 DOGRAY(24) DOALLGRAY6
#             define DOALLGRAY8 DOALLGRAY7 DOGRAY(28) DOALLGRAY7
#             define DOALLGRAY  CAT(DOALLGRAY, ENUM_LATTICE_UNROLL)
              uint64_t ii = GF2X_to_uint64_t(i);
              uint8_t  dd = scaledegp;
             /* __asm volatile( DOALLGRAY
                            : [i]   "+r" (ii)
                            : [S]    "r" (S+start),
                              [B]    "r" (gothp->basis),
                              [degp] "r" (dd)
                            : "memory");*/
              if((ii >> 32) != 0)
              {
                cout<<"ii>>32 is not zero\n";
                exit(EXIT_FAILURE);
              }
              i = uint64_t_to_GF2X(ii);
            } 
            else 
            {
                for (unsigned k = k0; k < ngray; ++k)
               {
                    i=i+gothp->basis[gray[k]];
                    pos = start + ijvec_get_offset(i, I);
                    sieve_hit(S, scaledegp, pos, pos0);
                }
            }
            k0 = 0;


            // Outer-level Gray code enumeration: using ij_set_next, the degree
            // of the difference with the previous one indicates which basis
            // vector should be added to the current lattice point.
            // rc is set to 0 when all vectors have been enumerated.
/*******************************************************************************************/
            //need to change
           s=t;
            rc = rc && ij_set_next_return(t, t, I-L-ENUM_LATTICE_UNROLL);
            if (rc)
           {
              s=s+t;
              add(i, i, gothp->basis[deg(s)+ENUM_LATTICE_UNROLL]);
              pos = start + ijvec_get_offset(i, I);
              sieve_hit(S, scaledegp, pos,pos0);
            }
          } while (rc);

          jj=j;
          rcj = ij_monic_set_next_return(j, j, J);
          if (rcj) 
         {
            sub(jj, jj, j);
            degjj = deg(jj);
            add(gothp->current, gothp->current,
                   gothp->basis[I-L+degjj]);
            if (degjj > degj) {
              sub(gothp->current, gothp->current,
                     gothp->adjustment_basis[degjj-1]);
              degj = degjj;
            }
          }
        }
    }
}

/*******************************************************************************************/
