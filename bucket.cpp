/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include <iostream>
#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include<assert.h>
#include"gray.h"
#include"bucketInfo.h"
#include"ij_vector.h"
#include"macros.h"

using namespace std;
using namespace NTL;


// Pack an update.
static inline update_packed_t update_pack(update_t update)
{ 
  update_conv_t conv = { .update = update }; 
  return conv.packed;
}

// Unpack an update.
static inline update_t update_unpack(update_packed_t packed)
{ 
  update_conv_t conv = { .packed = packed };
  return conv.update; 
}




static inline update_t update_set_hint(hint_t hint)
{ 
   ASSERT(hint < (hint_t)1<<UPDATE_HINT_BITS);
   return (update_t)hint << UPDATE_POS_BITS;
}


static inline update_t update_adj_pos(update_t update, pos_t pos)
{
  ASSERT(pos < (pos_t)1<<UPDATE_POS_BITS);
  return update | pos; }


static inline update_t update_set(pos_t pos, hint_t hint)
{ 
 return update_adj_pos(update_set_hint(hint), pos); 
}

// Retrieve the position from an update.
static inline pos_t update_get_pos(update_t update)
{ 
 return update & (((pos_t)1<<UPDATE_POS_BITS)-1); 
}


static inline hint_t update_get_hint(update_t update)
{  
 return (update >> UPDATE_POS_BITS) & (((hint_t)1<<UPDATE_HINT_BITS)-1); 
}


// Initialize structure and allocate buckets.
void buckets_init(buckets_t& buckets, unsigned I, unsigned J,double hits, unsigned min_degp, unsigned max_degp)
{

    if( UPDATE_POS_BITS != 16)
    {
      cout<<"UPDATE_POS_BITS are not equal to 16\n";
      exit(EXIT_FAILURE);
    }
 
    unsigned n= 1 + ((ijvec_get_max_pos(I, J)-1) >> UPDATE_POS_BITS);
    buckets.n=n;

  // The size of the buckets is the number of hits divided by the number
  // of buckets. But we take of margin of 10%, and an offset of 1000 to
  // allow for deviations.

  unsigned max_size = (unsigned)((1.1*hits / n)+1000);
  buckets.max_size = max_size;
  buckets.min_degp = min_degp;
  buckets.max_degp = max_degp;
  buckets.start    = (update_packed_t **)malloc(n * sizeof(update_packed_t *));


  if(buckets.start == NULL)
  {
   cout<<"Bucket memory is not allocated\n";
   exit(EXIT_FAILURE);
  }
 
  for (unsigned k = 0; k < n; ++k) 
 {
    buckets.start[k]=(update_packed_t *)malloc(max_size * sizeof(update_packed_t));
    
    if(buckets.start[k] == NULL)
    {
      cout<<"Bucket start[k] memory is not allocated\n";
      exit(EXIT_FAILURE);
    }
 }	
 
  unsigned ndegp = max_degp - min_degp;

  buckets.degp_end =(update_packed_t **)malloc(ndegp * n * sizeof(update_packed_t *));
   
  if(buckets.degp_end == NULL)
  {
   cout<<"Bucket memory for degp_end field is not allocated\n";
   exit(EXIT_FAILURE);
  }
  
}

// Clean up memory.
void buckets_clear(buckets_t& buckets)
{
  for (unsigned k = 0; k < buckets.n; ++k)
    free(buckets.start[k]);
  free(buckets.start);
  free(buckets.degp_end);
}

// Return the size of a bucket region.
unsigned bucket_region_size()
{
  return 1u << UPDATE_POS_BITS;
}

// Print information about the buckets.
void print_bucket_info(buckets_t& buckets0,  buckets_t& buckets1)
{

  cout<<"#   nb of buckets   = "<<buckets0.n<<endl;
 cout<<"#   size of buckets on side 0 = "<<buckets0.max_size<<endl;
  cout<<"#   size of buckets on side 1 = "<<buckets1.max_size<<endl;
  cout<<"#   size of bucket-region (ie, 2^UPDATE_POS_BITS) = "<<bucket_region_size()<<endl;
  cout<<"#   number of bits for the hint = "<<UPDATE_HINT_BITS<<endl;
  cout<<"#   bit-size of a bucket-update = "<<UPDATE_BITS<<" (rounded to "<<8*sizeof(update_packed_t)<<")"<<endl;
        
}

/* Array of buckets: filling buckets.
 *****************************************************************************/

// Push an update into the corresponding bucket.
static inline void buckets_push_update(buckets_t& buckets,update_packed_t **ptr, hint_t hint,
                         ij_vec& v, unsigned I, unsigned J)
{
  ijpos_t  pos = ijvec_get_pos(v, I, J);
  size_t   k   = pos >> UPDATE_POS_BITS;
  pos_t    p   = pos & (((pos_t)1<<UPDATE_POS_BITS)-1);
  
  if((ptr[k]-buckets.start[k])>=buckets.max_size)
  {
   cout<<"Difference between ptr[k] and buckets.start[k] is greater than or equal to bucket max_size\n";
   exit(EXIT_FAILURE);

  }
  *ptr[k]++ = update_pack(update_set(p, hint));
}


// Fill the buckets with updates corresponding to divisibility by elements of
// the factor base.
void buckets_fill(buckets_t& buckets,  large_factor_base_t& FB,
         sub_lattice *sublat, unsigned I, unsigned J,  q_lattice& q_lat)
{
  printf("\nBucket filling\n");
  printf("\nBucket 1111111111111\n");
  printf("Arrive Here");
  unsigned hatI = I + sublat->deg, hatJ = J + sublat->deg;
  
  Vec<ij_vec> basis;
  basis.SetLength(I+J);      
  //ij_vec *basis  = (ij_vec *)malloc((I+J)       * sizeof(ij_vec));
  
  Vec<ij_vec> euclid;
  basis.SetLength(I+J);
  //ij_vec *euclid = (ij_vec *)malloc((hatI+hatJ) * sizeof(ij_vec));
   
  cout<<"\nlength of basis="<<basis.length();
  cout<<"\nlength of ecluid="<<basis.length(); 

   if(buckets.min_degp < I)
  {
    cout<<"The bucket sieve requires all considered ideals to be of degree larger than I.\n";
    exit(EXIT_FAILURE);
  }

    

  // Pointers to the current writing position in each bucket.

  update_packed_t **ptr =(update_packed_t **)malloc(buckets.n * sizeof(update_packed_t *));
   if(ptr==NULL)
  {
    cout<<"Memory is not allocated for update_packed_t ptr\n";
    exit(EXIT_FAILURE);
  }

 
  for (unsigned k = 0; k < buckets.n; ++k)
    ptr[k] = buckets.start[k];
 


  // Go through the factor base by successive deg(gothp).
  // We should have no small prime in this factor base.

  large_ideal_t *gothp = &(FB.elts[0]);

  if(fbideal_deg(gothp)<buckets.min_degp)
  {
    cout<<"fbideal_deg(gothp)<buckets.min_degp\n";
    exit(EXIT_FAILURE);
  }

  unsigned i = 0;
  for (unsigned degp = buckets.min_degp; degp < buckets.max_degp; ++degp)
 {
    if(degp < I)
    {
      cout<<"degp < I\n";
      exit(EXIT_FAILURE);
    }
     
   
    // Go through each prime ideal of degree degp.

    for (; i < FB.n && fbideal_deg(gothp) == degp; ++i, ++gothp)
   {
      GF2X lambda;
      if (sublat->n == 0)
     {
        fb_lambda_compute(lambda, gothp->p, gothp->r, q_lat);

      }       

      if (lambda==gothp->p)
     {
          // This is a projective root. For the moment, we skip them.
          continue;
      }
      
      hint_t hint = 0;


      unsigned dim, euclid_dim;
      if (sublat->n == 0)
     {
        ijbasis_compute_large(basis, &dim, I, J,euclid, &euclid_dim, hatI, hatJ,
            gothp, lambda);
      }
     
     // else 
     //{
        // recover basis from the saved vectors
     //   ijbasis_complete_large(basis, &dim, I, J, gothp->euclid,hatI, hatJ);
    // }
 
      ij_vec v;
      ijvec_set_zero(v);

      //Unrolled p-ary Gray code of size ENUM_LATTICE_UNROLL.

      static const uint8_t gray[] = { GRAY(ENUM_LATTICE_UNROLL) };
      unsigned             ngray  = GRAY_LENGTH(ENUM_LATTICE_UNROLL);

      // We only need the "monic" Gray code for the first iteration. Jump
      // directly there in the array.

      unsigned gray_dim = MIN(dim, ENUM_LATTICE_UNROLL);
      unsigned i0       = ngray - GRAY_LENGTH(gray_dim);

      GF2X s, t;
      clear(t);
      int rc = dim > ENUM_LATTICE_UNROLL;
      do {
        // Inner-level Gray code enumeration: just go through the Gray code
        // array, each time adding the indicated basis vector.
/*****************************************************************************************/
        //need to modify

        for (unsigned ii = i0; ii < ngray; ++ii)
       {
          v=ijvec_add( v, basis[gray[ii]]);
          buckets_push_update(buckets, ptr, hint, v, I, J);
        }
        i0 = 0;

        // Outer-level Gray code enumeration: using ij_monic_set_next, the
        // degree of the difference with the previous one indicates which basis
        // vector should be added to the current lattice point.
        // rc is set to 0 when all vectors have been enumerated.

        s=t;
        rc = rc && ij_monic_set_next_return(t, t, dim-ENUM_LATTICE_UNROLL);
        if (rc) 
       {
          s=s-t;
          v=ijvec_add(v, basis[deg(s)+ENUM_LATTICE_UNROLL]);
          buckets_push_update(buckets, ptr, hint, v, I, J);
        }
      } while (rc);
    }
/**********************************************************************************************/
    // Mark the last position for this degree in the degp_end array.

    for (unsigned i = degp - buckets.min_degp, k = 0; k < buckets.n; ++k)
      buckets.degp_end[i*buckets.n + k] = ptr[k];
  }
 
  free(ptr);
 // free(basis);
//  free(euclid);
   printf("\nBucket filling Completed\n");
}

// Apply all the updates from a given bucket to the sieve region S.
void bucket_apply(uint8_t *S,  buckets_t& buckets, unsigned k)
{
 
  // Pointer to the current reading position in the bucket.

  update_packed_t *ptr = buckets.start[k];

  MAYBE_UNUSED ijpos_t pos0 = k*bucket_region_size();

  // Iterate through each group of updates having same deg(gothp).

  update_packed_t **end = buckets.degp_end+k;
  unsigned dd = buckets.min_degp>>SCALE;
  for (unsigned degp = buckets.min_degp; degp < buckets.max_degp;
             ++degp, dd = degp>>SCALE, end += buckets.n)
  {
    while (ptr != *end) 
    {
      update_t update = update_unpack(*ptr++);
      pos_t    pos    = update_get_pos(update);


      if (S[pos] < dd) 
        cout<<"faulty pos is "<<pos0+pos<<endl;
      if(S[pos] < dd)
      {
         cout<<"S[pos]< dd\n";
         exit(EXIT_FAILURE); 
      }
      S[pos] -= dd;
    }
  }
}









