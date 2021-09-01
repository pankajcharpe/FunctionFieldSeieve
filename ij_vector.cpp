/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include"ij_vector.h"

  void ijvec_set_i(ij_vec& v, GF2X& vec_i)
{
  v=GF2X_to_uint64_t(vec_i);    
}

  void ijvec_set_j(ij_vec& v, GF2X& vec_j, long I)
{
  GF2X f=vec_j;
  f=LeftShift(f,I);
  v=GF2X_to_uint64_t(vec_j);
}

  void ijvec_set_i_j(ij_vec& v, GF2X& vec_i, GF2X& vec_j, long I)
{
   ij_vec ii,jj;
   ijvec_set_i(ii,vec_i);
   ijvec_set_j(jj,vec_j,I);
   v=ii+jj;
}

  void ijvec_set_zero(ij_vec& v)
{
   v=0;
}

/*  void ijvec_get_i(GF2X& vec_i, ij_vec v, long I)
{
   for future purpose
}

  void ijvec_get_j(GF2X& vec_j, ij_vec v, long I)
{
   for future purpose
}


  void ijvec_get_i_j(GF2X& vec_i, GF2X& vec_j, ij_vec v v, long I)
{
   for future purpose
}*/

// Return a strict higher bound on the position, given the degree bounds
// I and J.
  ijpos_t ijvec_get_max_pos(long I, long J)
{
  GF2X f;
  set(f);
  f=LeftShift(f,I+J);
  return GF2X_to_uint64_t(f);
}

// Return the position corresponding to an (i,j)-vector.
  ijpos_t ijvec_get_pos(ij_vec& v, long I, long J)
{
  return v;
}

// Return the starting position of the line corresponding to j.
// Once i is known, just add the offset ijvec_get_offset(i, I) to compute the
// full position of the vector (i,j).
  ijpos_t ijvec_get_start_pos(GF2X& vec_j, long I, long J)
{
  GF2X f;
  f=vec_j;
  f=LeftShift(f,I);
  return GF2X_to_uint64_t(f); 
}

// Return the position offset corresponding to i.
  ijpos_t ijvec_get_offset(GF2X& vec_i, unsigned I)
{
  GF2X f;
  f=vec_i;
  return GF2X_to_uint64_t(f); 
}



 ij_vec ijvec_add(ij_vec v,ij_vec j)
{
   GF2X vec=uint64_t_to_GF2X(v);
   GF2X b=uint64_t_to_GF2X(j);
   GF2X a;
   add(a,vec,j);   	 
   ij_vec res=GF2X_to_uint64_t(a); 
   return res; 
}
// Convert a position to an (i,j)-vector.
// Return 1 if successful.
/*  int ijvec_set_pos(ij_vec& v, ijpos_t pos, long I, long J)
{
   for future purpose
}*/

  ij_vec ijvec_mul_ti(ij_vec v,long I )
{
  ij_vec res;
 
  GF2X f=uint64_t_to_GF2X(v);
  f=LeftShift(f,I);
  res=GF2X_to_uint64_t(f);
  return res;
}

  void ij_set_ti(GF2X& f,long i)
{
  GF2X t;
  set(t);
  t=LeftShift(t,i);
  f=t;
  t.kill();
} 

  int ij_monic_set_next_return(GF2X& f,GF2X& vec_j,long j)
{
   uint32_t temp=GF2X_to_uint32_t(vec_j);
   temp=temp+1;
   f=uint32_t_to_GF2X(temp); 
   return (int)temp; 
}

  void ij_monic_set_next(GF2X& f,GF2X& vec_j,long j)
{
   uint32_t temp=GF2X_to_uint32_t(vec_j);
   temp=temp+1;
   f=uint32_t_to_GF2X(temp); 
}
 
  int ij_set_next_return (GF2X& f,GF2X& vec_i,long i)
{
   uint32_t temp=GF2X_to_uint32_t(vec_i);
   temp=temp+1;
   f=uint32_t_to_GF2X(temp);
   return (int)temp;

}
 
   void ij_set_next(GF2X& f,GF2X& vec_i,long i)
{
   uint32_t temp=GF2X_to_uint32_t(vec_i);
   temp=temp+1;
   f=uint32_t_to_GF2X(temp);
}

  int ij_in_fp(GF2X& vec_i)
{
  long d=deg(vec_i);

  if(d<=0)
   return 1;
  else
   return 0;
}


/* Basis of the p-lattice seen as a GF(p)-vector space of (i,j)-vectors.
 *****************************************************************************/


// Fill the basis with the vectors (i*t^k, j*t^k) as long as their degrees
// stay below the given bounds.
  long fill_gap(ij_vec *v, GF2X& vec_i, GF2X& vec_j,long max_deg_i, long max_deg_j, long I)
{
  long deg_i = deg(vec_i);
  long deg_j = deg(vec_j);
  long d_i   = max_deg_i - deg_i;
  long d_j   = max_deg_j - deg_j;
  long n    = deg_i < 0 ? d_j : min(d_i, d_j);
  if (n <= 0) 
    return 0;
 
  GF2X ii, jj; 
  ii=vec_i;
  jj=vec_j;
  ijvec_set_i_j(v[0], ii, jj, I);
  for (int k = 1; k < n; ++k)
    v[k]=ijvec_mul_ti(v[k-1], 1);
  return n;
}

  unsigned fill_euclid(ij_vec *v, GF2X& vec_i, GF2X& vec_j,long max_deg_i,long max_deg_j, long I)
{
  if ((deg(vec_i) < max_deg_i) && (deg(vec_j) < max_deg_j)) 
  {
    GF2X ii, jj;
    ii=vec_i;
    jj=vec_j;
    ijvec_set_i_j(v[0], ii, jj, I);
    return 1;
  }
  return 0;
}


// Compute the (i,j)-basis of a given p-lattice.
// Small case.
void ijbasis_compute_small(Vec<GF2X>& basis, Vec<GF2X>& adjustment_basis,
        small_ideal_t& gothp, GF2X lambda,long I, long J)
{
  long L = gothp.degq;

  // First the canonical basis:
  // Basis is { (     q*t^k,       0  ) : k in [0..I-L-1] } join
  //          { (lambda*t^k mod q, t^k) : k in [0..J-1]   }.
  // The J vectors are not stored, however.

  long k = 0;
  basis[k]=gothp.q;
 
  while (++k < I-L)
    basis[k]=LeftShift(basis[k-1], 1);
  basis[k]=lambda;

  while (++k < I+J-L)
    basis[k]=(LeftShift(basis[k-1],1)) % basis[0];

  // Transform the second part into triangular instead of diagonal (on
  // the j part) and construct the adjustment part.

  GF2 one, two;
  one=1;
  two=0;
  for (unsigned k = 0; k < J; ++k)
    adjustment_basis[k]=basis[I-L+k]*two;
  for (unsigned k = I-L+1; k < I+J-L; ++k)
    basis[k]=basis[k]+basis[k-1];
}



 void specific_euclid_char2(Vec<ij_vec>& basis,  unsigned *basis_dim,unsigned I,      unsigned  J,
                           Vec<ij_vec>& euclid, unsigned *euclid_dim,unsigned hatI,unsigned  hatJ,
                            GF2X& alpha_0,  GF2X& beta_0,
                             GF2X& alpha_1,  GF2X& beta_1)
{ 
    uint64_t v0, v1;
    uint32_t alpha0,alpha1,beta0,beta1;
    
    alpha0=GF2X_to_uint32_t(alpha_0);
    alpha1=GF2X_to_uint32_t(alpha_1);
    beta0 =GF2X_to_uint32_t(beta_0);
    beta1 =GF2X_to_uint32_t(beta_1);
 
    v0= alpha0 | ((uint64_t)beta0 << 32);
    v1 = alpha1 | ((uint64_t)beta1 << 32);

    long da0, da1;
    da0 = deg(alpha_0);
    da1 = deg(alpha_1);

    GF2X f=uint64_t_to_GF2X(v1);

    while (deg(f) < 32+J && da1 >= 0) 
    {

        // alpha0 = alpha0 mod alpha1  (and betas follow)

        uint64_t mask = ((1U<<(da0+1))-1) ^ ((1U<<da1)-1);
        uint64_t shiftv1 = v1 << (da0-da1);
        uint64_t mask1 = -(uint64_t)1;
        uint64_t maskbit = 1U<<da0;
        do {
            v0 ^= shiftv1 & mask1;
            if (!(v0 & mask))
                break;
            maskbit >>= 1;
            mask1 = (v0 & maskbit) ? (-(uint64_t)1) : 0;
            shiftv1 >>= 1;
        } while (1);

        uint64_t_swap(&v0, &v1);

        da0 = da1;
        f=uint64_t_to_GF2X(v1); 
        da1 =deg(f);

        // fill gap
        uint32_t al1, be1;
        al1 =(uint32_t)v1;
        be1 = (uint32_t)(((uint32_t)v1)+1);
  
        GF2X al_1,be_1;
        al_1=uint32_t_to_GF2X(al1);
        be_1=uint32_t_to_GF2X(be1);

        Vec<ij_vec>::iterator bs = basis.begin() ;
        Vec<ij_vec>::iterator eu = euclid.begin() ;
         
        *basis_dim  += fill_gap   (bs +*basis_dim,  al_1, be_1,
                                   min(da0, (long)I), J, I);
        //if (euclid != NULL)
            *euclid_dim += fill_euclid(eu+*euclid_dim, al_1, be_1,
                hatI, hatJ, hatI);
    }
}

// Compute the (i,j)-basis of a given p-lattice.
// Large case.

void ijbasis_compute_large(Vec<ij_vec>& basis,  unsigned *basis_dim,
                           unsigned I,      unsigned  J,
                          Vec<ij_vec>& euclid, unsigned *euclid_dim,
                           unsigned hatI,   unsigned  hatJ,
                           large_ideal_t *gothp, GF2X lambda)
{
 
  // Basis is obtained from an Euclidian algorithm on
  // (p, 0) and (lambda, 1). 

  GF2X alpha_0, beta_0, alpha_1, beta_1;
  alpha_0=gothp->p;
  clear(beta_0);
  alpha_1=lambda;
  set(beta_1);

  *basis_dim  = fill_gap(basis.begin(),  alpha_1, beta_1, I,    J,    I);
  *euclid_dim = fill_euclid(euclid.begin(), alpha_1, beta_1, hatI, hatJ, hatI);


  specific_euclid_char2(basis,  basis_dim,  I,    J,
                        euclid, euclid_dim, hatI, hatJ,
                        alpha_0, beta_0, alpha_1, beta_1);

}







