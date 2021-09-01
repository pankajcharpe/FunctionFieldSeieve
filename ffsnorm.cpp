/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


#include"ffsnorm.h"
#include"ffstools.h"

#define MAX_PREC_N 32
#define SCALE 0 


void mul_high(GF2X& r,GF2X& p,GF2X& q,unsigned int N)
{
  if(N > MAX_PREC_N && N<=0) 
   exit(EXIT_FAILURE);
  

  if(deg(p)>=(int)N || deg(q)>=(int)N)
    exit(EXIT_FAILURE);
  
  if(MAX_PREC_N > 32) // otherwise, implement another case.
    {
      cout<<"The value of MAX_PREC_N is greater than 32\n";
      exit(EXIT_FAILURE);
    }
 
 
    GF2X pp, qq;
    GF2X rr;
    pp=p;
    qq=q;
    rr=pp*qq;

    if(N>16)
   {
     r=RightShift(rr, N-1);
    }
    else
   {
     rr=RightShift(rr, N-1);
     r=rr;
   }
  

  if(deg(r)>=(int)N)
  {
    cout<<"degree of r is greater than N\n";
    exit(EXIT_FAILURE);
  }
}


/* Function computing the norm of ffs_poly at (a,b)
   norm = b^d * ffs_poly(a/b), d = deg(ffs_poly) */

void ffs_poly_norm(GF2X& norm, ffs_poly& ffspol, GF2X& a, GF2X& b)
{
   
  Vec<GF2X> pow_b;
  GF2X pow_a;
  GF2X pol_norm_i;
  GF2X tmp_norm;

  clear(pol_norm_i);
  clear(tmp_norm);
  set(pow_a);

  /* pow_b contains b^d, b^{d-1}, ... , b^2, b, 1 */
  //pow_b = (poly_t *)malloc((ffspol->degree + 1) * sizeof(poly_t));

   pow_b.SetLength(ffspol.deg+1); 
   set(pow_b[ffspol.deg]);

  for (int i = ffspol.deg - 1; i > -1; i--) 
  {
    mul(pow_b[i], pow_b[i+1], b);
  }
  for (int i = 0; i < ffspol.deg + 1; i++)
 {
    mul(pol_norm_i,ffspol.coeffs[i],pow_b[i]);
    mul(pol_norm_i,pol_norm_i,pow_a);
    add(tmp_norm, tmp_norm, pol_norm_i);
    mul(pow_a,pow_a, a);
  }
  norm=tmp_norm;
  pow_a.kill();
  pol_norm_i.kill();
  tmp_norm.kill();
  pow_b.kill();
}

/* Function computing ffs_poly_ij, a polynomial such that
   norm(ffs_poly,a,b) = norm(ffs_poly_ij, i, j), i.e.
   it is possible to apply the function norm directly on (i,j)
   with the transformed polynomial ffspol_ij, without having 
   to compute a and b with ij2ab().
*/

void ffs_poly_2ij( ffs_poly& ffspol_ij,  ffs_poly& ffspol_ab,q_lattice& q_lat)
{
   
  int d = ffspol_ab.deg;
  Vec<GF2X> powb_ij;
  GF2X tmp1, tmp2;

  powb_ij.SetLength(d+1);
  
  /* Step 0 */
  ffspol_ij.coeffs[d]=ffspol_ab.coeffs[d];
  ffspol_ij.deg = d;
  set(powb_ij[0]);
  
  for (int k = d - 1; k >= 0; --k) {
    /* For hh(i,j) * (a0 i + a1 j) */
    
      mul(ffspol_ij.coeffs[k], ffspol_ij.coeffs[k + 1], q_lat.a1);
    
    for (int l = k + 1; l < d; ++l)
    {
        mul(tmp1, ffspol_ij.coeffs[l], q_lat.a0);
        mul(tmp2,ffspol_ij.coeffs[l + 1], q_lat.a1);    
        add(ffspol_ij.coeffs[l],tmp1,tmp2);
    }
     
     mul(ffspol_ij.coeffs[d],ffspol_ij.coeffs[d], q_lat.a0);
   
  
    /* For (b0 i + b1 j)^{d-k} */
     
      mul(powb_ij[d - k], powb_ij[d - k - 1], q_lat.b0);
   
    for (int l = d - k - 1; l > 0; --l)
    {
     
        mul(tmp1,powb_ij[l - 1], q_lat.b0);
        mul(tmp2, powb_ij[l], q_lat.b1);
      
        add(powb_ij[l],tmp1,tmp2);
    }
    
      mul(powb_ij[0],powb_ij[0],q_lat.b1);
   

    /* Multiply (b0 i + b1 j)^{d-k} by f_k and add it to hh(i,j) (a0 i + a1 j) we have computed */
    for (int l = k; l <= d; ++l) {
      mul(tmp1,powb_ij[l - k],ffspol_ab.coeffs[k]);
      add(ffspol_ij.coeffs[l],ffspol_ij.coeffs[l],tmp1);
    }
  }
  powb_ij.kill();
  tmp1.kill();
  tmp2.kill();
}


int max_special(int prev_max, int j, int *repeated)
{
  
  if (prev_max == j) {
    (*repeated)++;
    return prev_max;
  } 
  else if (prev_max > j) {
      return prev_max;
  }
  else {
    *repeated = 1;
    return j;
  }
}

int deg_norm_prec_0(ffs_poly& ffspol_ij, int deg_i, int deg_j, int *gap)
{
  int degree, max_deg = -1;
  int repeated = 1;
  
  for (int k = 0; k < ffspol_ij.deg + 1; ++k) 
  {
    degree = deg(ffspol_ij.coeffs[k]) + k * deg_i + (ffspol_ij.deg - k) * deg_j;
    max_deg = max_special(max_deg, degree, &repeated);
  }
  
  /* If there is only one monomial of maximal degree */
  if (repeated == 1) 
    *gap = -1;
  else if (repeated != 1 && (repeated & 1u)) //to find cancellation
   {
    *gap = 0;
   } 

  else
    *gap = max_deg + 1;
  return max_deg;
}



/* function which takes as input a poly_t and returns a poly_t with
   N bits of precision containing the N monomials of higher degree,
   with N <= MAX_PREC_N */

void to_prec_N(GF2X& r,GF2X& p, unsigned int N)
{
  
  if(N > MAX_PREC_N && N<=0) 
   exit(EXIT_FAILURE);
  
  int shift = deg(p) - N + 1;
  if (shift == 0)
   { 
     r=p;
    }
  else 
 {
    GF2X tmp;
    tmp.SetLength(N);

    if (shift > 0) 
   {
      tmp=RightShift(p, (unsigned int) shift);
      r=tmp;
    }
    else 
    { /* we need to align the most significant bit on the left */
     tmp=LeftShift(p, (unsigned int) -shift);
     r=tmp;
    }
    tmp.kill();
  }
}


int deg_norm_prec_N(ffs_poly &ffspol_ij, int degi, Vec<GF2X>& pow_i, int degj,Vec<GF2X>& pow_j, int *gap, int max_deg)
{

  
  /* monomials contains the poly_t monomials of the norm in precision
     0 < N <= 32. The monomials and their degrees are sorted by increasing
     order of their degrees */

  int deg_norm = -1;
  unsigned int N = 0;
  const int d = ffspol_ij.deg;
  int tab_size;
  
  int degree;
  int *degrees;
  degrees = (int *)malloc((d+1) * sizeof(int));

  GF2X monomial;
  Vec<GF2X> monomials;
  
  monomials.SetLength(d+1);
  GF2X coeff_prec_N;
  GF2X pow_i_prec_N;
  GF2X pow_j_prec_N;
        
#define PREC_GAP 8
 
  do{
     N = N + PREC_GAP;

    /* Computing and sorting the degrees and the monomials in precision N.
       Also align the truncated monomials to be able to add them without
       further shifts. */   

    tab_size = 0;
    for (int k = 0; k < d+1; ++k)
    {
      degree = deg(ffspol_ij.coeffs[k]) + k * degi
        + (d - k) * degj;
      if (degree <= max_deg - (int)N) 
        continue;
      to_prec_N(coeff_prec_N,ffspol_ij.coeffs[k], N);
      to_prec_N(pow_i_prec_N, pow_i[k], N);
      to_prec_N(pow_j_prec_N, pow_j[k], N);
      mul_high(monomial, pow_i_prec_N,pow_j_prec_N, N);
      mul_high(monomial, monomial,coeff_prec_N, N);
      
      // insertion sort
      int l = tab_size;
      while (l > 0 && degrees[l - 1] > degree) {
	degrees[l] = degrees[l - 1];
	monomials[l]=monomials[l - 1];
	--l;
      }
      degrees[l] = degree;
      monomials[l]=RightShift(monomial, (max_deg - degree));
      tab_size++;
    }
    if(tab_size < 2) // otherwise, we shouldn't be here!
    {
      cout<<"\nTab_size is less than two.\n";
      exit(EXIT_FAILURE);
    }

    /* We set tab_size to the index of the last value in degrees[]
       and monomials[]. This will always be the case until the end 
       of this loop. */

    --tab_size;
    
    /* Computing the sum of all monomials of maximal degree until we
       find only one monomial (term) of maximal degree */
    
    GF2X sum;
    int deg_sum_prec; 
    do {
        sum=monomials[tab_size];
      
      /* we make the sum of all monomials of maximal degree */
      while (tab_size > 0 && degrees[tab_size] == degrees[tab_size - 1]) 
     {
	--tab_size;
	add(sum, sum,monomials[tab_size]);
      }
      
      deg_sum_prec =deg(sum);

      /* If we still have information, we put the 
	 sum in monomials[] and degree in degrees[] and keep
         the tables sorted */

      if (deg_sum_prec >= 0) 
      {
	degree = deg_sum_prec + 1 + max_deg - N;  // exact degree of the sum.
	int l = tab_size;
	while (l > 0 && degrees[l - 1] > degree) 
        {
	  degrees[l] = degrees[l - 1];
	  monomials[l]=monomials[l - 1];
	  --l;
	}
	degrees[l] = degree;
	monomials[l]=sum; 
      } 
   else {
        /* If the drop in degree is greater than the precision, we throw
           away the new monomial sum */

	--tab_size;
      }

      /* We do this until there is only one monomial of maximal degree */

    } while (tab_size > 0 && degrees[tab_size] == degrees[tab_size - 1]);

    /* Finished or have to increase the precision ? */
    if ((tab_size > 0 && degrees[tab_size] != degrees[tab_size - 1])
        || tab_size == 0) 
   {
      deg_norm = degrees[tab_size];
      *gap = max_deg - deg_norm;
      free(degrees);
      monomials.kill();
      return deg_norm;
    }

  } while (tab_size < 0 && N + PREC_GAP <= MAX_PREC_N);

  // Failed, even at maximum allowed precision. Mark it in gap.
  *gap = max_deg + 1;
  free(degrees);
  monomials.kill();
  return deg_norm;
}

#undef PREC_GAP
#undef MAX_PREC_N

int deg_norm_full(ffs_poly& ffspol_ij, Vec<GF2X>& pow_i,Vec<GF2X>& pow_j, int *gap, int max_deg)
{
   
  GF2X pol_norm_k;
  GF2X norm;
  int deg_norm = -1;
 

  clear(pol_norm_k);
  clear(norm);

  for (int k = 0; k < ffspol_ij.deg + 1; ++k) 
 {
    mul(pol_norm_k, ffspol_ij.coeffs[k], pow_j[k]);
    mul(pol_norm_k,pol_norm_k, pow_i[k]);
    add(norm,norm,pol_norm_k);
  }
  deg_norm =deg(norm);
  *gap = max_deg - deg_norm;
  
  pol_norm_k.kill();
  norm.kill();
  return deg_norm;
}

int deg_norm_ij(ffs_poly& ffspol_ij, GF2X& i, GF2X& j, int *gap)
{
   
  int degi = deg(i);
  int degj = deg(j);
  int max_deg;
  int deg_norm = -1;

  if (degi == -1)
    return (ffspol_ij.deg)*degj + deg(ffspol_ij.coeffs[0]);
  if (degj == -1)
    return (ffspol_ij.deg)*degi + deg(ffspol_ij.coeffs[ffspol_ij.deg]);

  
  /* Computation of the degree of the norm for precision 0 
   i.e. the computation works only if deg_norm is the maximal degree
   of the monomials */

  max_deg =  deg_norm_prec_0(ffspol_ij, degi, degj, gap);
  
  if (*gap != max_deg + 1)
    return max_deg;
  else {
    Vec<GF2X> pow_i;
    Vec<GF2X> pow_j;
    GF2X ii, jj;   
    ii=i;
    jj=j;  

    /* pow_j contains the powers of j in DECREASING order */

    pow_j.SetLength(ffspol_ij.deg+1);
    set(pow_j[ffspol_ij.deg]);

    for (int k = ffspol_ij.deg - 1; k > -1; --k) {
      mul(pow_j[k],pow_j[k + 1],jj);
    }

    /* pow_i contains the powers of i in INCREASING order */
    pow_i.SetLength(ffspol_ij.deg+1);
    set(pow_i[0]);

    for (int k = 1; k < ffspol_ij.deg + 1; ++k) {
      mul(pow_i[k],pow_i[k - 1],ii);
    }
    
    deg_norm = deg_norm_prec_N(ffspol_ij, degi, pow_i, degj, pow_j, gap, max_deg);
    if (*gap == max_deg + 1)
   {
      deg_norm = deg_norm_full(ffspol_ij, pow_i, pow_j, gap, max_deg);
    } 
  
     ii.kill();
     jj.kill();
     pow_i.kill();
     pow_j.kill();
    return deg_norm;
  }
}


/* Function init_norms 
   Compute the degree of the norm at each valid position in the given
   j-range.
   
   The sqside parameter is a boolean that tells whether we are on the
   side of the special q. If so, then the degree of q must be subtracted
   from the norm.
   */

void init_norms(uint8_t * S, ffs_poly& ffspol, unsigned I, unsigned J,
                GF2X& j0, ijpos_t& pos0, ijpos_t& size, q_lattice& q_lat,
                int sqside, sub_lattice* sub_lat,int side)
{
   
  ffs_poly ffspol_ij;
  ffspol_ij.coeffs.SetLength(ffspol_ij.deg+1);
 
  ffs_poly_2ij(ffspol_ij, ffspol, q_lat);  

  int degq = 0;
  if (sqside)
 { 
    degq = deg(q_lat.q);
 }

  GF2X i, j, hati, hatj;
  int gap;
  int rci, rcj = 1;
/**************************************************************/
  //need to modified  
for (j=j0; rcj; rcj = ij_monic_set_next_return(j, j, J)) {
    ijpos_t start = ijvec_get_start_pos(j, I, J) - pos0;
    if (start >= size)
      break;

    rci = 1;
    for (clear(i); rci; ) {
      ijpos_t pos = start + ijvec_get_offset(i, I);
      if (S[pos] == 255) {
        rci = ij_set_next_return(i, i, I);
        continue;
      }

      // Sublat conversion
      ij_convert_sublat(hati, hatj, i, j, sub_lat);



      // Compute the degree of the norm, and the gap information.
      int degree = deg_norm_ij(ffspol_ij, hati, hatj, &gap);
      
      // Deduce the next i for which we have to compute the norm.
      GF2X i_next;
      {
        int degi = deg(i);
        if (gap == -1) { 
          ij_set_ti(i_next, degi+1);
        } 
   else {
          int s = max(degi - gap, 0);
          RightShift(i_next, i, (long)s);
          ij_set_next(i_next, i_next, I+1); // we don't care for overflow, here
          LeftShift(i_next, i_next, (long)s);
        }
      }
/***************************************************************/
      // Fast loop with constant degree of norm.
      degree -= degq;
      degree >>= SCALE;
      if (degree > 254) {
          cout<<"Error: the scaling of norms is not enough.\n";
          exit(EXIT_FAILURE);
      }
      if (degree == 0)
        degree = 255  ;  
      
        uint8_t * Sptr = S + pos;
        uint32_t ii = GF2X_to_uint32_t(i);
        uint32_t iinext = GF2X_to_uint32_t(i_next);
        do {
          *Sptr++ |= degree;
          ii++;
        } while (ii!=iinext);
        i = uint32_t_to_GF2X(ii);
        i_next = uint32_t_to_GF2X(iinext);
        rci = iinext < (1U<<I);
      
    }
  }
  ffspol_ij.coeffs.kill(); 
}

