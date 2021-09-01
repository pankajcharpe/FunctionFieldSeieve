/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ij_vector.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 12:28 PM
 */

#ifndef IJ_VECTOR_H
#define IJ_VECTOR_H

#include <iostream>
#include<NTL/tools.h>
#include<NTL/vector.h>
#include<NTL/matrix.h>
#include<NTL/GF2X.h>
#include<stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include"factorBase.h"
#include"ffstools.h"
using namespace std;
using namespace NTL;


/* Vector in the reduced q-lattice as a pair of polynomials (i,j).
 *  j should be kept monic.
 *****************************************************************************/
// Corresponding ij vector as an unsigned int.
typedef uint64_t ij_vec;

// Corresponding position as an unsigned int.
typedef uint64_t ijpos_t;

  void ijvec_set_i(ij_vec& v, GF2X& vec_i);
  void ijvec_set_j(ij_vec& v, GF2X& vec_j, long I);
  void ijvec_set_i_j(ij_vec& v, GF2X& vec_i, GF2X& vec_j, long I);
  void ijvec_set_zero(ij_vec& v);
  ijpos_t ijvec_get_max_pos(long I, long J);
  ijpos_t ijvec_get_pos(ij_vec& v, long I, long J);
  ijpos_t ijvec_get_start_pos(GF2X& vec_j, long I, long J);
  ijpos_t ijvec_get_offset(GF2X& vec_i, unsigned I);
  ij_vec ijvec_mul_ti(ij_vec v,long I );
  ij_vec ijvec_add(ij_vec v,ij_vec j); 
  void ij_set_ti(GF2X& f,long i);
  void ij_monic_set_next(GF2X& f,GF2X& vec_j,long j);
  int ij_monic_set_next_return(GF2X& f,GF2X& vec_j,long j);
  int ij_set_next_return(GF2X& f,GF2X& vec_i,long i);
  void ij_set_next(GF2X& f,GF2X& vec_i,long i);
  int ij_in_fp(GF2X& vec_i);
  long fill_gap(Vec<ij_vec>& v, GF2X& vec_i, GF2X& vec_j,long max_deg_i, long max_deg_j, long I);
  unsigned fill_euclid(Vec<ij_vec>& v, GF2X& vec_i, GF2X& vec_j,long max_deg_i,long max_deg_j, long I);
void ijbasis_compute_small(Vec<GF2X>& basis, Vec<GF2X>& adjustment_basis,small_ideal_t& gothp, GF2X lambda,long I, long J);
void specific_euclid_char2(Vec<ij_vec>& basis,  unsigned *basis_dim,unsigned I,unsigned  J,Vec<ij_vec>& euclid, unsigned *euclid_dim,unsigned hatI,unsigned  hatJ,GF2X& alpha_0,  GF2X& beta_0,GF2X& alpha_1,  GF2X& beta_1);
void ijbasis_compute_large(Vec<ij_vec>& basis,  unsigned *basis_dim,unsigned I,unsigned  J,Vec<ij_vec>& euclid, unsigned *euclid_dim,             unsigned hatI,unsigned  hatJ,large_ideal_t *gothp, GF2X lambda);

#endif /* IJ_VECTOR_H */

