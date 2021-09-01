/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ffsnorm.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:54 AM
 */

#ifndef FFSNORM_H
#define FFSNORM_H


#include"ffspol.h"
#include"q_lattice.h"
#include"sub_lattice.h"
#include"ij_vector.h"

void ffs_poly_norm(GF2X& norm, ffs_poly &ffspol, GF2X& a, GF2X& b);
void ffs_poly_2ij( ffs_poly& ffspol_ij,  ffs_poly& ffspol_ab,q_lattice& q_lat);
int max_special(int prev_max, int j, int *repeated);
int deg_norm_prec_0(ffs_poly& ffspol_ij, int deg_i, int deg_j, int *gap);
void to_prec_N(GF2X& r,GF2X& p, unsigned int N);
int deg_norm_prec_N(ffs_poly &ffspol_ij, int degi, Vec<GF2X>& pow_i, int degj,Vec<GF2X>& pow_j, int *gap, int max_deg);
int deg_norm_full(ffs_poly& ffspol_ij, Vec<GF2X>& pow_i,Vec<GF2X>& pow_j, int *gap, int max_deg);
int deg_norm_ij(ffs_poly& ffspol_ij, GF2X& i, GF2X& j, int *gap);
void init_norms(uint8_t * S, ffs_poly& ffspol, unsigned I, unsigned J,GF2X& j0, ijpos_t& pos0, ijpos_t& size, q_lattice& q_lat, int sqside, sub_lattice* sub_lat,int side);


#endif /* FFSNORM_H */

