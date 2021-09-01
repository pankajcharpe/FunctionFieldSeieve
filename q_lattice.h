/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   q_lattice.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:56 AM
 */

#ifndef Q_LATTICE_H
#define Q_LATTICE_H

#include<NTL/GF2X.h>
#include"ffspol.h"
#include"sub_lattice.h"
using namespace std;
using namespace NTL;


/******************************************************************************************/

/*----------------------------------q-lattice---------------------------------------------*/

/******************************************************************************************/

struct q_lattice
{
  GF2X q;
  GF2X r;
  GF2X a0;
  GF2X a1;
  GF2X b0;
  GF2X b1;
  int side;
  int want_long_q;
  GF2X long_q;
  GF2X long_r;
  GF2X long_a0;
  GF2X long_a1;
  GF2X long_b0;
  GF2X long_b1;
};

int skewness_Gaussian(q_lattice& q_lat, unsigned int skewness);
void print_q_lattice_info(q_lattice& q_lat);
int is_valid_special_q(q_lattice& q_lat, ffs_poly& F);
void ij_to_ab(GF2X& a, GF2X& b, GF2X& i, GF2X& j, q_lattice& q_lat);
void ab_to_ij(GF2X& i, GF2X& j, GF2X& a, GF2X& b, q_lattice& q_lat);
void fb_lambda_compute(GF2X& lambda,GF2X& p, GF2X& r, q_lattice& q_lat);

#endif /* Q_LATTICE_H */

