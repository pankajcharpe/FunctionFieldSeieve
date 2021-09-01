/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   factorBase.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 12:30 PM
 */

#ifndef FACTORBASE_H
#define FACTORBASE_H


#include <iostream>
#include<NTL/tools.h>
#include<NTL/vector.h>
#include<NTL/matrix.h>
#include<NTL/GF2X.h>
#include<stdint.h>


using namespace std;
using namespace NTL;

#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include"q_lattice.h"
/******************************************************************************************/

/*----------------------------------FactorBase---------------------------------------------*/
/*
Factor Base contains prime ideals of the form (p,r).
These prime ideal can be small or large.The large ideal are to be bucket sieved.
For small ideal we can afford heavy structure,whereas for large ideals we have memory consideration  	

/******************************************************************************************/

typedef uint64_t ij_vec;
//Large Ideals

struct large_ideal
{
  GF2X p;
  GF2X r;
  GF2X lambda;
  ij_vec euclid[3];
  long data;
};

typedef struct large_ideal large_ideal_t;

struct large_factor_base
{
   unsigned alloc;
   long n;                                        //no of entries in factor base
   Vec<large_ideal_t>  elts;						
};
   
typedef struct large_factor_base large_factor_base_t;

static inline unsigned fbideal_deg(large_ideal_t *ele) {
  return ele->data & 31U;
} 


//Small Ideals

 struct small_ideal{
    GF2X q;
    GF2X r;
    GF2X lambda;
    Vec< GF2X> basis;
    Vec< GF2X> adjustment_basis;
    Vec< GF2X> projective_basis;
    GF2X       current;
    GF2X       tildep;
    long       degp;
    long       degq;
    int        proj;
    int        power;
};

typedef struct small_ideal small_ideal_t;

struct small_factor_base
{
   unsigned alloc;
   long n;                       //no of entries in factor base
   Vec<small_ideal_t>  elts;						
};
   
typedef struct small_factor_base small_factor_base_t; 


void small_factor_base_precomputation(small_factor_base_t& FB,
                               unsigned I, unsigned J, q_lattice& qlat);
 void normalized_echelon_multiples( Vec< GF2X>& basis, GF2X& p, long degp, int J);
 void push_small_ideal(small_factor_base_t& FB, GF2X& p, GF2X& r,
    unsigned degp, int power, unsigned I, unsigned J);
 void push_large_ideal(large_factor_base_t& FB, GF2X& p, GF2X& r,
    unsigned degp);
 void push_ideal(large_factor_base_t& LFB, small_factor_base_t& SFB,
    GF2X& p, GF2X& r, unsigned degp, int power, unsigned min_degp,
    unsigned I, unsigned J);
long factor_base_max_degp(large_factor_base_t& FB);
double expected_hit_number(large_factor_base_t& LFB,
    unsigned I, unsigned J);

int factor_base_init(large_factor_base_t& LFB, small_factor_base_t& SFB,
    const char *filename, long sorted_min_degp, long max_degp,
    unsigned I, unsigned J);

void factor_base_clear(large_factor_base_t& LFB, small_factor_base_t& SFB);


#endif /* FACTORBASE_H */

