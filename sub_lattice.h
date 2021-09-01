/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   sub_lattice.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 8:38 PM
 */

#ifndef SUB_LATTICE_H
#define SUB_LATTICE_H

#include <iostream>
#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>

using namespace std;
using namespace NTL;

#define MAX_SUBLAT 9

#ifdef __cplusplus
extern "C" {
#endif
 
    
typedef struct 
{
    int nb;             // the number of valid sublattices
    int deg;            // the degree of the modulus
                        // The following fields are meaningless if nb = 1, deg = 0.
    int n;              // the index of the current sublattice (in .lat)
    uint16_t modulus;  // the modulus used for sublatticing
    uint16_t lat[MAX_SUBLAT][2]; // a description of all sublattices.
} sublat_struct;

typedef sublat_struct sub_lattice;

static sub_lattice no_sublat = {
    1,
    0,
    0,
    0,
    {{0},{0}}
};

static sub_lattice nine_sublat = {
    9,
    2,
    0,
     6 ,   // t^2+t = t*(t+1)
    {
        {0,1},
        {1,0}, {1,1}, {1,2}, {1,3},
        {2,1}, {2,3},
        {3,1}, {3,2},
    }
};

#ifdef __cplusplus
}
#endif

// Convert an (i,j)-vector inside a sublattice into the
// (hat i, hat j)-vector in the traditional q-lattice coordinates.
// The conversion to (a,b) is then the usual ij2ab() function.
//   hati = i*modulus + i0
//   hatj = j*modulus + j0

static inline void ij_convert_sublat(GF2X& hati, GF2X& hatj, GF2X& i, GF2X& j,
        sub_lattice *sublat)
{
      hati=i;
      hatj=j;
}

#endif /* SUB_LATTICE_H */

