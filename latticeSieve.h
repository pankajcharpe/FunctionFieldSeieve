/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   latticeSieve.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 3:13 PM
 */

#ifndef LATTICESIEVE_H
#define LATTICESIEVE_H

#include"ij_vector.h"
static int next_projective_j(GF2X& rj, GF2X& j, Vec<GF2X>& basis, int degp, int J);
static inline void sieve_hit(uint8_t *S, uint8_t scaledegp, ijpos_t pos,ijpos_t pos0);

void sieveSFB(uint8_t *S, unsigned int *thr,small_factor_base_t& FB, unsigned I, unsigned J,
    GF2X&  j0, ijpos_t pos0, ijpos_t size, sub_lattice *sublat);


#endif /* LATTICESIEVE_H */

