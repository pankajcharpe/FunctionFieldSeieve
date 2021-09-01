/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   smoothness.h
 * Author: pankaj
 *
 * Created on 21 February, 2016, 7:23 PM
 */

#ifndef SMOOTHNESS_H
#define SMOOTHNESS_H

#include <iostream>
#include<string.h>
#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include"macros.h"
#include"ffstools.h"
using namespace std;
using namespace NTL;


static int calculate_newton_indices(int *ind, int k);
void GF2X_msb_preinverse(GF2X& inv_f, GF2X& f, int k);
void GF2X_rem_precomp(GF2X& r, GF2X& pq, 
        GF2X& m, GF2X& invm);
int GF2X_is_smooth(GF2X PP, int B);



#endif /* SMOOTHNESS_H */

