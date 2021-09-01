/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ffspol.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:35 AM
 */

#ifndef FFSPOL_H
#define FFSPOL_H

#include<NTL/GF2X.h>
#include<NTL/vector.h>
#include<string>
#include<cstring>

using namespace std;
using namespace NTL;

//definition of ffs bivariate polynomial
struct poly
{
   long deg;
   Vec<GF2X> coeffs;
};

typedef struct poly ffs_poly;

void read_ffs_poly(ffs_poly& f, string pol);
void print_ffs_poly(const ffs_poly& f);
void ffs_poly_evaluation(GF2X& y,ffs_poly& ffspol, GF2X& x);
void ffs_poly_evaluation_diff(GF2X& y, ffs_poly& ffspol, GF2X& x);
void ffs_poly_mul(ffs_poly& ffspol,ffs_poly& x,ffs_poly& y);
void ffs_poly_add(ffs_poly& ffspol,ffs_poly& x,ffs_poly& y);
void ffs_poly_scalar_mul(ffs_poly& ffspol,ffs_poly& x,GF2X& y);
void print_GF2X(GF2X& f);

#endif /* FFSPOL_H */

