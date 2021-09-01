/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   poly_q.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:37 AM
 */

#ifndef POLY_Q_H
#define POLY_Q_H

#include <iostream>
#include"ffspol.h"
#include"ffstools.h"
#include <cstdint>
using namespace std;
using namespace NTL;

//Declaration of poly_q_info which stores the information about polynomial q 
struct poly_q_info
{
   GF2X poly_q;
   long deg_q;
   uint64_t order ;  //q^deg
 };
 
 typedef  struct poly_q_info poly_q_info;
 
 //Declaration of polynomial q
 struct polynomial_q
 {     
    long degree;
    long alloc;
    Vec<GF2X> coeff;
  };
  
 typedef struct polynomial_q poly_q;                 

/*********************************************************************************************************/

//         Information of a polynomial q

/*********************************************************************************************************/  
void poly_q_info_init(poly_q_info& q_info ,GF2X& pol_q);
void poly_q_info_add(GF2X& res,GF2X& pol_p,GF2X& pol_q);
void poly_q_info_sub(GF2X& res,GF2X& pol_p,GF2X& pol_q);
void poly_q_info_opp(GF2X& res,GF2X& pol_p);
void poly_q_info_mul(GF2X& res,GF2X& pol_p,GF2X& pol_q,poly_q_info& q_info);
void poly_q_info_square(GF2X& res,GF2X& pol_p,poly_q_info& q_info);
void poly_q_info_inv(GF2X& res,GF2X& pol_p,poly_q_info& q_info);
void poly_q_info_reduce(GF2X& res,GF2X& pol_p,poly_q_info& q_info);
void poly_q_info_multi_precision(GF2X& res,GF2X& pol_t,poly_q_info& q_info); 
void poly_q_info_set_ti(GF2X& res,long i,poly_q_info& q_info);
void print_q(poly_q& q);
void print_q_info(poly_q_info& q_info);

/*********************************************************************************************************/

//         Manipulating functions for poynomial q

/*********************************************************************************************************/

//Initialize the polynomial
void poly_q_init(poly_q& pol_q);

//Free the storage space of coefficients
void poly_q_clear(poly_q& pol_q);

//Reallocate space for n coefficients
inline void poly_q_realloc(poly_q& pol_q,long n);

//Update the degree of the polynomial
void poly_q_update_degree(poly_q& pol_q);

//set the polynomial to zero
void poly_q_set_zero(poly_q& pol_q);

//check whether the function is zero
int poly_q_is_zero(poly_q& pol_q);

//set polynomial res to pol_p
void poly_q_set(poly_q& res,poly_q& pol_p,poly_q_info& q_info);

//Set random polynomial of degree atmost d
void poly_q_set_random(poly_q& pol_q,long d);

//set polynoomial to ffs polynomial
void poly_q_set_ffs_poly(poly_q& pol_q,ffs_poly& ffspol,poly_q_info& q_info);

//Set ith coefficient of the polynomial  
void poly_q_set_ti(poly_q& pol_q,long n,poly_q_info& q_info);

//Addition of two polynomials
void poly_q_add(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info);

//Subtraction of two polynomials
void poly_q_sub(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info); 

//Multiplicatiob by ti  
void poly_q_mul_ti(poly_q& res,poly_q& p,long i,poly_q_info& q_info);

//Polynomial multiplication
void poly_q_mul(poly_q& z,poly_q& x,poly_q& y,poly_q_info& q_info);

//squaring of a polynomial
void poly_q_square(poly_q& res,poly_q& p,poly_q_info& q_info);

//Scalar multiplication of a polynomial with s
void poly_q_scalar_mul(poly_q& res,poly_q& pol_p,GF2X& s,poly_q_info& q_info);

//Scalar division of a polynomial with s
void poly_q_scalar_div(poly_q& res,poly_q& pol_p,GF2X& s,poly_q_info& q_info);

//Get the ith coefficient of a polynomial
void poly_q_get_coeff(GF2X& coeff,poly_q& pol_p,long i,poly_q_info& q_info);

//Set the ith coefficient of a polynomial
void poly_q_set_coeff(poly_q& res,GF2X& temp,long i,poly_q_info& q_info);

//degree of the polynomial
long poly_q_degree(poly_q& pol_p);

//Find the remainder and the quotient of the polynomail
int poly_q_divrem(poly_q& que,poly_q& rem,poly_q& x,poly_q& y,poly_q_info& q_info);

//Return 1 if division is possible and store quotient 
int poly_q_div(poly_q& que,poly_q& x,poly_q& y,poly_q_info& q_info);

//Return 1 if division is possible and store remainder    
int poly_q_rem(poly_q& rem,poly_q& x,poly_q& y,poly_q_info& q_info);

//Gcd of two polynomials
void poly_q_gcd(poly_q& gcd,poly_q& x,poly_q& y,poly_q_info& q_info);

//power modulo
void poly_q_powerMod(poly_q& res,poly_q& pol_p,uint64_t power,poly_q& pol_f,poly_q_info& q_info);

//Testing for irreducibility
int poly_q_is_irreducible(poly_q& pol_p,poly_q_info& q_info);

//power X modulo
void poly_q_powerXmod(poly_q& res,uint64_t power,poly_q& f,poly_q_info& q_info);

//Splitting of a polynomial in a linear terms
void poly_q_split_linear(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info,long term);
 
//This function returns the number of roots,if there is a multiple root,it is returned only once
int poly_q_roots(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info);

//split
int poly_q_split(Vec<GF2X>& roots,poly_q& f,poly_q_info& q_info);

//polynomial evaluation
void poly_q_evaluation(GF2X& y,poly_q& f,GF2X& x,poly_q_info& q_info);

//multiplicity of the root
int poly_q_root_multiplicity(poly_q& pol_p,GF2X& root,poly_q_info& q_info);  

/*---------------------------------------END------------------------------------------------*/


#endif /* POLY_Q_H */

