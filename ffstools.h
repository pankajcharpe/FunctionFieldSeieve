/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ffstools.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 11:36 AM
 */

#ifndef FFSTOOLS_H
#define FFSTOOLS_H
#include<NTL/GF2X.h>
#include<NTL/ZZ.h>
#include<NTL/vector.h>
#include<string>
#include <cstdint>
#include<fstream>

using namespace std;
using namespace NTL;

void GF2X_set_ti(GF2X& pol,long i);
void GF2X_submul(GF2X& f,GF2X& g,GF2X& h);
ZZ hex_to_decimal(char hex[]);
ZZ hex_to_binary(char hex[]);
 GF2X set_GF2X(ZZ num);
uint64_t GF2X_to_uint64_t(GF2X& f);
GF2X uint64_t_to_GF2X(uint64_t num);
uint32_t GF2X_to_uint32_t(GF2X& f);
 void uint64_t_swap(uint64_t *x,uint64_t *y);
GF2X uint32_t_to_GF2X(uint32_t num);
uint32_t hex_to_uint32_t(string hex);


#endif /* FFSTOOLS_H */

