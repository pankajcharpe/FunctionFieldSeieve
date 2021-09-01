/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bucketInfo.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 3:02 PM
 */

#ifndef BUCKETINFO_H
#define BUCKETINFO_H


#include <iostream>
#include<NTL/tools.h>
#include<NTL/GF2X.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <limits.h>
#include"ij_vector.h"
using namespace std;
using namespace NTL;

// Size of Gray codes to use for the inner loop.
#define ENUM_LATTICE_UNROLL 8
#define SCALE 0

#define BUCKET_REGION_BITS 16
#define UPDATE_HINT_BITS  0 //for time being
#define UPDATE_POS_BITS  BUCKET_REGION_BITS	
#define UPDATE_BITS     (UPDATE_POS_BITS + UPDATE_HINT_BITS)
#define UPDATE_ALIGN      1



//defined in buckets.h(forward eclaration)
typedef struct update_packed_struct update_packed_t;

typedef struct {
  unsigned          n;
  unsigned          max_size;
  unsigned          min_degp, max_degp;
  update_packed_t **start;
  update_packed_t **degp_end;
} buckets_struct;

typedef buckets_struct  buckets_t;

typedef uint16_t update_t;

struct update_packed_struct
 {
  uint8_t d[((UPDATE_BITS-1) / (8*UPDATE_ALIGN) + 1) * UPDATE_ALIGN];
};


typedef union 
{
  update_t        update;
  update_packed_t packed;
} update_conv_t;

typedef unsigned pos_t;
typedef unsigned hint_t;

static inline update_packed_t update_pack(update_t update);
static inline update_t update_unpack(update_packed_t packed);
static inline update_t update_set_hint(hint_t hint);
static inline update_t update_adj_pos(update_t update, pos_t pos);
static inline update_t update_set(pos_t pos, hint_t hint);
static inline pos_t update_get_pos(update_t update);
static inline hint_t update_get_hint(update_t update);
void buckets_init(buckets_t& buckets, unsigned I, unsigned J,double hits, unsigned min_degp, unsigned max_degp);
void buckets_clear(buckets_t& buckets);
unsigned bucket_region_size();
void print_bucket_info(buckets_t& buckets0,  buckets_t& buckets1);
static inline void buckets_push_update(buckets_t& buckets,update_packed_t **ptr, hint_t hint,
                         ij_vec& v, unsigned I, unsigned J);
void buckets_fill(buckets_t& buckets,  large_factor_base_t& FB,
         sub_lattice *sublat, unsigned I, unsigned J,  q_lattice& q_lat);
void bucket_apply(uint8_t *S,  buckets_t& buckets, unsigned k);

#endif /* BUCKETINFO_H */

