/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   gray.h
 * Author: pankaj
 *
 * Created on 20 February, 2016, 8:41 PM
 */

#ifndef GRAY_H
#define GRAY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "metaInfo.h"

// Binary Gray code.

# define __GRAY_1            0
# define __GRAY_2  __GRAY_1, 1, __GRAY_1
# define __GRAY_3  __GRAY_2, 2, __GRAY_2
# define __GRAY_4  __GRAY_3, 3, __GRAY_3
# define __GRAY_5  __GRAY_4, 4, __GRAY_4
# define __GRAY_6  __GRAY_5, 5, __GRAY_5
# define __GRAY_7  __GRAY_6, 6, __GRAY_6
# define __GRAY_8  __GRAY_7, 7, __GRAY_7
# define __GRAY_9  __GRAY_8, 8, __GRAY_8
# define __GRAY_10 __GRAY_9, 9, __GRAY_9

  // Length of size-n Gray code is 2^n-1.
# define GRAY_LENGTH(n) ((1u<<(n))-1)

// Ternary Gray code.

#define GRAY(n) CAT(__GRAY_, n)



#ifdef __cplusplus
}
#endif

#endif /* GRAY_H */

