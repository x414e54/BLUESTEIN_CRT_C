#ifndef BLUESTEIN_Z_H
#define BLUESTEIN_Z_H

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

///////////////////////////////
/// Debug macros 
//////////////////////////////
#if DEBUG
#define DEBUG_LOG(X) X; DEBUG_LOG_NL
#define DEBUG_LOG_LAST(X, i, n) if (i < n - 1) { X; }
#define DEBUG_LOG_NL printf("\n");
#define DEBUG_LOG_LOOP(x, n, p) for (size_t __dd = 0; __dd < n; ++__dd) {\
    p; DEBUG_LOG_LAST(printf(", ");, __dd, n) }; DEBUG_LOG_NL
#else
#define DEBUG_LOG(X)
#define DEBUG_LOG_LAST(X, i, n)
#define DEBUG_LOG_NL
#define DEBUG_LOG_LOOP(x, n, t)
#endif

#define DEBUG_LOG_LOOP_Z(x, n) DEBUG_LOG_LOOP(x, n, printf("%" PRIu64, x[__dd]))

///////////////////////////////
/// "int" types
//////////////////////////////
#define __BS_INLINE static inline
typedef uint64_t z_t; // NWORD e.g. dword
typedef __uint128_t z2_t; // NWORD e.g. dword

const size_t Z_NUM_BITS = 64;
const z_t Z_LSB_MASK = 1;
const z_t Z_MSB_MASK = (z_t)1 << (Z_NUM_BITS - 1);

const z_t z_t_zero = 0;
const z_t z_t_one = 1;

typedef struct
{
    z_t q;
    z_t r;
    size_t r_shift;
} bred_params;

typedef struct
{
    z_t q;
    z_t qi; // q*q = -1 mod r=2^k
    z_t r2;
    z_t ri;
    size_t r_shift;
    z_t r_mask;
} mred_params;

//////////////////////////////
/// Reduction Math functions NP
//////////////////////////////
__BS_INLINE
void z_addred(const z_t *a, const z_t *b, z_t *c, const z_t *q)
{
    *c = *a + *b;
    
    if (*q <= *c) {
        *c = *c - *q;
    }
}

__BS_INLINE
void z_subred(const z_t *a, const z_t *b, z_t *c, const z_t *q)
{
    if (*a < *b) {
        *c = (*q - *b) + *a;
    } else {
        *c = *a - *b;
    }
}

__BS_INLINE
void z_bred(const z2_t *a, z_t *c, const bred_params *params)
{
    // Barret reduction
    z2_t t = *a;
    t >>= params->r_shift;
    t *= params->r;
    t >>= params->r_shift;
    t = t * params->q;
    *c = *a - t;

    if (params->q <= *c) {
        *c = *c - params->q;
    }
}

__BS_INLINE
void z_mred(const z2_t *a, z_t *c, const mred_params *params)
{
    // Montgomery reduction
    *c = ((z_t)*a * params->qi) & params->r_mask;
    
    *c = (*a + (z2_t)*c * (z2_t)params->q) >> params->r_shift;
    
    if (params->q < *c) {
        *c = *c - params->q;
    }
}

__BS_INLINE
void z_mulred(const z_t *a, const z_t *b, z_t *c,
              const bred_params *params)
{
    z2_t tmp = (z2_t)*a;
    tmp *= *b;
    z_bred(&tmp, c, params);
}

__BS_INLINE
void z_mulred2(const z_t *a, const z_t *b, z_t *c,
               const mred_params *params)
{
    z2_t tmp = (z2_t)*b;
    tmp *= params->r2;
    z_t tmp2;
    z_mred(&tmp, &tmp2, params);
    tmp = (z2_t)*a;
    tmp *= tmp2;
    z_mred(&tmp, c, params);
}

__BS_INLINE
void z_powred(const z_t *a, const z_t *b, z_t *c, const mred_params *params)
{
    size_t bits = Z_NUM_BITS;
    
    // Find MSB
    z_t mask = Z_MSB_MASK;
    for (size_t i = 0; i < Z_NUM_BITS; ++i) {
        if (*b & mask) { break; }
        mask >>= 1;
        --bits;
    }
    
    z_t s = z_t_one;
    z_t t = *b;
    *c = *a;
    
    if (bits == 1) {
        return;
    }
    
    if (bits == 0) {
        *c = s;
        return;
    }
 
    mask = Z_LSB_MASK;
    for (size_t i = 0; i < bits - 1; ++i) {
        if ((t & mask) == 0) {
            z_mulred2(c, c, c, params);
        } else {
            z_mulred2(c, &s, &s, params);
            z_mulred2(c, c, c, params);
        }
        
        mask <<= 1;
    }
    
    z_mulred2(c, &s, c, params);
}

__BS_INLINE
size_t z_np2(const z_t *a)
{
    z_t tmp = *a;
    ssize_t last_1 = 0;
    
    for (ssize_t i = 0; i < Z_NUM_BITS; ++i) {
        if (tmp & Z_MSB_MASK) {
            last_1 = i;
            break;
        }
        tmp <<= 1;
    }
    
    last_1 = (Z_NUM_BITS) - 1 - last_1;
    
    ++last_1;
    
    return last_1;
}

__BS_INLINE
z_t z_str(const char* str)
{
    const z_t base = 10;
    
    z_t tmp = 1;
    z_t z = 0;
    
    size_t len = strlen(str);
    for (ssize_t i = len - 1; i >= 0; --i) {
        z_t tmp2 = str[i] - '0';
        tmp2 = tmp * tmp2;
        z = z + tmp2;
        tmp = tmp * base;
    }
    return z;
}
#endif
