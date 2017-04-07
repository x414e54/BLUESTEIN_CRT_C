#ifndef BLUESTEIN_ZZ_H
#define BLUESTEIN_ZZ_H

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

#define DEBUG_LOG_LOOP_ZZ(x, n) DEBUG_LOG_LOOP(x, n, zz_print(&x[__dd]))

///////////////////////////////
/// "big int" types
//////////////////////////////
#if !defined(ZZ_MAX_NWORDS)
#define ZZ_MAX_NWORDS 10
#endif

#define __ZZ_INLINE static inline
typedef uint64_t _z_t; // NWORD/2 e.g. dword
typedef __uint128_t z_t; // NWORD e.g. qword

const size_t Z_NUM_BITS = 128;
const z_t Z_LSB_MASK = 1;
const z_t Z_MSB_MASK = (z_t)1 << (Z_NUM_BITS - 1);

const size_t Z_2NWORDS = 2;
const size_t ZZ_NWORDS = ZZ_MAX_NWORDS;
const size_t ZZ_2NWORDS = Z_2NWORDS * ZZ_NWORDS;
typedef union
{
    z_t z[ZZ_NWORDS];
    _z_t _z[ZZ_2NWORDS];
} zz_t;

__ZZ_INLINE
zz_t zz_uint64(uint64_t val)
{
    zz_t zz;
    for (size_t i = 0; i < ZZ_NWORDS - 1; ++i) {
        zz.z[i] = 0;
    }
    zz.z[ZZ_NWORDS - 1] = val;
    return zz;
}

const z_t z_t_zero = 0;
const z_t z_t_one = 1;
#define zz_t_zero zz_uint64(0)
#define zz_t_one zz_uint64(1)

// Converts index to correct hi/lo byte order
// Could be lookup based.
__ZZ_INLINE
size_t zz_2norder(size_t i)
{
    size_t x = i % 2;
    size_t z = ZZ_2NWORDS - i - 1;
    size_t y = z % 2;
    z += x;
    z -= y;
    return z;
}

typedef struct
{
    zz_t q;
    zz_t r;
    size_t r_shift;
} bred_params;

typedef struct
{
    zz_t q;
    zz_t qi; // q*q = -1 mod r=2^k
    zz_t r2;
    zz_t ri;
    size_t r_shift;
    z_t r_mask;
    size_t r_idx;
} mred_params;

//////////////////////////////
/// "big int" Math functions NP
//////////////////////////////
// TODO these math functions can be more efficient and more complete.
__ZZ_INLINE
void zz_print(const zz_t *a);
__ZZ_INLINE
void zz_printb(const zz_t *a);

__ZZ_INLINE
ssize_t zz_cmp(const zz_t *a, const zz_t *b)
{
    for (size_t i = 0; i < ZZ_NWORDS; ++i) {
        if (a->z[i] < b->z[i]) { return -1; }
        if (a->z[i] > b->z[i]) { return 1; }
    }
    return 0;
}

__ZZ_INLINE
ssize_t zz_cmp_zero(const zz_t *a)
{
    for (size_t i = 0; i < ZZ_NWORDS; ++i) {
        if (a->z[i] > 0) { return 1; }
    }
    return 0;
}

__ZZ_INLINE
ssize_t zz_cmp_one(const zz_t *a)
{
    for (size_t i = 0; i < ZZ_NWORDS - 1; ++i) {
        if (a->z[i] != 0) { return 1; }
    }
    
    if (a->z[ZZ_NWORDS - 1] > 1) { return 1; }
    if (a->z[ZZ_NWORDS - 1] == 0) { return -1; }
    return 0;
}

__ZZ_INLINE
void zz_lshift(zz_t *a, size_t b)
{
    size_t ob = b / Z_NUM_BITS;
    size_t shift = b % Z_NUM_BITS;
    for (size_t i = 0; i < ZZ_NWORDS; ++i) {
        if ((i + ob) < ZZ_NWORDS) {
            a->z[i] = a->z[i + ob];
        } else {
            a->z[i] = z_t_zero;
        }
        if (shift != 0) {
            a->z[i] <<= shift;
            if ((i + ob + 1) < ZZ_NWORDS) {
                a->z[i] ^= (a->z[i + ob + 1] >> (Z_NUM_BITS - shift));
            }
        }
    }
}

__ZZ_INLINE
void zz_rshift(zz_t *a, size_t b)
{
    ssize_t ob = b / Z_NUM_BITS;
    ssize_t shift = b % Z_NUM_BITS;
    for (ssize_t i = ZZ_NWORDS - 1; i >= 0; --i) {
        if ((i - ob) >= 0) {
            a->z[i] = a->z[i - ob];
        } else {
            a->z[i] = z_t_zero;
        }
        if (shift != 0) {
            a->z[i] >>= shift;
            if ((i - ob - 1) >= 0) {
                a->z[i] ^= (a->z[i - ob - 1] << (Z_NUM_BITS - shift));
            }
        }
    }
}

__ZZ_INLINE
void zz_add(const zz_t *a, const zz_t *b, zz_t *c)
{
    z_t carry = z_t_zero;
    for (ssize_t i = ZZ_NWORDS - 1; i >= 0; --i) {
        z_t tmp = a->z[i] + b->z[i] + carry;
        carry = (tmp < a->z[i] || (carry == 1 && tmp == a->z[i]));
        c->z[i] = tmp;
    }
}

__ZZ_INLINE
void zz_sub(const zz_t *a, const zz_t *b, zz_t *c)
{
    z_t borrow = 0;
    for (ssize_t i = ZZ_NWORDS - 1; i >= 0; --i) {
        z_t tmp = a->z[i] - b->z[i] - borrow;
        borrow = (tmp > a->z[i] || (borrow == 1 && tmp == a->z[i]));
        c->z[i] = tmp;
    }
}

__ZZ_INLINE
void zz_mul(const zz_t *a, const zz_t *b, zz_t *c)
{
    zz_t _c = zz_t_zero;
    
    union __z_t {
        z_t z;
        _z_t _z[Z_2NWORDS];
    };
    
    union __z_t tmp = {z_t_zero};
    zz_t tmpzz = zz_t_zero;
    
    for (size_t i = 0; i < ZZ_2NWORDS; ++i) {
        size_t x = zz_2norder(i);
        for (size_t j = 0; j < ZZ_2NWORDS - i; ++j) {
            size_t y = zz_2norder(j);
            tmp.z = (z_t)a->_z[y] * (z_t)b->_z[x];
            tmpzz = zz_t_zero;
            for (size_t k = 0;
                 k < Z_2NWORDS && i + j + k < ZZ_2NWORDS; ++k) {
                size_t z = zz_2norder(i + j + k);
                tmpzz._z[z] = tmp._z[k];
            }
            zz_add(&_c, &tmpzz, &_c);
        }
    }
    
    *c = _c; // For now incase a==c
}

__ZZ_INLINE
void zz_addred(const zz_t *a, const zz_t *b, zz_t *c, const zz_t *q)
{
    zz_add(a, b, c);
    
    if (zz_cmp(q, c) <= 0) {
        zz_sub(c, q, c);
    }
}

__ZZ_INLINE
void zz_bred(const zz_t *a, zz_t *c, const bred_params *params)
{
    // Barret reduction
    zz_t t;
    zz_mul(a, &params->r, &t);
    zz_rshift(&t, params->r_shift * 2);
    zz_mul(&t, &params->q, &t);
    zz_sub(a, &t, c);

    if (zz_cmp(&params->q, c) <= 0) {
        zz_sub(c, &params->q, c);
    }
}

__ZZ_INLINE
void zz_mred(const zz_t *a, zz_t *c, const mred_params *params)
{
    // Montgomery reduction
    zz_t m; // a == c
    zz_mul(a, &params->qi, &m);
    m.z[params->r_idx] &= params->r_mask;
    for (size_t i = 0; i < params->r_idx; ++i) {
        m.z[i] = z_t_zero;
    }
    
    zz_mul(&m, &params->q, &m);
    zz_add(a, &m, c);
    zz_rshift(c, params->r_shift);

    if (zz_cmp(&params->q, c) <= 0) {
        zz_sub(c, &params->q, c);
    }
}

__ZZ_INLINE
void zz_mulred(const zz_t *a, const zz_t *b, zz_t *c,
                const bred_params *params)
{
    zz_t tmp; // TODO zz_2
    zz_mul(a, b, &tmp); // TODO zz_2
    zz_bred(&tmp, c, params);
}

__ZZ_INLINE
void zz_mulred2(const zz_t *a, const zz_t *b, zz_t *c,
               const mred_params *params)
{
    zz_t tmp; // a == c
    zz_mul(b, &params->r2, &tmp); // TODO zz_2
    zz_mred(&tmp, &tmp, params);
    zz_mul(a, &tmp, c); // TODO zz_2
    zz_mred(c, c, params);
}

__ZZ_INLINE
void zz_div_rem(const zz_t *a, const zz_t *b, zz_t *q, zz_t *rem)
{
    const zz_t one = zz_t_one;
    
    zz_t t = *b;
    zz_t r = *a;
    *q = zz_t_zero;
    
    ssize_t i = 0;
    
    while (zz_cmp_zero(&t) != 0 &&
           zz_cmp(&t, &r) < 1 &&
           ((t.z[0] & Z_MSB_MASK) == 0)) {
        zz_lshift(&t, 1);
        ++i;
    }
    
    for (;i >= 0; --i) {
        zz_lshift(q, 1);
        if (zz_cmp(&r, &t) >= 0) {
            zz_sub(&r, &t, &r);
            zz_add(q, &one, q);
        }
        zz_rshift(&t, 1);
    }
    
    if (rem != NULL) {
        *rem = r;
    }
}

__ZZ_INLINE
void zz_powred(const zz_t *a, const zz_t *b, zz_t *c, const mred_params *params)
{
    size_t bits = (ZZ_NWORDS * Z_NUM_BITS);
    
    // Find MSB
    z_t mask = Z_MSB_MASK;
    for (size_t i = 0; i < ZZ_NWORDS * Z_NUM_BITS; ++i) {
        if (i % Z_NUM_BITS == 0) {
            mask = Z_MSB_MASK;
        }
        if (b->z[i / Z_NUM_BITS] & mask) { break; }
        mask >>= 1;
        --bits;
    }
    
    zz_t s = zz_t_one;
    zz_t t = *b;
    *c = *a;
    
    if (bits == 1) {
        return;
    }
    
    if (bits == 0) {
        *c = s;
        return;
    }
 
    size_t z_i = ZZ_NWORDS;
    for (size_t i = 0; i < bits - 1; ++i) {
        if (i % Z_NUM_BITS == 0) {
            mask = Z_LSB_MASK;
            --z_i;
        }
        
        // Even
        if ((t.z[z_i] & mask) == 0) {
            zz_mulred2(c, c, c, params);
        } else {
            zz_mulred2(c, &s, &s, params);
            zz_mulred2(c, c, c, params);
        }
        
        mask <<= 1;
    }
    
    zz_mulred2(c, &s, c, params);
}

__ZZ_INLINE
size_t zz_np2(const zz_t *a)
{
    zz_t tmp = *a;
    ssize_t last_1 = 0;
    
    for (ssize_t i = 0; i < Z_NUM_BITS * ZZ_MAX_NWORDS; ++i) {
        if (tmp.z[0] & Z_MSB_MASK) {
            last_1 = i;
            break;
        }
        zz_lshift(&tmp, 1);
    }
    
    last_1 = (Z_NUM_BITS * ZZ_MAX_NWORDS) - 1 - last_1;
    
    ++last_1;
    
    return last_1;
}

__ZZ_INLINE
zz_t zz_str(const char* str)
{
    const zz_t base = zz_uint64(10);
    
    zz_t tmp = zz_t_one;
    zz_t zz = zz_t_zero;
    
    size_t len = strlen(str);
    for (ssize_t i = len - 1; i >= 0; --i) {
        zz_t tmp2 = zz_uint64(str[i] - '0');
        zz_mul(&tmp, &tmp2, &tmp2);
        zz_add(&zz, &tmp2, &zz);
        zz_mul(&tmp, &base, &tmp);
    }
    return zz;
}

__ZZ_INLINE
void zz_printb(const zz_t *a)
{
    for (int i = 0; i < ZZ_NWORDS; ++i) {
        z_t tmp = a->z[i];
        for (int j = 0; j < sizeof(z_t) * 8; ++j) {
            printf("%c", '0' + (int)((tmp & Z_MSB_MASK) > 0));
            tmp <<= 1;
        }
    }
}

__ZZ_INLINE
size_t zz_toa(const zz_t *a, char *arr, size_t max_chars)
{
    const zz_t base = zz_uint64(10);
    
    zz_t _a = *a;
    
    size_t i = max_chars - 1;
    arr[i] = 0;
    
    do {
        zz_t rem = zz_t_zero;
        zz_div_rem(&_a, &base, &_a, &rem);
        arr[--i] = '0' + (int)rem.z[ZZ_NWORDS - 1];
    } while (zz_cmp_zero(&_a) > 0);
    
    return i;
}

__ZZ_INLINE
void zz_print(const zz_t *a)
{
    const size_t max_chars = (ZZ_NWORDS * sizeof(z_t) * 8) + 1;
    char tmp[max_chars];
    
    size_t i = zz_toa(a, tmp, max_chars);
    
    printf("%s", &tmp[i]);
}
#endif
