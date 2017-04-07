#ifndef BLUESTEIN_SUBRING_H
#define BLUESTEIN_SUBRING_H

#include "bluestein_z.h"

///////////////////////////////
/// FFT Parameters 
//////////////////////////////
typedef struct
{
    bred_params p_b; // Barett reduction parameters
    
    z_t Ni; // N*Ni mod p = 1
    
    z_t *y_table;
    z_t *yi_table;
    
    z_t *w_table;
    z_t *wi_table;
    
    z_t *y_fft;
    z_t *yi_fft;
} prime_parameters;

typedef struct
{
    size_t N; // Ring N (Subring g)
    size_t M; // FFT for M > 2N - 1
    size_t k; // M = 2^k
    
    z_t *tmp_fft;
    
    // Could be elsewhere
    size_t num_primes;
    prime_parameters *primes;
    //
    
    size_t *bitrev_table;
    
    prime_parameters *cur_prime;
} parameters;

__BS_INLINE
z_t *table(int64_t k, const parameters *params, size_t tbl_idx)
{
    if (tbl_idx == 0) {
        return &params->cur_prime->y_table[k];
    } else if (tbl_idx == 1) {
        return &params->cur_prime->yi_table[k];
    } else if (tbl_idx == 2) {
        return &params->cur_prime->w_table[k];
    } else if (tbl_idx == 3) {
        return &params->cur_prime->wi_table[k];
    } else {
        printf("Invalid table");
        abort();
    }
}

//////////////////////////////
/// FFT/NTT Butterfly
//////////////////////////////
__BS_INLINE
void radix2dit_ip(z_t *a, z_t *b, const z_t *w, const parameters *params)
{
    z_t tmp = *b;
    printf("[ tmp %llu * %llu ]", tmp, *w);
    z_mulred(&tmp, w, &tmp, &params->cur_prime->p_b);
    printf("[ b %llu - %llu ]", *a, tmp);
    z_subred(a, &tmp, b, &params->cur_prime->p_b.q);
    printf("[ a %llu + %llu ]", *a, tmp);
    z_addred(a, &tmp, a, &params->cur_prime->p_b.q);
}

__BS_INLINE
void radix2dif_ip(z_t *a, z_t *b, const z_t *w, const parameters *params)
{
    z_t tmp = *a;
    z_addred(a, b, a, &params->cur_prime->p_b.q);
    z_subred(b, &tmp, b, &params->cur_prime->p_b.q);
    z_mulred(b, w, b, &params->cur_prime->p_b);
}

/* 
Butterfly implementations currently from
https://www.microsoft.com/en-us/research/wp-content/uploads/2016/05/RLWE-1.pdf
*/

/* Gentleman-Sande butterfly FFT rev input */
__BS_INLINE
void conv_fftdif(z_t *v, const parameters *params2, size_t tbl_idx)
{
    parameters params3 = *params2;
    parameters *params = &params3;
    params->k = 2;
    params->M = 4;
    
    DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL;
    
    size_t d = 1;
    for (size_t s = params->M; s > 1; s >>= 1) {
        size_t start = 0;
        size_t h = s >> 1;
        for (size_t i = 0; i < h; ++i) {
            z_t *w = table(i + h - 1, params, tbl_idx);
            for (size_t j = start; j <= start + d - 1; ++j) {
                printf("j:%ld, d:%ld, s:%ld, i:%ld\n", j, d, s, i);
                radix2dif_ip(&v[j], &v[j + d], w, params);
                DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL;
            }
            start += d << 1;
        }
        d <<= 1;
    }
    
    DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL; abort();
}

/* Cooley-Tukey butterfly FFT non rev input */
__BS_INLINE
void conv_fftdit(z_t *v, const parameters *params, size_t tbl_idx)
{
    //DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL;

    size_t d = params->M;
    for (size_t s = (size_t)1; s < params->M; s <<= 1) {
        d >>= 1;
        for (size_t i = 0; i < s; ++i) {
            z_t *w = table(d * i, params, tbl_idx);
            size_t start = (d * i) << 1;
            for (size_t j = start; j <= start + d - 1; ++j) {
                //printf("j:%ld, d:%ld, s:%ld, i:%ld\n", j, d, s, i);
                radix2dit_ip(&v[j], &v[j + d], w, params);
                //DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL;
            }
        }
    }
    
    //DEBUG_LOG_LOOP_Z(v, params->M) DEBUG_LOG_NL; abort();
}

__BS_INLINE
void bitrev(const z_t *in, z_t *out, const parameters *params)
{
    for (size_t i = 0; i < params->M / 2; ++i) {
        size_t rev = params->bitrev_table[i];
        if (i < rev) {
            z_t tmp = in[rev];
            out[rev] = in[i];
            out[i] = tmp;
        }
    }
}

__BS_INLINE
void conv_copy(const z_t *in, z_t *out, const parameters *params)
{
    for (size_t i = 0; i < (2 * params->N - 1); ++i) {
        size_t idx = abs((int)params->N - 1 - (int)i);
        out[i] = in[idx];
    }
}

__BS_INLINE
void multable_bitrev(const z_t *in, z_t *out, const parameters *params,
                     size_t tbl_idx)
{
    for (size_t i = 0; i < params->N; ++i) {
        size_t rev = params->bitrev_table[i];
        if (i < rev) {
            z_t tmp = in[rev]; // Incase in == out
            z_mulred(&in[i], table(i, params, tbl_idx),
                     &out[rev], &params->cur_prime->p_b);
            z_mulred(&tmp, table(rev, params, tbl_idx),
                     &out[i], &params->cur_prime->p_b);
        }
    }
}

__BS_INLINE
void multable(const z_t *in, z_t *out, const parameters *params,
              size_t tbl_idx)
{
    for (size_t i = 0; i < params->N; ++i) {
        z_mulred(&in[i], table(i, params, tbl_idx),
                 &out[i], &params->cur_prime->p_b);
    }
}

//////////////////////////////
/// FFT Parameters Setup
//////////////////////////////
static
void destroy_prime_param_table(prime_parameters *params)
{
    if (params->y_table != NULL) {
        free(params->y_table);
    }
    params->y_table = NULL;
    
    if (params->yi_table != NULL) {
        free(params->yi_table);
    }
    params->yi_table = NULL;
    
    if (params->w_table != NULL) {
        free(params->w_table);
    }
    params->w_table = NULL;
    
    if (params->wi_table != NULL) {
        free(params->wi_table);
    }
    params->wi_table = NULL;
    
    if (params->y_fft != NULL) {
        free(params->y_fft);
    }
    params->y_fft = NULL;
    
    if (params->yi_fft != NULL) {
        free(params->yi_fft);
    }
    params->yi_fft = NULL;
}

static
void destroy_param_table(parameters *params)
{
    params->N = 0;
    
    params->M = 0;
    params->k = 0;
    
    if (params->bitrev_table != NULL) {
        free(params->bitrev_table);
    }
    params->bitrev_table = NULL;
    
    if (params->tmp_fft != NULL) {
        free(params->tmp_fft);
    }
    params->tmp_fft = NULL;
    
    if (params->primes != NULL) {
        for (size_t i = 0; i < params->num_primes; ++i) {
            destroy_prime_param_table(&params->primes[i]);
        }
        free(params->primes);
    }
    params->primes = NULL;
}

static
void build_prime_param_table(
            size_t i,                   // prime index
            z_t *ni,                    // ni*n mod p = 1
            z_t *w, z_t *wi ,           // w, wi  // (y^2)^m mod p = 1
            z_t *y, z_t *yi,            // y, yi  // (y^2)^n mod p = 1
            bred_params *p_b,
            mred_params *p_m,
            parameters *params)
{
    if (i >= params->num_primes) {
        return;
    }
    
    prime_parameters *prime_params = &params->primes[i];
    destroy_prime_param_table(prime_params);
    
    // Main Bluestein params
    prime_params->Ni = *ni;
    prime_params->p_b = *p_b;
    
    prime_params->y_table = (z_t*)malloc(sizeof(z_t) * params->N);
    prime_params->yi_table = (z_t*)malloc(sizeof(z_t) * params->N);
    //
    
    // Secondary FFT params
    prime_params->w_table = (z_t*)malloc(sizeof(z_t) * params->M / 2);
    prime_params->wi_table = (z_t*)malloc(sizeof(z_t) * params->M / 2);
    //

    // Main Bluestein tables
    prime_params->y_table[0] = z_t_one;
    prime_params->y_table[1] = *y;
    prime_params->yi_table[0] = z_t_one;
    prime_params->yi_table[1] = *yi;
    
    for (uint64_t i = 2; i < params->N; ++i) {
        z_t i2 = (i * i) % (2 * params->N);
        z_powred(&prime_params->y_table[1], &i2,
                 &prime_params->y_table[i], p_m);
        z_powred(&prime_params->yi_table[1], &i2,
                 &prime_params->yi_table[i], p_m);
    }
    DEBUG_LOG(printf("y_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->y_table, params->N) DEBUG_LOG_NL
    DEBUG_LOG(printf("yi_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->yi_table, params->N) DEBUG_LOG_NL
    //
    
    // Secondary FFT tables
    prime_params->w_table[0] = z_t_one;
    prime_params->w_table[1] = *w;
    prime_params->wi_table[0] = z_t_one;
    prime_params->wi_table[1] = *wi;
    
    for (uint64_t i = 2; i < params->M / 2; ++i) {
        z_t _i = i;
        z_powred(&prime_params->w_table[1], &_i,
                 &prime_params->w_table[i], p_m);
        z_powred(&prime_params->wi_table[1], &_i,
                 &prime_params->wi_table[i], p_m);
    }
    //DEBUG_LOG_LOOP_Z(prime_params->w_table, params->M / 2) DEBUG_LOG_NL
    //bitrev(prime_params->w_table, prime_params->w_table, params);
    DEBUG_LOG(printf("w_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->w_table, params->M / 2) DEBUG_LOG_NL
    //bitrev(prime_params->wi_table, prime_params->wi_table, params);
    DEBUG_LOG(printf("wi_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->wi_table, params->M / 2) DEBUG_LOG_NL
    
    prime_params->y_fft = (z_t*)malloc(params->M * sizeof(z_t));
    prime_params->yi_fft = (z_t*)malloc(params->M * sizeof(z_t));
    
    // TODO dont need bitrev?
    params->cur_prime = prime_params;
    conv_copy(prime_params->y_table, prime_params->y_fft, params);
    conv_fftdit(prime_params->y_fft, params, 2);
    DEBUG_LOG(printf("y_fft:");)
    DEBUG_LOG_LOOP_Z(prime_params->y_fft, params->M) DEBUG_LOG_NL
    
    conv_copy(prime_params->yi_table, prime_params->yi_fft, params);
    conv_fftdif(prime_params->yi_fft, params, 2);
    DEBUG_LOG(printf("yi_fft:");)
    DEBUG_LOG_LOOP_Z(prime_params->yi_fft, params->M) DEBUG_LOG_NL
    //
}

static
void build_param_table(size_t n, size_t k, size_t num_primes,
                       parameters *params)
{
    // Main Bluestein params
    params->N = n;
    
    params->primes = (prime_parameters*)calloc(num_primes,
                                               sizeof(prime_parameters));
    params->num_primes = num_primes;
    //
    
    // Secondary FFT params
    params->M = (size_t)1 << k;
    params->k = k;
    
    params->bitrev_table = (size_t*)malloc(sizeof(size_t) * params->M / 2);
    //
    
    // Secondary Global FFT tables
    k -= 1;
    for (size_t i = 0; i < params->M  / 2; ++i) {
        size_t rev = 0;
        for (size_t j = 0; j < k / 2; ++j) {
            rev |= ((i >> j) & 1) << (k - j - 1);
            rev |= ((i >> (k - j - 1)) & 1) << j;
        }
        params->bitrev_table[i] = rev;
    }
    
    params->tmp_fft = (z_t*)malloc(params->M * sizeof(z_t));
    //
}

//////////////////////////////
/// Bluestein FFT 
//////////////////////////////
__BS_INLINE
void mulni(z_t *x, const parameters *params)
{
    for (size_t i = 0 ; i < params->N; ++i) {
        z_mulred(&x[i], &params->cur_prime->Ni, &x[i], &params->cur_prime->p_b);
    }
}

__BS_INLINE
void muly(const z_t *v1, z_t *v2, const parameters *params)
{
    for (size_t i = 0; i < params->M; ++i) {
        z_mulred(&v2[i], &v1[i], &v2[i], &params->cur_prime->p_b);
    }
}

__BS_INLINE
void bluestein(const z_t *in, z_t *out, const parameters *params,
               size_t tbl_idx, size_t inv_tbl_idx, const z_t *y_fft)
{
    // radix 2 butterfly dit fft
    multable(in, params->tmp_fft, params, inv_tbl_idx);
    conv_fftdit(params->tmp_fft, params, 2);
    muly(y_fft, params->tmp_fft, params);
    conv_fftdif(params->tmp_fft, params, 3);
    multable(params->tmp_fft, out, params, inv_tbl_idx);
}

////////////////////////////////////////////////////////////
/// CRT polynomial FFT code (could be elsewhere)
////////////////////////////////////////////////////////////
__BS_INLINE
void pointwise_mul(const z_t *v1, z_t *v2, parameters *params)
{
    for (size_t i = 0; i < params->num_primes; ++i) {
        params->cur_prime = &params->primes[i];
        for (size_t j = 0; j < params->N; ++j) {
            z_mulred(&v2[j], &v1[j], &v2[j], &params->cur_prime->p_b);
        }
        
        v1 += params->N;
        v2 += params->N;
    }
}

__BS_INLINE
void ifft(const z_t *p1, z_t *p2, parameters *params)
{
    for (size_t i = 0; i < params->num_primes; ++i) {
        params->cur_prime = &params->primes[i];

// Inverse Bluestein
        bluestein(p1, p2, params, 1, 0, params->cur_prime->yi_fft);
    
        DEBUG_LOG_LOOP_Z(p2, params->N)
    
        mulni(p2, params);
    
        DEBUG_LOG_LOOP_Z(p2, params->N)
//
        
        p1 += params->N;
        p2 += params->N;
    }
}

__BS_INLINE
void fft(const z_t *p1, z_t *p2, parameters *params)
{
    for (size_t i = 0; i < params->num_primes; ++i) {
        params->cur_prime = &params->primes[i];
        
// Foward Bluestein
        bluestein(p1, p2, params, 0, 1, params->cur_prime->y_fft);
     
        DEBUG_LOG_LOOP_Z(p2, params->N)
//
        
        p1 += params->N;
        p2 += params->N;
    }
}
#endif
