#ifndef BLUESTEIN_FFT_H
#define BLUESTEIN_FFT_H

#include "bluestein_z.h"

///////////////////////////////
/// FFT Parameters 
//////////////////////////////
typedef struct
{
    bred_params p_b; // Barett reduction parameters
    
    z_t Ni; // N*Ni mod p = 1
    z_t Mi; // M*Mi mod p = 1
    
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
    
    size_t *bitrev_table;
} parameters;

__BS_INLINE
z_t *table(int64_t k, const parameters *params, size_t tbl_idx,
           const prime_parameters *p)
{
    if (tbl_idx == 0) {
        return &p->y_table[k];
    } else if (tbl_idx == 1) {
        return &p->yi_table[k];
    } else if (tbl_idx == 2) {
        return &p->w_table[k];
    } else if (tbl_idx == 3) {
        return &p->wi_table[k];
    } else {
        printf("Invalid table");
        abort();
    }
}

//////////////////////////////
/// FFT/NTT Butterfly
//////////////////////////////
__BS_INLINE
void radix2dit_ip(z_t *a, z_t *b, const z_t *w,
                  const parameters *params, const prime_parameters *p)
{
    z_t tmp = *b;
    z_mulred(&tmp, w, &tmp, &(p->p_b));
    z_subred(a, &tmp, b, &(p->p_b.q));
    z_addred(a, &tmp, a, &(p->p_b.q));
}

/* Cooley-Tukey butterfly FFT rev input */
__BS_INLINE
void conv_fftdit(z_t *v, const parameters *params,
                 size_t tbl_idx, const prime_parameters *p)
{
    size_t d = 1;
    size_t m = params->M;
    for (size_t s = 0; s < params->k; ++s) { // Stage
        m >>= 1;
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < d; ++j) {
                z_t *w = table(m * j, params, tbl_idx, p);
                size_t idx = j + (2 * d * i);
                radix2dit_ip(&v[idx], &v[idx + d], w, params, p);
            }
        }
        d <<= 1;
    }
}

__BS_INLINE
void bitrev(z_t *out, const parameters *params)
{
    for (size_t i = 0; i < params->M; ++i) {
        size_t rev = params->bitrev_table[i];
        if (i < rev) {
            z_t tmp = out[rev];
            out[rev] = out[i];
            out[i] = tmp;
        }
    }
}

__BS_INLINE
void conv_bitrev(const z_t *in, z_t *out, const parameters *params)
{
    for (size_t i = 0; i < params->M; ++i) {
        size_t rev = params->bitrev_table[i];
        if (rev < (2 * params->N - 1)) {
            size_t rev_idx = abs((int)params->N - 1 - (int)rev);
            out[i] = in[rev_idx];
        }
    }
}

__BS_INLINE
void multable(const z_t *in, z_t *out, const parameters *params,
              size_t tbl_idx, const prime_parameters *p)
{
    for (size_t i = 0; i < params->N; ++i) {
        z_mulred(&in[i + params->N - 1], table(i, params, tbl_idx, p),
                 &out[i], &(p->p_b));
        z_mulred(&out[i], &p->Mi,
                 &out[i], &(p->p_b));
    }
}

__BS_INLINE
void multable_bitrev(const z_t *in, z_t *out, const parameters *params,
                     size_t tbl_idx, const prime_parameters *p)
{
    for (size_t i = 0; i < params->M; ++i) {
        size_t rev = params->bitrev_table[i];
        if (i < params->N) {
            z_mulred(&in[i], table(i, params, tbl_idx, p),
                     &out[rev], &(p->p_b));
        } else {
            out[rev] = 0;
        }
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
}

static
void build_prime_param_table(
            z_t *ni, z_t *mi,           // ni*n mod p = 1 // mi*m mod p = 1
            z_t *w, z_t *wi,            // w, wi  // (y^2)^m mod p = 1
            z_t *y, z_t *yi,            // y, yi  // (y^2)^n mod p = 1
            bred_params *p_b,
            mred_params *p_m,
            parameters *params, prime_parameters *prime_params)
{
    destroy_prime_param_table(prime_params);
    
    // Main Bluestein params
    prime_params->Ni = *ni;
    prime_params->p_b = *p_b;
    
    prime_params->y_table = (z_t*)malloc(sizeof(z_t) * params->N);
    prime_params->yi_table = (z_t*)malloc(sizeof(z_t) * params->N);
    //
    
    // Secondary FFT params
    prime_params->Mi = *mi;
    
    prime_params->w_table = (z_t*)malloc(sizeof(z_t) * params->M / 2);
    prime_params->wi_table = (z_t*)malloc(sizeof(z_t) * params->M / 2);
    //

    // Main Bluestein tables
    prime_params->y_table[0] = z_t_one;
    prime_params->y_table[1] = *y;
    prime_params->yi_table[0] = z_t_one;
    prime_params->yi_table[1] = *yi;
    
    // TODO include Ni and Mi in these tables?
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
    
    DEBUG_LOG(printf("w_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->w_table, params->M / 2) DEBUG_LOG_NL
    DEBUG_LOG(printf("wi_table:");)
    DEBUG_LOG_LOOP_Z(prime_params->wi_table, params->M / 2) DEBUG_LOG_NL
    
    prime_params->y_fft = (z_t*)malloc(params->M * sizeof(z_t));
    prime_params->yi_fft = (z_t*)malloc(params->M * sizeof(z_t));
    
    // TODO need bitrev?
    conv_bitrev(prime_params->y_table, prime_params->y_fft, params);
    conv_fftdit(prime_params->y_fft, params, 2, prime_params);
    DEBUG_LOG(printf("y_fft:");)
    DEBUG_LOG_LOOP_Z(prime_params->y_fft, params->M) DEBUG_LOG_NL
    
    conv_bitrev(prime_params->yi_table, prime_params->yi_fft, params);
    conv_fftdit(prime_params->yi_fft, params, 2, prime_params);
    DEBUG_LOG(printf("yi_fft:");)
    DEBUG_LOG_LOOP_Z(prime_params->yi_fft, params->M) DEBUG_LOG_NL
    //
}

static
void build_param_table(size_t n, size_t k, parameters *params)
{
    // Main Bluestein params
    params->N = n;
    //
    
    // Secondary FFT params
    params->M = (size_t)1 << k;
    params->k = k;
    
    params->bitrev_table = (size_t*)malloc(sizeof(size_t) * params->M);
    //
    
    // Secondary Global FFT tables
    for (size_t i = 0; i < params->M; ++i) {
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
void mulni(z_t *x, const parameters *params, const prime_parameters *p)
{
    for (size_t i = 0 ; i < params->N; ++i) {
        z_mulred(&x[i], &p->Ni, &x[i], &(p->p_b));
    }
}

__BS_INLINE
void muly(const z_t *v1, z_t *v2, const parameters *params,
          const prime_parameters *p)
{
    for (size_t i = 0; i < params->M; ++i) {
        z_mulred(&v2[i], &v1[i], &v2[i], &(p->p_b));
    }
}

__BS_INLINE
void bluestein(const z_t *in, z_t *out, const parameters *params,
               size_t inv_tbl_idx, const z_t *y_fft, const prime_parameters *p)
{
    // radix 2 butterfly dit fft
    multable_bitrev(in, params->tmp_fft, params, inv_tbl_idx, p);
    conv_fftdit(params->tmp_fft, params, 2, p);
    muly(y_fft, params->tmp_fft, params, p);
    bitrev(params->tmp_fft, params);
    conv_fftdit(params->tmp_fft, params, 3, p);
    multable(params->tmp_fft, out, params, inv_tbl_idx, p);
}
#endif
