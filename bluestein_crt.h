#ifndef BLUESTEIN_CRT_H
#define BLUESTEIN_CRT_H

#include "bluestein_fft.h"

///////////////////////////////
/// CRT Parameters
//////////////////////////////
typedef struct
{
    size_t length;
    z_t *Pi; // P*Pi mod Prod(p_i) = 1
} crt_prod_parameters;

typedef struct
{
    size_t num_primes;
    prime_parameters *primes;
    
    size_t num_prods;
    crt_prod_parameters *prods;
} crt_parameters;

static
void destroy_crt_prod_params(crt_prod_parameters *params)
{
    if (params->Pi != NULL) {
        free(params->Pi);
    }
    params->Pi = NULL;
}

static
void destroy_crt_params(crt_parameters *params)
{
    if (params->primes != NULL) {
        for (size_t i = 0; i < params->num_primes; ++i) {
            destroy_prime_param_table(&params->primes[i]);
        }
        free(params->primes);
    }
    params->primes = NULL;
    
    if (params->prods != NULL) {
        for (size_t i = 0; i < params->num_prods; ++i) {
            destroy_crt_prod_params(&params->prods[i]);
        }
        free(params->prods);
    }
    params->prods = NULL;
    
}

static
void build_crt_prod_params(size_t length, crt_prod_parameters *params)
{
    params->Pi = (z_t*)calloc(length,
                              sizeof(z_t));
                              
    params->length = length;
}

static
void build_crt_params(size_t num_primes, size_t num_prods,
                      crt_parameters *params)
{
    params->primes = (prime_parameters*)calloc(num_primes,
                                               sizeof(prime_parameters));
    params->prods = (crt_prod_parameters*)calloc(num_primes,
                                                 sizeof(crt_prod_parameters));
                              
    params->num_primes = num_primes;
    params->num_prods = num_prods;
}

////////////////////////////////////////////////////////////
/// CRT polynomial FFT code
////////////////////////////////////////////////////////////
__BS_INLINE
void pointwise_mul(const z_t *v1, z_t *v2, const parameters *fft_params,
                   const crt_parameters *crt_params, size_t idx)
{
    DEBUG_LOG(printf("pointwise_mul");)
    for (size_t i = 0; i < crt_params->prods[idx].length; ++i) {
        
        DEBUG_LOG_LOOP_Z(v2, fft_params->N)
        DEBUG_LOG_LOOP_Z(v1, fft_params->N)
        
        for (size_t j = 0; j < fft_params->N; ++j) {
            z_mulred(&v2[j], &v1[j], &v2[j], &crt_params->primes[i].p_b);
        }
        
        DEBUG_LOG_LOOP_Z(v2, fft_params->N)
        
        v1 += fft_params->N;
        v2 += fft_params->N;
    }
}

__BS_INLINE
void ifft(const z_t *p1, z_t *p2, const parameters *fft_params,
          const crt_parameters *crt_params, size_t idx)
{
    DEBUG_LOG(printf("ifft");)
    for (size_t i = 0; i < crt_params->prods[idx].length; ++i) {
    
// Inverse Bluestein
        DEBUG_LOG_LOOP_Z(p1, fft_params->N);
        
        bluestein(p1, p2, fft_params, 0, crt_params->primes[i].yi_fft,
                  &crt_params->primes[i]);
    
        DEBUG_LOG_LOOP_Z(p2, fft_params->N)
    
        mulni(p2, fft_params, &crt_params->primes[i]);
    
        DEBUG_LOG_LOOP_Z(p2, fft_params->N)
//
        
        p1 += fft_params->N;
        p2 += fft_params->N;
    }
}

__BS_INLINE
void fft(const z_t *p1, z_t *p2, const parameters *fft_params,
         const crt_parameters *crt_params, size_t idx)
{
    DEBUG_LOG(printf("fft");)
    for (size_t i = 0; i < crt_params->prods[idx].length; ++i) {
    
// Foward Bluestein
        DEBUG_LOG_LOOP_Z(p1, fft_params->N);
        
        bluestein(p1, p2, fft_params, 1, crt_params->primes[i].y_fft,
                  &crt_params->primes[i]);

        DEBUG_LOG_LOOP_Z(p2, fft_params->N);
//
        
        p1 += fft_params->N;
        p2 += fft_params->N;
    }
}
#endif
