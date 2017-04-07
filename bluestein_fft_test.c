#include "bluestein_fft.h"
#include "bluestein_crt.h"
#include "bluestein_utils.h"

//////////////////////////////
/// Test
//////////////////////////////
ssize_t polycmp(const z_t *p1, const z_t *p2, size_t d)
{
    ssize_t same = 0;
    for (size_t i = 0; i < d; ++i) {
        if (p1[i] != p2[i]) {
            same = -1;
        }
    }
    return same;
}

void test_polys(z_t *polys, parameters* fft_params, crt_parameters* crt_params,
                size_t idx)
{
    /* FFT */
    z_t *A_fft = malloc(fft_params->N * sizeof(z_t)
                        * crt_params->prods[idx].length);
    z_t *B_fft = malloc(fft_params->N * sizeof(z_t)
                        * crt_params->prods[idx].length);
    {
        fft(polys, A_fft, fft_params, crt_params, idx);
        fft(polys + (fft_params->N * crt_params->prods[idx].length),
            B_fft, fft_params, crt_params, idx);
    }
    
    pointwise_mul(B_fft, A_fft, fft_params, crt_params, idx);
    
    /* Inverse FFT */
    z_t *AB = malloc(fft_params->N * sizeof(z_t)
                     * crt_params->prods[idx].length);
    {
        ifft(A_fft, AB, fft_params, crt_params, idx);
        
        /* Compare */
        assert(polycmp(polys + (2 * fft_params->N
                                * crt_params->prods[idx].length),
               AB, fft_params->N * crt_params->prods[idx].length) == 0);
    }
    
    /* Cleanup */
    free(A_fft);
    free(B_fft);
    free(AB);
}

int main()
{
    typedef struct {
        size_t n;
        size_t r;
    } test_case;
    
    const size_t PATH_MAX = 256;
    const char* base = "./data.bs";
    
    parameters params;
    crt_parameters crt_params;
    char path[PATH_MAX];
    
    test_case tests[] = {{18, 162}, {7710, 242}}; // TODO get from args or read folder
    
    for (size_t i = 0; i < sizeof(tests) / sizeof(test_case); ++i) {
        build_path(base, "params", tests[i].n, tests[i].r, path, PATH_MAX);
        if (load_params(path, tests[i].n, tests[i].r,
                        &params, &crt_params) != 0) {
            printf("Test skipped"); continue;
        }
        
        z_t polys[3 * tests[i].n * crt_params.num_primes];
        build_path(base, "polys.crt", tests[i].n, tests[i].r, path, PATH_MAX);
        // Currently only testing single precision params
        while (next_path(path, PATH_MAX) == 0 &&
               load_ints(path, polys,
                         3 * crt_params.prods[0].length, params.N) == 0) {
            test_polys(polys, &params, &crt_params, 0);
            printf("Test completed successfully\n");
        }
        
        printf("Tests for this parameter set completed\n");
        
        destroy_param_table(&params);
        destroy_crt_params(&crt_params);
    }
    
    printf("Tests completed\n");
}
