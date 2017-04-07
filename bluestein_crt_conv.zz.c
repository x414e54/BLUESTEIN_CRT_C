#include "bluestein_utils.zz.h"

void conv_ints(zz_t *params, zz_t *polys, zz_t *polys_out, size_t num_primes,
               size_t num_polys, size_t N)
{
    for (size_t i = 0; i < num_polys; ++i) {
        for (size_t j = 0; j < num_primes; ++j) {
            zz_t *prime = &params[(j * NUM_PARAMS) + 2];
            for (size_t d = 0; d < N; ++d) {
                zz_t q;
                zz_div_rem(&polys[(i * N) + d],
                           prime, &q,
                           &polys_out[(i * N * num_primes) + (j * N) + d]);
            }
        }
    }
}

int main()
{
    typedef struct {
        size_t n;
        size_t r;
    } test_case;
    
    const size_t PATH_MAX = 256;
    const char* base = "./data.bs";
    
    char path[PATH_MAX];
    char path_out[PATH_MAX];
    
    test_case tests[] = {{18, 162}, {7710, 242}}; // TODO get from args or read folder
    
    for (size_t i = 0; i < sizeof(tests) / sizeof(test_case); ++i) {
        const size_t num_primes_d = ceil((float)(4.0 * tests[i].r)
                                       / ((Z_NUM_BITS / 2.0) - 4.0));
        const size_t num_primes_s = ceil((float)(2.0 * tests[i].r)
                                       / ((Z_NUM_BITS / 2.0) - 4.0));
        
        build_path(base, "params", tests[i].n, tests[i].r, path, PATH_MAX);
        zz_t params[NUM_PARAMS * num_primes_d];
        
        if (load_ints(path, params, NUM_PARAMS * num_primes_d, 1) != 0) {
            printf("Test skipped"); continue;
        }
        
        zz_t *polys = (zz_t*)calloc(3 * tests[i].n, sizeof(zz_t));
        zz_t *polys_out = (zz_t*)calloc(3 * num_primes_d * tests[i].n,
                                        sizeof(zz_t));
                                        
        build_path(base, "polys", tests[i].n, tests[i].r, path, PATH_MAX);
        build_path(base, "polys.crt", tests[i].n,
                   tests[i].r, path_out, PATH_MAX);
        
        while (next_path(path, PATH_MAX) == 0 &&
               next_path(path_out, PATH_MAX) == 0 &&
               load_ints(path, polys, 3, tests[i].n) == 0) {
            // Currently only testing/generating single precision parameters.
            conv_ints(params, polys, polys_out, num_primes_s, 3, tests[i].n);
            save_ints(path_out, polys_out, 3 * num_primes_s, tests[i].n);
        }
        
        free(polys);
        free(polys_out);
    }
}
