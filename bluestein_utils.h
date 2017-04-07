const size_t NUM_PARAMS = 13;

__BS_INLINE
size_t load_ints(const char* name, z_t *polys, size_t num, size_t N)
{
    //a_0, a_1, [....], a_(N-1)\n
    const size_t MAX = (N * (ceil(Z_NUM_BITS / 3.0) + 2)) - 1 + 1;
    DEBUG_LOG(printf("Max poly chars per line (digit, ) ~%" PRIuPTR, MAX))
    
    char *tmp = (char *)malloc(MAX);
    if (tmp == NULL) { printf("Unable to allocate\n"); return 1; }
    
    FILE *f;
    f = fopen(name,"r");
    if (f != NULL) {
        for (size_t i = 0; i < num; ++i) {
            if (fgets(tmp, MAX, f) == NULL) {
                printf("Test file does not contain enough ints/polys\n");
                fclose(f);
                free(tmp);
                return 1;
            }
            char *end = tmp;
            for (size_t j = 0; j < N; ++j) {
                char *start = end;
                // Scan till end of number
                while (*end != 0 && *end != ',' && *end != '\n') { ++end; }
                
                // Set end as null char and goto next char
                if (*end != 0) { *end = 0; ++end; }
                
                // Scan till start of number
                while (*start != 0 && *start == ' ') { ++start; }
                
                // Zero pad if there are not enough coefficients
                if (*start == 0) {
                    polys[(N * i) + j] = z_t_zero;
                } else {
                    polys[(N * i) + j] = z_str(start);
                }
                DEBUG_LOG(printf("%" PRIu64, polys[(N * i) + j]);)
            }
        }
        fclose(f);
    } else {
        printf("File not found: %s\n", name); return 1;
    }
    free(tmp);
    
    return 0;
}

__BS_INLINE
size_t load_params(const char *file, size_t n, size_t r,
                   parameters *fft_params, crt_parameters *crt)
{
    // FFT primes p_i as a CRT split of large enough to fit polynomial
    // multiplication mod 2^(r * 2).
    // Such that (p_i - 1) % n = 0 or gcd(p_i, n) = 1 and (p_i - 1) % m = 0
    // params generated from gen_bs_params.c

    // Num primes = ceil(4*r/(k-4))
    size_t num_primes = ceil((float)(4.0 * r) / ((float)Z_NUM_BITS - 4));
    size_t single_num_primes = ceil((float)(2.0 * r) / ((float)Z_NUM_BITS - 4));
    
    build_crt_params(num_primes, 2, crt);
    
    build_crt_prod_params(single_num_primes, &crt->prods[0]); // Single
    build_crt_prod_params(num_primes, &crt->prods[1]); // Double
    
    z_t tmp_m = (2 * n) - 1;
    size_t k = z_np2(&tmp_m);
    
    build_param_table(
        n, k,               // n, k
        fft_params
    );

    z_t ints[NUM_PARAMS * num_primes];
    if (load_ints(file, ints, NUM_PARAMS * num_primes, 1) != 0) { return 1; }
    
    for (size_t i = 0; i < num_primes; ++i) {
        size_t offset = (NUM_PARAMS * i) + 2;
        
        if (i < single_num_primes) {
            crt->prods[0].Pi[i] = ints[offset - 2];
        }
        crt->prods[1].Pi[i] = ints[offset - 1];
        
        bred_params b_params;
        b_params.q = ints[offset + 0];
        b_params.r_shift = z_np2(&b_params.q) - 1;
        b_params.r = ints[offset + 1];
    
        mred_params m_params;
        m_params.q = b_params.q;
        m_params.qi = ints[offset + 2];
        m_params.r_shift = z_np2(&m_params.q);
        m_params.r_mask = (z_t_one << (m_params.r_shift % Z_NUM_BITS)) - 1;
        m_params.r2 = ints[offset + 3];
        m_params.ri = ints[offset + 4];
        
        build_prime_param_table(
            &ints[offset + 5], &ints[offset + 6],   // ni*n mod p = 1 // mi*m mod p = 1
            &ints[offset + 7], &ints[offset + 8],   // w, wi  // (y^2)^m mod p = 1
            &ints[offset + 9], &ints[offset + 10],  // y, yi  // (y^2)^n mod p = 1
            &b_params,
            &m_params,
            fft_params, crt->primes + i
        );
    }
    
    return 0;
}

__BS_INLINE
void build_path(const char* base, const char* type, size_t n, size_t r,
                char * path, size_t max)
{
    snprintf(path, max, "%s/%s.%" PRIuPTR".%" PRIuPTR, base, type, n, r);
}

__BS_INLINE
size_t next_path(char *path, size_t max)
{
    size_t idx = strlen(path) - 1;
    
    while (path[idx] != '_' && path[idx--] != '.');

    size_t i = 0;
    if (path[idx] == '_') {
        sscanf(path + idx, "_%" PRIuPTR, &i);
        ++i;
        if (i > 20) { return 1; }
    } else {
        while (path[idx] != 0) { idx++; }
    }

    snprintf(path + idx, max - idx, "_%" PRIuPTR, i);
    return 0;
}
