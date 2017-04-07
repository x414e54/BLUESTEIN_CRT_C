#include "bluestein_utils.zz.h"

//////////////////////////////
/// Bluestein FFT Parameter generation
//////////////////////////////

// Hard code for now
#define MAX_PRIMES 11
size_t primes[MAX_PRIMES] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 257};

// Find a better way to do this
zz_t calc_omega(size_t n, bred_params *p_b, mred_params *p_m)
{
    const zz_t one = zz_t_one;
    printf("Omega: n = %"PRIuPTR"\n", n);
    
    // Find an n that divides q - 1 = phi(q)
    zz_t q_minus1;
    zz_sub(&p_m->q, &one, &q_minus1);
    
    zz_t rem;
    zz_t _n = zz_uint64(n);
    zz_div_rem(&q_minus1, &_n, &_n, &rem);
    assert(zz_cmp_zero(&rem) == 0);
    
    zz_t z = zz_t_one;
    
    // Split n into prime factors p_i^{c_i}
    uint64_t tmp = n;
    uint64_t exponent = 0;
    for (size_t i = 0; i < MAX_PRIMES; ++i) {
        exponent = 0;
        
        while (tmp != 1 && tmp % primes[i] == 0) {
            tmp /= primes[i];
            ++exponent;
        }
        
        if (exponent == 0) { continue; }
        
        zz_t x = zz_t_one;
        zz_t y = zz_t_zero;
        
        // Find an x such that x^{(q - 1)/p_i} != 1 mod q
        zz_t p = zz_uint64(primes[i]);
        zz_t q_minus1_div_p_i;
        zz_div_rem(&q_minus1, &p, &q_minus1_div_p_i, NULL);
        do {
            zz_add(&x, &one, &x);
            zz_powred(&x, &q_minus1_div_p_i, &y, p_m);
        } while (zz_cmp_zero(&y) == 0 || zz_cmp_one(&y) == 0);
        
        // z = z * x_i^{(q - 1)/p_i^{c_i}} mod q
        zz_t q_minus1_div_p_i_c_i;
        zz_t exp = zz_uint64(exponent);
        zz_powred(&p, &exp, &q_minus1_div_p_i_c_i, p_m);
        zz_div_rem(&q_minus1, &q_minus1_div_p_i_c_i, &q_minus1_div_p_i_c_i, NULL);
        zz_powred(&x, &q_minus1_div_p_i_c_i, &y, p_m);
        
        zz_mulred(&z, &y, &z, p_b);
    }
    
    return z;
}

zz_t calc_psi(size_t n, const zz_t *w, bred_params *p_b, mred_params *p_m)
{
    printf("Psi: n = %"PRIuPTR"\n", n);
    
    zz_t tmp = calc_omega(2 * n, p_b, p_m);
    zz_t y = zz_t_one;
    zz_t y2 = zz_t_one;
    for (size_t i = 1; i < n; ++i) {
        zz_mulred(&y, &tmp, &y, p_b);
        zz_mulred(&y, &y, &y2, p_b);
        if (zz_cmp(w, &y2) == 0) { break; }
    }
    return tmp;
}

zz_t inverse(const zz_t *a, const zz_t *p)
{
    const zz_t one = zz_t_one;
    
    zz_t max_zz = zz_t_one;
    zz_lshift(&max_zz, (Z_NUM_BITS * ZZ_NWORDS) - 1);
    zz_sub(&max_zz, &one, &max_zz);
    
    zz_t n = zz_t_zero;
    zz_t n2 = zz_t_one;
    zz_t r = *p;
    zz_t r2 = *a;
    
    zz_t q = zz_t_zero;
    zz_t tmp = zz_t_zero;
    zz_t tmp2 = zz_t_zero;
    while (zz_cmp_zero(&r2) != 0) {
        tmp = r2;
        zz_div_rem(&r, &r2, &q, &r2);
        r = tmp;
        tmp = n2;
        zz_mul(&q, &n2, &tmp2);
        zz_sub(&n, &tmp2, &n2);
        n = tmp;
    }
    
    if (zz_cmp(&max_zz, &n) < 0) {
        zz_add(&n, p, &n);
    }
    return n;
}

void calc_barret_params(bred_params *p_b)
{
    printf("==Barret params==\n");
    p_b->r_shift = zz_np2(&p_b->q) - 1;
    printf("%"PRIdPTR"\n" , p_b->r_shift);
    
    p_b->r = zz_t_one;
    zz_lshift(&p_b->r, p_b->r_shift * 2);
    zz_div_rem(&p_b->r, &p_b->q, &p_b->r, NULL);
    zz_print(&p_b->r); printf("\n");
}

void calc_montgomery_params(mred_params *p_m)
{
    printf("==Montgomery params==\n");
    p_m->r_shift = zz_np2(&p_m->q);
    printf("%"PRIdPTR"\n" , p_m->r_shift);
    p_m->r_idx = ZZ_NWORDS - 1 - (p_m->r_shift / Z_NUM_BITS);
    printf("%"PRIdPTR"\n" , p_m->r_idx);
    p_m->r_mask = ((z_t)1 << (p_m->r_shift % Z_NUM_BITS)) - 1;
    zz_t r = zz_t_one;
    zz_lshift(&r, p_m->r_shift);
    zz_print(&r); printf("\n");
    printf("Inverse: r\n");
    p_m->ri = inverse(&r, &p_m->q);
    zz_print(&p_m->ri); printf("\n");
    printf("Square: r\n");
    zz_mul(&r, &r, &p_m->r2);
    zz_t tmp;
    zz_div_rem(&p_m->r2, &p_m->q, &tmp, &p_m->r2);
    zz_print(&p_m->r2); printf("\n");
    printf("Inverse: q\n");
    zz_t one = zz_t_one;
    zz_mul(&p_m->ri, &r, &p_m->qi);
    zz_sub(&p_m->qi, &one, &p_m->qi);
    zz_div_rem(&p_m->qi, &p_m->q, &p_m->qi, NULL);
    zz_print(&p_m->qi); printf("\n");
}

void gen_params_prime(uint64_t g, uint64_t r, uint64_t m, zz_t *ints,
                      const zz_t *q)
{
    printf("==Generating Prime Params==\n");
    
    zz_t _g = zz_uint64(g);
    zz_t _m = zz_uint64(m);
    
    bred_params p_b;
    p_b.q = *q;
    calc_barret_params(&p_b);
        
    mred_params p_m;
    p_m.q = p_b.q;
    calc_montgomery_params(&p_m);
    
    printf("==NTT/FFT params==\n");
    
    printf("Inverse: %"PRIu64" (n)\n", g);
    zz_t ni = inverse(&_g, &p_b.q);
    zz_print(&ni); printf("\n");
    
    zz_t w = calc_omega(g, &p_b, &p_m);
    zz_print(&w); printf("\n");
        
    printf("Inverse: Omega (w)\n");
    zz_t wi = inverse(&w, &p_b.q); // inverse w = w^(g-1)
    zz_print(&wi); printf("\n");
    
    zz_t y = calc_psi(g, &w, &p_b, &p_m);
    zz_print(&y); printf("\n");
        
    printf("Inverse: Psi (y)\n");
    zz_t yi = inverse(&y, &p_b.q); // inverse y = y^(g-1)
    zz_print(&yi); printf("\n");
    
    printf("Padding FFT Params: \n");

    printf("Inverse: %"PRIu64" (m)\n", m);
    zz_t mi = inverse(&_m, &p_b.q);
    zz_print(&mi); printf("\n");
    
    zz_t w2 = calc_omega(m, &p_b, &p_m);
    zz_print(&w2); printf("\n");
        
    printf("Inverse: Omega (w)\n");
    zz_t w2i = inverse(&w2, &p_b.q); // inverse w2 = w2^(g-1)
    zz_print(&w2i); printf("\n");
    
    ints[0] = p_b.q;
    ints[1] = p_b.r;
    ints[2] = p_m.qi;
    ints[3] = p_m.r2;
    ints[4] = p_m.ri;
    ints[5] = ni;
    ints[6] = mi;
    ints[7] = w2;
    ints[8] = w2i;
    ints[9] = y;
    ints[10] = yi;
}

size_t calc_crt_primes(zz_t **primes, zz_t *product,
                       zz_t *double_product, size_t *s, size_t r, size_t k)
{
    printf("==Generating CRT split==\n");
    *product = zz_t_one;
    *double_product = zz_t_one;
    
    // Hard code for now - just a random split
    if (r == 162 && k == 64) {
        // Num primes = ceil(2*r/(k-4))
        *primes = (zz_t*)malloc(11 * sizeof(zz_t));
        (*primes)[0] = zz_uint64(1152921504606850369);
        (*primes)[1] = zz_uint64(1152921504606864193);
        (*primes)[2] = zz_uint64(1152921504606875713);
        (*primes)[3] = zz_uint64(1152921504606913729);
        (*primes)[4] = zz_uint64(1152921504606920641);
        (*primes)[5] = zz_uint64(1152921504606924673);
        for (size_t i = 0; i < 6; ++i) {
            zz_mul(product, &(*primes)[i], product);
        }
        *s = 6;
        // Double precision ceil(4*r/(k-4))
        (*primes)[6] = zz_uint64(1152921504606932161);
        (*primes)[7] = zz_uint64(1152921504606933889);
        (*primes)[8] = zz_uint64(1152921504606940801);
        (*primes)[9] = zz_uint64(1152921504606944833);
        (*primes)[10] = zz_uint64(1152921504606948289);
        for (size_t i = 0; i < 11; ++i) {
            zz_mul(double_product, &(*primes)[i], double_product);
        }
        return 11;
    } else if (r == 242 && k == 64) {
        // Num primes = ceil(2*r/(k-4))
        *primes = (zz_t*)malloc(9 * sizeof(zz_t));
        (*primes)[0] = zz_uint64(1152921504843694081);
        (*primes)[1] = zz_uint64(1152921505159495681);
        (*primes)[2] = zz_uint64(1152921506927984641);
        (*primes)[3] = zz_uint64(1152921507054305281);
        (*primes)[4] = zz_uint64(1152921507117465601);
        (*primes)[5] = zz_uint64(1152921507180625921);
        (*primes)[6] = zz_uint64(1152921508001710081);
        (*primes)[7] = zz_uint64(1152921508380672001);
        (*primes)[8] = zz_uint64(1152921508443832321);
        for (size_t i = 0; i < 9; ++i) {
            zz_mul(product, &(*primes)[i], product);
        }
        *s = 9;
        // Double precision ceil(4*r/(k-4))
        (*primes)[9] = zz_uint64(1152921511033405441);
        (*primes)[10] = zz_uint64(1152921512865054721);
        (*primes)[11] = zz_uint64(1152921513054535681);
        (*primes)[12] = zz_uint64(1152921513117696001);
        (*primes)[13] = zz_uint64(1152921513686138881);
        (*primes)[14] = zz_uint64(1152921515644108801);
        (*primes)[15] = zz_uint64(1152921515896750081);
        (*primes)[16] = zz_uint64(1152921517412597761);
        for (size_t i = 0; i < 17; ++i) {
            zz_mul(double_product, &(*primes)[i], double_product);
        }
        return 17;
    } else {
        printf("%ld", k);
        abort();
    }
    return 0;
}

void calc_crt_prime_inverse(const zz_t *q, const zz_t *prod_q, zz_t *qi)
{
    printf("Inverse: (Prod(p_i)/q)\n");
    zz_t tmp;
    zz_div_rem(prod_q, q, &tmp, NULL);
    *qi = inverse(&tmp, q); // inverse (Prod(p_i)/q)
    zz_print(qi); printf("\n");
}

void gen_params(size_t g, size_t r, size_t k)
{
    size_t m = (2 * g) - 1; // 2^np2(2*g - 1)
    zz_t _m = zz_uint64(m);
    m = (size_t)1 << zz_np2(&_m);
    
    printf("==Generating Bluestein CRT Params==\n");
    zz_t *small_primes;
    zz_t prod_small_primes;
    zz_t double_prod_small_primes;
    size_t num_single_primes = 0;
    size_t num_primes = calc_crt_primes(&small_primes,
                                        &prod_small_primes,
                                        &double_prod_small_primes,
                                        &num_single_primes, r, k);
    zz_t ints[NUM_PARAMS * num_primes];
    
    printf("Product of CRT Primes: \n");
    zz_print(&prod_small_primes); printf("\n");
    printf("Double Product of CRT Primes: \n");
    zz_print(&double_prod_small_primes); printf("\n");
    
    for (size_t i = 0; i < num_primes; ++i) {
        printf("\n==CRT Prime==\n");
        zz_print(&small_primes[i]); printf("\n");
        if (i < num_single_primes) {
            calc_crt_prime_inverse(&small_primes[i], &prod_small_primes,
                                   ints + (i * NUM_PARAMS));
        } else {
            ints[(i * NUM_PARAMS)] = zz_t_zero;
        }
        calc_crt_prime_inverse(&small_primes[i], &double_prod_small_primes,
                               ints + (i * NUM_PARAMS) + 1);
        gen_params_prime(g, r, m, ints + (i * NUM_PARAMS) + 2,
                         &small_primes[i]);
    }
    
    printf("\n==Save==\n");
    const size_t PATH_MAX = 256;
    const char* base = "./data.bs";
    
    char path[PATH_MAX];
    build_path(base, "params", g, r, path, PATH_MAX);
    save_ints(path, ints, NUM_PARAMS * num_primes, 1);
    
    free(small_primes);
}

int main()
{
    size_t k = Z_NUM_BITS / 2; // WORD bit length
    // Parameters for r=162 and g=18 subring
    {
        size_t r = 162;
        size_t g = 18;
        gen_params(g, r, k);
    }
    // Parameters for r=242 and g=7710 subring
    {
        size_t r = 242;
        size_t g = 7710;
        gen_params(g, r, k);
    }
}
