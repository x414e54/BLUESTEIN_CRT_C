#include "timer.h"

#include "bluestein_z.h"

int main()
{
    {
        z_t add = z_t_zero;
        z_t mul = z_t_zero;
        
        DEBUG_LOG(printf("Add-red-test small"))
        bred_params params;
        params.q = (127);
        DEBUG_LOG(printf("%" PRIu64, params.q))
        params.r_shift = 7;
        params.r = (129);
        DEBUG_LOG(printf("%" PRIu64, params.r))
        z_t a = (120);
        z_t b = (80);
        DEBUG_LOG(printf("%" PRIu64, a))
        DEBUG_LOG(printf("%" PRIu64, b))
        z_addred(&a, &b, &add, &params.q);
        DEBUG_LOG(printf("%" PRIu64, add))
        
        DEBUG_LOG(printf("Mul-red-test small"))
        z_mulred(&a, &b, &mul, &params);
        DEBUG_LOG(printf("%" PRIu64, mul))
        /*z_mulred3(&a, &b, &mul, &params);
        DEBUG_LOG(printf("%" PRIu64, mul))*/
        
        DEBUG_LOG(printf("Add-red-test large"))
        params.q = (1152921504606875713);
        DEBUG_LOG(printf("%" PRIu64, params.q))
        params.r_shift = 61;
        params.r = (4611686018427272956);
        DEBUG_LOG(printf("%" PRIu64, params.r))
        a = (1152921504606875712);
        b = (1152921504606875712);
        DEBUG_LOG(printf("%" PRIu64, a))
        DEBUG_LOG(printf("%" PRIu64, b))
        z_addred(&a, &b, &add, &params.q);
        DEBUG_LOG(printf("%" PRIu64, add))
        
        DEBUG_LOG(printf("Mul-red-test large"))
        z_mulred(&a, &b, &mul, &params);
        DEBUG_LOG(printf("%" PRIu64, mul))
        /*z_mulred3(&a, &b, &mul, &params);
        DEBUG_LOG(printf("%" PRIu64, mul))*/

        DEBUG_LOG(printf("Powred-test"))
        mred_params params2;
        params2.q = params.q;
        params2.qi = (1947693923657605183);
        params2.r_shift = 61;
        params2.r_mask = ((z_t)1 << (params2.r_shift % (2*Z_NUM_BITS))) - 1;
        params2.r2 = (3303260676);
        params2.ri = (973846961828826865);
        z_t pow;
        z_t x = (1152921504606875712);
        z_t y = (1152921504606933187);
        z2_t xy = (z2_t)x*(z2_t)y;
        z_mred(&xy, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        x = (10);
        z_powred(&x, &x, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        x = (974839422359456298);
        y = (6);
        z_mulred(&x, &x, &pow, &params);
        DEBUG_LOG(printf("%" PRIu64, pow))
        z_powred(&x, &y, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        z_powred(&x, &x, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        z_powred(&pow, &pow, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        x = (1074853835389015833);
        y = (2);
        z_powred(&x, &y, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        z_mulred2(&x, &x, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        
        DEBUG_LOG(printf("Pow mulred-test"))
        x = (10748538353);
        z_mulred(&x, &x, &pow, &params);
        DEBUG_LOG(printf("%" PRIu64, pow))
        z_mulred2(&x, &x, &pow, &params2);
        DEBUG_LOG(printf("%" PRIu64, pow))
        /*z_mulred3(&x, &x, &pow, &params);
        DEBUG_LOG(printf("%" PRIu64, pow))*/
        
        DEBUG_LOG(printf("Add-red-test large 2"))
        params.q = (1152921504606850369);
        DEBUG_LOG(printf("%" PRIu64, params.q))
        params.r_shift = 60;
        params.r = (1152921504606843583);
        DEBUG_LOG(printf("%" PRIu64, params.r))
        a = (596820845941714346);
        b = (1152921504606850368);
        DEBUG_LOG(printf("%" PRIu64, a))
        DEBUG_LOG(printf("%" PRIu64, b))
        z_mulred(&a, &b, &mul, &params);
        DEBUG_LOG(printf("%" PRIu64, mul))
    }
}
