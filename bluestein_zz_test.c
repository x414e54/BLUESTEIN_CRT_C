#include "timer.h"

#define ZZ_MAX_NWORDS 10
#include "bluestein_zz.h"

int main()
{
    {
        zz_t zero = zz_uint64(0);
        DEBUG_LOG(zz_print(&zero));
        zz_t one = zz_uint64(1);
        DEBUG_LOG(zz_print(&one));
        
        zz_t q = zz_str("5846006549323611672814739330865132078623730171904");
        zz_t qi = zz_t_one;
        DEBUG_LOG(zz_printb(&q))
        DEBUG_LOG(zz_print(&q))
        
        DEBUG_LOG(printf("Sub-test"))
        zz_t minus = zz_t_zero;
        zz_sub(&minus, &one, &minus);
        DEBUG_LOG(zz_printb(&minus))
        
        DEBUG_LOG(printf("Add-test"))
        zz_t add;
        zz_add(&minus, &minus, &add);
        DEBUG_LOG(zz_printb(&add))
        
        DEBUG_LOG(printf("Mul-test"))
        zz_t mul;
        zz_mul(&minus, &minus, &mul);
        DEBUG_LOG(zz_printb(&mul))
        
        DEBUG_LOG(printf("Div-test"))
        zz_t rem = zz_t_zero;
        zz_t d = zz_uint64(10);
        zz_t qt;
        zz_div_rem(&q, &d, &qt, &rem);
        DEBUG_LOG(zz_printb(&qt))
        DEBUG_LOG(zz_print(&qt))
        DEBUG_LOG(zz_print(&rem))
        
        rem = zz_t_zero;
        d = zz_uint64(12315263);
        zz_div_rem(&q, &d, &qt, &rem);
        DEBUG_LOG(zz_printb(&qt))
        DEBUG_LOG(zz_print(&qt))
        DEBUG_LOG(zz_print(&rem))
        
        DEBUG_LOG(printf("Add-red-test small"))
        bred_params params;
        params.q = zz_str("127");
        DEBUG_LOG(zz_print(&params.q))
        params.r_shift = 8;
        params.r = zz_str("516");
        DEBUG_LOG(zz_print(&params.r))
        zz_t a = zz_uint64(120);
        zz_t b = zz_uint64(80);
        DEBUG_LOG(zz_print(&a))
        DEBUG_LOG(zz_print(&b))
        zz_addred(&a, &b, &add, &params.q);
        DEBUG_LOG(zz_print(&add)) 
        
        DEBUG_LOG(printf("Add-red-test large"))
        params.q = zz_str("6864797660130609714981900799081393217269435300143305"
                          "4093944634591855431833976560521225596406614545549772"
                          "9631139148085803712198799971664381257402829111505715"
                          "1");
        DEBUG_LOG(zz_print(&params.q))
        params.r_shift = 1042;
        params.r = zz_str("6864797660130609714981900799081393217269435300143305"
                          "4093944634591855431833976560521225596406614545549772"
                          "9631139148085803712198799971664381257402829111505715"
                          "3");
        DEBUG_LOG(zz_print(&params.r))
        a = zz_str("68647976601306097149819007990813932172694353001433054093944"
                   "63459185543183397656052122559640661454554977296311391480858"
                   "037121987999716643812574028291115055889");
        b = zz_str("68647976601306097149819007990813932172694353001433054093944"
                   "63459185543183397656052122559640661454554977296311391480858"
                   "037121987999716643812574028291115056811");
        DEBUG_LOG(zz_print(&a))
        DEBUG_LOG(zz_print(&b))
        zz_addred(&a, &b, &add, &params.q);
        
        /* 686479766013060971498190079908139321726943530014330540939446345918554
           318339765605212255964066145455497729631139148085803712198799971664381
           2574028291115055549*/
        
        DEBUG_LOG(zz_print(&add))
        
        DEBUG_LOG(printf("Mul-test"))
        b = zz_str("29261548461302889838782688711419692380627239651802068370568"
                   "5592329307783526392215844889498905607517009");
        zz_mul(&params.q, &b, &mul);
        DEBUG_LOG(zz_print(&mul))
        b = zz_str("20087460940895052122097302704039145229845989819056718811231"
                   "68264577292018875343486642277013790596455087124928929262638"
                   "18736815870545389921173230988934429265956947171170509594020"
                   "16848816668571554507153727657392759040184978585872359605068"
                   "95267326520450239581359");
        DEBUG_LOG(zz_print(&b))
        
        DEBUG_LOG(printf("Powred-test"))
        mred_params params2;
        params2.q = params.q;
        params2.qi = zz_str("1");
        params2.r_shift = 521;
        params2.r_idx = ZZ_MAX_NWORDS - 1 - (params2.r_shift / Z_NUM_BITS);
        params2.r_mask = ((z_t)1 << (params2.r_shift % Z_NUM_BITS)) - 1;
        params2.r2 = zz_str("1");
        //params2.r = zz_str("6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057152");
        params2.ri = zz_str("1");
        zz_t x = zz_uint64(10);
        zz_t pow;
        zz_powred(&x, &x, &pow, &params2);
        DEBUG_LOG(zz_print(&pow))
        x = zz_str("540939446345918554318339765605212255964066145455497");
        zz_t y = zz_str("6");
        zz_mulred(&x, &x, &pow, &params);
        DEBUG_LOG(zz_print(&pow))
        zz_powred(&x, &y, &pow, &params2);
        DEBUG_LOG(zz_print(&pow))
        zz_powred(&x, &x, &pow, &params2);
        DEBUG_LOG(zz_print(&pow))
        zz_powred(&pow, &pow, &pow, &params2);
        DEBUG_LOG(zz_print(&pow))
        x = zz_str("30836930154696239327703235465270660545807034070493490334445"
                   "24938043738981036623154754158014413980206698113937156552031"
                   "240595342613675854473494562025624744146");
        y = zz_str("2");
        zz_powred(&x, &y, &pow, &params2);
        DEBUG_LOG(zz_print(&pow))
        
        DEBUG_LOG(printf("Pow mulred-test"))
        x = zz_str("3700648187917522963662432569620338463673921942095410804685581680264981047797315044018355331808868197855302649864208491985624725732995545526307274016399168855");
        zz_mulred(&x, &x, &pow, &params);
        DEBUG_LOG(zz_print(&pow))
        
        DEBUG_LOG(printf("Shift-test 1"))
        zz_t shift = minus;
        zz_t shift2 = minus;
        DEBUG_LOG(zz_printb(&shift))
        for (int i = 0; i < ZZ_MAX_NWORDS * Z_NUM_BITS; ++i) {
            zz_lshift(&shift, 1);
            printf("Single: \n"); DEBUG_LOG(zz_printb(&shift));
            shift2 = minus;
            zz_lshift(&shift2, i + 1);
            printf("Full: \n"); DEBUG_LOG(zz_printb(&shift2))
        }
        
        shift = minus;
        shift2 = minus;
        for (int i = 0; i < ZZ_MAX_NWORDS * Z_NUM_BITS; ++i) {
            zz_rshift(&shift, 1);
            printf("Single: \n"); DEBUG_LOG(zz_printb(&shift));
            shift2 = minus;
            zz_rshift(&shift2, i + 1);
            printf("Full: \n"); DEBUG_LOG(zz_printb(&shift2))
        }
        
        DEBUG_LOG(printf("Shift Test 2"))
        shift = zz_t_one;
        nano_time_t t1 = getCurrentTime();
        zz_lshift(&shift, (ZZ_MAX_NWORDS * Z_NUM_BITS) - 1);
        nano_time_t t2 = getCurrentTime();
        DEBUG_LOG(zz_printb(&shift))
        nano_time_t diff = t2 - t1;
        printf("Time - %"PRIu64"\n", diff);
        
    }
}
