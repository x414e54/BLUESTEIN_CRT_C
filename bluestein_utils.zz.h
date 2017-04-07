#include "bluestein_zz.h"

//////////////////////////////
/// Bigint write to file
//////////////////////////////
const size_t NUM_PARAMS = 13;

size_t save_ints(const char* path, zz_t *ints, size_t num, size_t N)
{
    const size_t max_chars = (ZZ_NWORDS * sizeof(zz_t) * 8) + 1;
    char tmp[max_chars];
    
    FILE *f;
    f = fopen(path, "w");
    if (f != NULL) {
        for (size_t i = 0; i < num; ++i) {
            for (size_t j = 0; j < N; ++j) {
                size_t sz = zz_toa(&ints[(i * N) + j], tmp, max_chars);
                fputs(&tmp[sz], f);
                if (j != N - 1) {
                    fputs(",", f);
                }
            }
            fputs("\n", f);
        }
        fclose(f);
    } else {
        printf("File not found: %s\n", path); return 1;
    }
    
    return 0;
}

void build_path(const char* base, const char* type, size_t n, size_t r,
                char * path, size_t max)
{
    snprintf(path, max, "%s/%s.%"PRIuPTR".%"PRIuPTR, base, type, n, r);
}

size_t load_ints(const char* name, zz_t *polys, size_t num, size_t N)
{
    //a_0, a_1, [....], a_(N-1)\n
    const size_t MAX = (N * (ZZ_NWORDS * ceil(Z_NUM_BITS / 3.0) + 2)) - 1 + 1;
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
                    polys[(N * i) + j] = zz_t_zero;
                } else {
                    polys[(N * i) + j] = zz_str(start);
                }
                DEBUG_LOG(zz_print(&polys[(N * i) + j]);)
            }
        }
        fclose(f);
    } else {
        printf("File not found: %s\n", name); return 1;
    }
    free(tmp);
    
    return 0;
}

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
