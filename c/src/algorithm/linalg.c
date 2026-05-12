#include "linalg.h"

#include <stdbool.h>
#include <stdio.h>

void mat_print_dim(const mat_t* A) { printf("(%d x %d)\n", A->m, A->n); }

void mat_print(const mat_t* A) {
    for (int j = 0; j < A->m; j++) {
        for (int i = 0; i < A->n; i++) {
            printf("%5.4f\t", *(A->p + A->ld * i + j));
        }
        printf("\n");
    }
    return;
}

static void _mat_dbg(const mat_t* A) {
    const size_t n = 2;
    size_t r = 0;

    while (r < A->m) {
        bool rellipsis = (r >= n) && (r + n < A->m);

        size_t c = 0;
        while (c < A->n) {
            bool cellipsis = (c >= n) && (c + n < A->n);

            printf("    ");
            if (rellipsis && !cellipsis) {
                printf("   :   ");
            } else if (!rellipsis && cellipsis) {
                printf(" ..... ");
            } else if (rellipsis && cellipsis) {
                printf("       ");
            } else {
                printf("%5.4f\t", *(A->p + A->ld * c + r));
            }

            if (cellipsis) {
                c = A->n - n;
            } else {
                ++c;
            }
        }

        printf("\n");

        if (rellipsis) {
            r = A->m - n;
        } else {
            ++r;
        }
    }
}

void mat_dbg(const char* name, const mat_t* A) {
    printf("%s:", name);
    mat_print_dim(A);
    _mat_dbg(A);
}
