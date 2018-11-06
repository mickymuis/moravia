#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <sys/types.h>
#endif

typedef size_t indx_t;   // Type large enough to hold an arbitrary index 
typedef ssize_t sindx_t; // Signed type large enough to hold an index

typedef struct {
    double* values;
    indx_t* col_ind;
    indx_t* row_ptr_begin;
    indx_t* row_ptr_end;
    size_t m,n;
    size_t nz;
} matrix_t;

void dump_nonzeros( const matrix_t* mat );

bool load_matrix_market(const char *filename,  matrix_t *mat );

void matrix_free( matrix_t* mat );

#endif /* __MATRIX_H__ */
