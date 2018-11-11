#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <stdlib.h>
#include <stddef.h>
#include <inttypes.h>
#include <stdbool.h>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <sys/types.h>
#endif

typedef int32_t idx_t;   // Type large enough to hold an arbitrary index 

typedef struct {
    idx_t  index;
    uint32_t  nz;
    idx_t* cols;
    float*  values;
} row_t;

typedef struct {
    row_t** rowPtr;
    uint32_t m,n;
    uint32_t nz;
} matrix_t;

void matrix_dump( const matrix_t* mat );

bool matrix_loadMM( const char *filename,  matrix_t *mat );

void matrix_free( matrix_t* mat );

#endif /* __MATRIX_H__ */
