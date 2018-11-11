

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "matrix.h"

#include "mmio.h"


/*
 * Load matrix market file
 */

typedef struct {
    size_t row, col;
    float val;
} mmelement_t;

int
compareMMElement( const void* a, const void* b ) {
    mmelement_t* elemA = (mmelement_t*)a;
    mmelement_t* elemB = (mmelement_t*)b;

    int i = elemA->row - elemB->row;

    return i ? i : elemA->col - elemB->col;
}

static bool
read_matrix_market(const char *filename,
        mmelement_t **elements_ptr, uint32_t* nz_ptr,
        uint32_t* n_rows_ptr, uint32_t* n_cols_ptr)
{
    FILE *fh = fopen(filename, "r");
    if (!fh)
    {
        perror("fopen");
        return false;
    }

    printf("Reading matrix '%s'...\n", filename);

    MM_typecode matcode;
    if (mm_read_banner(fh, &matcode) != 0)
    {
        fprintf(stderr, "(e) mm_read_banner failed\n");
        return false;
    }

    int M, N, nz;
    int ret_code;
    ret_code = mm_read_mtx_crd_size(fh, &M, &N, &nz);
    if (ret_code != 0)
    {
        fprintf(stderr, "(e) mm_read_mtx_crd_size failed\n");
        return false;
    }

    *n_rows_ptr = M;
    *n_cols_ptr = N;

    // Reserve space for all elements

    const size_t nmemb = mm_is_symmetric(matcode) ? 2 * nz : nz;
    mmelement_t* elements = calloc( nmemb, sizeof( mmelement_t ) );
    size_t off =0;

    for (int i = 0; i < nz; i++)
    {
        int row, col;
        double val;

        if (mm_is_pattern(matcode))
        {
            fscanf(fh, "%d %d\n", &row, &col);
            val = 1.0;
        }
        else
            fscanf(fh, "%d %d %lg\n", &row, &col, &val);

        row--; /* adjust from 1-based to 0-based */
        col--;

        mmelement_t e = { row, col, val };
        elements[off++] = e;
        if (mm_is_symmetric(matcode) && row != col) {
            e.row = col; e.col =row;
            elements[off++] = e;
        }
    }

    fclose(fh);
    qsort( elements, off, sizeof( mmelement_t ), compareMMElement );

    *elements_ptr = elements;
    *nz_ptr =off;

    return true;
}

/*
 * Transfer matrix elements into Compressed Row Storage structure
 */

static void
load_elements(const mmelement_t* elements, size_t n, size_t nz, row_t* rowPtr[] )
{
    idx_t idx = 0;
    idx_t row = elements[0].row;
    

    idx_t colsTemp[n];
    float  valuesTemp[n];

    for( size_t i =0; i < nz; i++ ) {

        const mmelement_t *it = &elements[i];

        valuesTemp[idx] = it->val;
        colsTemp[idx] = it->col;
        idx++;
        
        // We have complete one row, write it to the datastructure
        if ( i+1 == nz || elements[i+1].row != row)
        {
            rowPtr[row] = malloc( sizeof(row_t) );
            assert( rowPtr[row] != NULL );

            rowPtr[row]->nz = idx;
            rowPtr[row]->values = calloc( idx, sizeof(float) );
            assert( rowPtr[row]->values != NULL );
            rowPtr[row]->cols = calloc( idx, sizeof(idx_t) );
            assert( rowPtr[row]->cols != NULL );

            memcpy( rowPtr[row]->values, valuesTemp, idx * sizeof(float) );
            memcpy( rowPtr[row]->cols, colsTemp, idx * sizeof( idx_t) );
            rowPtr[row]->index =row;

            if( i+1 != nz )
                row = elements[i+1].row;
            idx =0;
        }
    }
}

void
matrix_dump( const matrix_t* mat ) {

    for (int i = 0; i < mat->m; ++i)
    {
        if( mat->rowPtr[i] == NULL ) continue;
        if( mat->rowPtr[i]->index != i ) continue; // contracted node
        for (idx_t idx =0; idx < mat->rowPtr[i]->nz; ++idx)
        {
            printf("(%d)->(%d) weight %f\n", i, mat->rowPtr[i]->cols[idx], mat->rowPtr[i]->values[idx] );
        }
    }
}


bool
matrix_loadMM(const char *filename, matrix_t *mat ) {

    mmelement_t *elements = NULL;

    if (!read_matrix_market(filename, &elements, &mat->nz, &mat->m, &mat->n))
        return false;

    mat->rowPtr = calloc( mat->m, sizeof(row_t*) );

    load_elements(elements, mat->n, mat->nz, mat->rowPtr);
    free( elements );

    printf("(i) Import ok: %d x %d matrix, %d non-zeroes\n",
            mat->m, mat->n, mat->nz);

    return true;
}

void
matrix_free( matrix_t* mat ) {
    for( size_t i =0; i < mat->m; i++ ) {
        if( mat->rowPtr[i] != NULL && mat->rowPtr[i]->index == i ) {
            free( mat->rowPtr[i]->values );
            free( mat->rowPtr[i]->cols );
            free( mat->rowPtr[i] );
        }
    }
    free( mat->rowPtr );
}

