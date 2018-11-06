#include <stdlib.h>
#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef size_t indx_t;   // Type large enough to hold an arbitrary index 
typedef ssize_t sindx_t; // Signed type large enough to hold an index

void dump_nonzeros(const size_t     n_rows,
                   const double  values[],
                   const indx_t     col_ind[],
                   const indx_t     row_ptr_begin[],
                   const indx_t     row_ptr_end[]);

bool load_matrix_market(const char *filename,
                        const size_t   max_n_elements,
                        const size_t   max_n_rows,
                        size_t        *nnz,
                        size_t        *n_rows,
                        size_t        *n_cols,
                        double        values[],
                        indx_t        col_ind[],
                        indx_t        row_ptr_begin[],
                        indx_t        row_ptr_end[]);

#endif /* __MATRIX_H__ */
