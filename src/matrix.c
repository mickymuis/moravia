

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "matrix.h"

#include "mmio.h"


/*
 * Load matrix market file
 */

typedef struct {
    size_t row, col;
    double val;
} mmelement_t;

int
compareMMElement( void* a, void* b ) {
    mmelement_t* elemA = (mmelement_t*)a;
    mmelement_t* elemB = (mmelement_t*)b;

    int i = a->row - b->row;

    return i ? i : a->col - b->col;
}

struct Element
{
  size_t row, col;
  double val;

  Element(int row, int col, double val) : row(row), col(col), val(val) { }
};

inline bool operator==(const struct Element &a, const struct Element &b)
{
  return a.row == b.row && a.col == b.col && a.val == b.val;
}

inline bool operator!=(const struct Element &a, const struct Element &b)
{
  return !operator==(a, b);
}

inline bool operator<(const struct Element &a, const struct Element &b)
{
  if (a.row < b.row)
    return true;

  if (a.row == b.row && a.col < b.col)
    return true;

  return false;
}



static bool
read_matrix_market(const char *filename,
                   mmelement_t **elements_ptr, size_t* nz_ptr
                   size_t* n_rows_ptr, size_t* n_cols_ptr)
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
      elements[offs++] = e;
      if (mm_is_symmetric(matcode) && row != col) {
          e.row = col; e.col =row;
        elements[offs++] = e;
      }
    }

  fclose(fh);
    qsort( elements, offs, sizeof( mmelement_t ), compareMMElement );

    *elements_ptr = elements;
    *nz_ptr =offs;

  return true;
}

/*
 * Transfer matrix elements into Compressed Row Storage structure
 */

static void
load_elements(const mmelement_t* elements, size_t nz,
              double values[],
              indx_t col_ind[],
              indx_t row_ptr_begin[], indx_t row_ptr_end[])
{
  indx_t current_val = 0;
  indx_t current_row = 0;
  row_ptr_begin[0] = 0;

  for( size_t i =0; i < nz; i++ ) {

    {
        mmelement_t *it = &elements[i];
      if (it->row != current_row)
        {
          if (current_row + 1 != it->row)
            {
              fprintf(stderr, "(e) Row skipping not implemented.\n");
              abort();
            }

          row_ptr_end[current_row] = current_val - 1;
          current_row++;
          row_ptr_begin[current_row] = current_val;
        }

      values[current_val] = it->val;
      col_ind[current_val] = it->col;
      current_val++;
    }

  row_ptr_end[current_row] = current_val - 1;
}

void
dump_nonzeros(const size_t n_rows,
              const double values[],
              const indx_t col_ind[],
              const indx_t row_ptr_begin[],
              const indx_t row_ptr_end[])
{
  for (size_t row = 0; row < n_rows; ++row)
    {
      for (indx_t idx = row_ptr_begin[row]; idx <= row_ptr_end[row]; ++idx)
        {
          printf("%ld %ld %f\n", row, col_ind[idx], values[idx]);
        }
    }
}

bool
load_matrix_market(const char *filename,
                   const size_t max_n_elements,
                   const size_t max_n_rows,
                   size_t nnz, size_t n_rows, size_t n_cols,
                   double values[],
                   indx_t col_ind[],
                   indx_t row_ptr_begin[], indx_t row_ptr_end[])
{

    mmelemet_t *elements = NULL;

  if (!read_matrix_market(filename, &elements, nnz, n_rows, n_cols))
    return false;

  if (nnz >= (size_t)max_n_elements || n_rows >= max_n_rows)
    return false;

  load_elements(elements, *nnz, values, col_ind, row_ptr_begin, row_ptr_end);
    free( *elements );

  printf("(i) Import ok: %ld x %ld matrix, %ld non-zeroes\n",
         *n_rows, *n_cols, *nnz);

  return true;
}
