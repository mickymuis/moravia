#include "matrix.h"
#include <stdio.h>

int
main( int argc, const char** argv ) {

    if( argc != 2 ) {
        fprintf( stderr, "(i) Usage: %s [matrix market file]\n", argv[0] );
        return -1;
    }
    matrix_t mat;

    if( !load_matrix_market( argv[1], &mat ) ) {
        fprintf( stderr, "(e) Could not load matrix market file `%s'\n", argv[1] );
        return -1;
    }
    
    dump_nonzeros( &mat );

    matrix_free( &mat );

    return 0;
}
