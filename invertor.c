/* A simple example of inverting a matrix using the gsl */


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>


gsl_matrix *invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
void randomize_mat_contents(gsl_matrix *matrix);


static size_t size = 5;


        /************************************************************
         * PROCEDURE: invert_a_matrix
         * 
         * DESCRIPTION: Invert a matrix using GSL.
         *
         * RETURNS:
         *      gsl_matrix pointer
         */
gsl_matrix *
invert_a_matrix(gsl_matrix *matrix)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}


        /************************************************************
         * PROCEDURE: print_mat_contents
         * 
         * DESCRIPTION: Print the contents of a gsl-allocated matrix
         * 
         * RETURNS: 
         *      None.
         */
void
print_mat_contents(gsl_matrix *matrix)
{
    size_t i, j;
    double element;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}


        /************************************************************
         * PROCEDURE: randomize_mat_contents
         * 
         * DESCRIPTION: Overwrite entries in matrix with randomly
         *              generated values.
         *
         * RETURNS:
         *      None.
         */
void
randomize_mat_contents(gsl_matrix *matrix)
{
    size_t i, j;
    double random_value;
    double range = 1.0 * RAND_MAX;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {

            // generate a random value
            random_value = rand() / range;
    
            // set entry at i, j to random_value
            gsl_matrix_set(matrix, i, j, random_value);

        }
    }
}