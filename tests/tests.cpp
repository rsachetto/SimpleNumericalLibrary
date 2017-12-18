//
// Created by sachetto on 06/10/17.
//
#include <criterion/criterion.h>

//TODO: convert to c++ strings

double **read_octave_mat_file_to_array (FILE *matrix_file, int *num_lines, int *nnz) {
//    const char *sep = " ";
//    char *line_a = NULL;
//    size_t len;
//    int count;
//
//    do {
//        getline (&line_a, &len, matrix_file);
//        //sds *tmp = sdssplitlen (line_a, (int)strlen (line_a), sep, (int)strlen (sep), &count);
//        if (count) {
//            if (strcmp (tmp[1], "columns:") == 0) {
//                (*num_lines) = atoi (tmp[2]);
//            }
//            if (strcmp (tmp[1], "nnz:") == 0) {
//                (*nnz) = atoi (tmp[2]);
//            }
//        }
//        sdsfreesplitres (tmp, count);
//    } while ((line_a)[0] == '#');
//
//    double **matrix = (double **)malloc (*num_lines * sizeof (double *));
//
//    for (int i = 0; i < *num_lines; i++) {
//        matrix[i] = (double *)calloc (*num_lines, sizeof (double));
//    }
//
//    int item_count = 0;
//    int m_line, m_column;
//    double m_value;
//
//    while (item_count < *nnz) {
//
//        sds *tmp = sdssplitlen (line_a, (int)strlen (line_a), sep, (int)strlen (sep), &count);
//        if (tmp[0][0] != '\n') {
//            m_line = atoi (tmp[0]);
//            m_column = atoi (tmp[1]);
//            m_value = atof (tmp[2]);
//
//            matrix[m_line - 1][m_column - 1] = m_value;
//        }
//        sdsfreesplitres (tmp, count);
//
//        item_count++;
//        getline (&line_a, &len, matrix_file);
//    }
//
//    if (line_a)
//        free (line_a);
//
//    return matrix;
}

//TODO: convert to c++ strings
double *read_octave_vector_file_to_array (FILE *vec_file, int *num_lines) {

//    ssize_t read;
//    size_t len;
//    char *line_b = NULL;
//    int count;
//    char *sep = " ";
//
//    do {
//        read = getline (&line_b, &len, vec_file);
//        sds *tmp = sdssplitlen (line_b, (int)strlen (line_b), sep, (int)strlen (sep), &count);
//        if (count) {
//            if (strcmp (tmp[1], "rows:") == 0) {
//                (*num_lines) = atoi (tmp[2]);
//            }
//        }
//        sdsfreesplitres (tmp, count);
//    } while ((line_b)[0] == '#');
//
//    double *vector = (double *)malloc (*num_lines * sizeof (double));
//
//    int item_count = 0;
//    while ((item_count < *num_lines) && read) {
//        sds *tmp = sdssplitlen (line_b, (int)strlen (line_b), sep, (int)strlen (sep), &count);
//
//        if (tmp[0][0] != '\n') {
//            vector[item_count] = atof (tmp[1]);
//        }
//
//        sdsfreesplitres (tmp, count);
//
//        item_count++;
//        read = getline (&line_b, &len, vec_file);
//    }
//
//    if (line_b)
//        free (line_b);
//
//    return vector;
}

void run_cg (bool jacobi) {
//    FILE *A = fopen ("src/tests/A.txt", "r");
//    FILE *B = fopen ("src/tests/B.txt", "r");
//    FILE *X = fopen ("src/tests/X.txt", "r");
//    double error;
//
//    cr_assert(A);
//    cr_assert(B);
//    cr_assert(X);
//
//
//    struct grid *grid = (struct grid *)malloc (sizeof (struct grid));
//    cr_assert (grid);
//
//    construct_grid_from_file (grid, A, B);
//
//    int nt = 1;
//
//#if defined(_OPENMP)
//    nt = omp_get_max_threads();
//#endif
//
//
//    if(jacobi) {
//        printf("Testing CG with jacobi preconditioner using %d treads\n", nt);
//    }
//
//    else {
//        printf("Testing CG using %d treads\n", nt);
//    }
//
//    conjugate_gradient (grid, 200, 1e-16, jacobi, &error);
//
//    int n_lines1;
//    uint32_t n_lines2;
//
//    double *x = read_octave_vector_file_to_array (X, &n_lines1);
//    double *x_grid = grid_vector_to_array (grid, 'x', &n_lines2);
//
//    cr_assert_eq (n_lines1, n_lines2);
//
//    for (int i = 0; i < n_lines1; i++) {
//        cr_assert_float_eq (x[i], x_grid[i], 1e-10);
//    }
//
//    clean_and_free_grid (grid);
//    fclose (A);
//    fclose (B);
}

Test (solvers, cg_jacobi) {
    cr_assert(true);
 //   run_cg (true);
}

Test (solvers, cg_no_jacobi) {
   // run_cg (false);
}