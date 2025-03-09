#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*Function to copy one matrix to another*/
void copyMat(double **matrix0, double **matrix1, int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            matrix1[i][j] = matrix0[i][j] + 0;
        }
    }
}

/* Function to multiply two matrices (with an option for transpose) */
double **multiplyMat(double **matrix0, double **matrix1, int n, int m, int t, int isItTranspose)
{
    double **matrix;
    int i, j, k;
    double sum = 0;
    matrix = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; i++)
    {
        matrix[i] = (double *)malloc(sizeof(double) * m);
    }

    /* Perform matrix multiplication with or without transposition */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            for (k = 0; k < t; k++)
            {
                if (isItTranspose > 0)
                {
                    sum += matrix0[i][k] * matrix1[j][k];
                }
                else
                {
                    sum += matrix0[i][k] * matrix1[k][j];
                }
            }
            matrix[i][j] = sum;
            sum = 0;
        }
    }
    return matrix;
}

/* Function to update H matrix in NMF */
void getUpdatedH(double **matrix, double **W, double **updatedH, int r, int k)
{
    double **helperMat, **mat1, **mat2;

    int i, j;
    /* Multiply matrices to get necessary components for the update*/
    mat1 = multiplyMat(matrix, matrix, r, r, k, 1);
    helperMat = multiplyMat(mat1, matrix, r, k, r, 0);
    mat2 = multiplyMat(W, matrix, r, k, r, 0);

    /* Update the values of H using the formula */
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < k; j++)
        {
            updatedH[i][j] = matrix[i][j] * (0.5 + 0.5 * (mat2[i][j] / helperMat[i][j]));
        }
    }

    /* Free the allocated memory for the intermediate matrices */
    for (i = 0; i < r; i++)
    {
        free(mat1[i]);
        free(helperMat[i]);
        free(mat2[i]);
    }
    free(helperMat);
    free(mat1);
    free(mat2);
}

/* Function to print a matrix to the console */
void printMat(double **matrix, int n, int m)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%.4f", matrix[i][j]);
            if (j < m - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

/* Function to compute the Frobenius norm between two matrices */
double frobeniusNorm(double **matrix0, double **matrix1, int n, int m)
{
    double sum = 0;
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            sum += pow(matrix0[i][j] - matrix1[i][j], 2);
        }
    }
    return sum;
}

/* Function to compute the Euclidean distance between two vectors */
double euclidianDist(double *x, double *y, int c)
{
    double dist = 0;
    int i;

    for (i = 0; i < c; i++)
    {
        dist += pow((x[i] - y[i]), 2.0);
    }
    return dist;
}

/* Function to compute the symmetric matrix based on the Euclidean distance */
double **sym(double **mat, double **symMat, int r, int c)
{
    int i, j;

    /* Loop through the matrix to calculate the symmetric matrix */
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            /*Set diagonal elements to 0*/
            if (i == j)
            {
                symMat[i][j] = 0;
                continue;
            }
            /* Use Euclidean distance to fill the symmetric matrix */
            symMat[i][j] = exp(-1 * euclidianDist(mat[i], mat[j], c) / 2);
        }
    }
    return symMat;
}

/* Function to compute the degree matrix (DDG) from the symmetric matrix */
double **ddg(double **mat, double **diag, double **symMat, int r, int c)
{
    int i, j;
    double d;

    /* Compute the symmetric matrix first */
    sym(mat, symMat, r, c);

    /* Loop through the rows to compute the degree matrix */
    for (i = 0; i < r; i++)
    {
        d = 0;
        /* Sum the values in each row of the symmetric matrix */
        for (j = 0; j < r; j++)
        {
            d += symMat[i][j];
            diag[i][j] = 0;
        }
        /* Set the diagonal element to the sum of the row */
        diag[i][i] = d;
    }
    return diag;
}

/* Function to normalize the symmetric matrix (NSM) */
double **norm(double **mat, double **diag, double **symMat, double **NSM, int r, int c)
{
    int i, j;

    /* Compute the degree matrix */
    ddg(mat, diag, symMat, r, c);

    /* Normalize the symmetric matrix using the degree matrix */
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            NSM[i][j] = pow(diag[i][i], -0.5) * symMat[i][j] * pow(diag[j][j], -0.5);
        }
    }

    return NSM;
}

/* Function to update H matrix using non-negative matrix factorization (NMF) */
double **symnmf(double **W, double **H, int r, int k)
{
    int i;
    double **updatedH;

    /* Allocate memory for the updated H matrix */
    updatedH = (double **)malloc(sizeof(double *) * r);
    for (i = 0; i < r; i++)
    {
        updatedH[i] = (double *)malloc(k * sizeof(double));
    }

    /* Get the updated H matrix using the current H and W matrices */
    getUpdatedH(H, W, updatedH, r, k);

    i = 0;
    /* Iterate until the Frobenius norm converges */
    while (i < 300 && frobeniusNorm(updatedH, H, r, k) >= 0.0001)
    {
        /* Copy the updated matrix to H and update again */
        copyMat(updatedH, H, r, k);
        getUpdatedH(H, W, updatedH, r, k);
        i++;
    }

    /* Copy the final updated H back */
    copyMat(updatedH, H, r, k);

    /* Free the allocated memory for updatedH */
    for (i = 0; i < r; i++)
    {
        free(updatedH[i]);
    }
    free(updatedH);

    return H;
}

/* Main function that reads the input data, processes the matrix and calls the appropriate function */
int main(int argc, char *argv[])
{
    char *goal, *fileName, *chars, buffer[1024], *token;
    FILE *file;
    int r, c, i, j;
    double val;
    double **X, **diag, **symMat, **NSM;

    /* Check for correct number of arguments */
    if (argc != 3)
    {
        perror("An Error Has Occurred");
        return 1;
    }
    goal = argv[1];
    fileName = argv[2];
    file = fopen(fileName, "r");

    if (file == NULL)
    {
        perror("An Error Has Occurred");
        return 1;
    }
    r = 0;
    c = 0;

    /*Read the matrix dimensions by analyzing the input file*/
    while (fgets(buffer, sizeof(buffer), file))
    {
        if (r == 0)
        {
            chars = buffer;
            while (sscanf(chars, "%lf", &val) == 1)
            {
                c++;
                while (*chars != ',' && *chars != '\t' && *chars != '\0')
                    chars++;
                while (*chars == ',' || *chars == '\t')
                    chars++;
            }
        }
        r++;
    }

    rewind(file);

    X = (double **)malloc(r * sizeof(double *));
    for (i = 0; i < r; i++)
    {
        X[i] = (double *)malloc(c * sizeof(double));

        if (fgets(buffer, sizeof(buffer), file) == NULL)
        {
            perror("An Error Has Occurred");
            fclose(file);
            return 1;
        }

        token = strtok(buffer, ",");
        for (j = 0; j < c; j++)
        {
            if (token == NULL)
            {
                perror("An Error Has Occurred");
                fclose(file);
                return 1;
            }
            X[i][j] = strtod(token, NULL); /* Convert string to double and store in X */
            token = strtok(NULL, ",");
        }
    }
    fclose(file);

    symMat = (double **)malloc(sizeof(double *) * r);
    diag = (double **)malloc(sizeof(double *) * r);
    NSM = (double **)malloc(sizeof(double *) * r);
    for (i = 0; i < r; i++)
    {
        symMat[i] = (double *)malloc(sizeof(double) * r);
        diag[i] = (double *)malloc(sizeof(double) * r);
        NSM[i] = (double *)malloc(sizeof(double) * r);
    }

    if (strcmp(goal, "sym") == 0)
    {
        sym(X, symMat, r, c);
        printMat(symMat, r, r);
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        ddg(X, diag, symMat, r, c);
        printMat(diag, r, r);
    }
    else if (strcmp(goal, "norm") == 0)
    {
        norm(X, diag, symMat, NSM, r, c);
        printMat(NSM, r, r);
    }
    else
    {
        fprintf(stderr, "An Error Has Occurred: Invalid goal '%s'\n", goal);
        return 1;
    }

    for (i = 0; i < r; i++)
    {
        free(symMat[i]);
        free(diag[i]);
        free(NSM[i]);
        free(X[i]);
    }
    free(symMat);
    free(diag);
    free(NSM);
    free(X);

    return 0;
}