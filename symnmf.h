#include <Python.h>

// symnmf.c functions
double **sym(double **mat, double **symMat, int r, int c);
double **ddg(double **mat, double **diag, double **symMat, int r, int c);
double **norm(double **mat, double **diag, double **symMat, double **NSM, int r, int c);
double **symnmf(double **W, double **H, int r, int k);
void copyMat(double **matrix0, double **matrix1, int n, int m);
double **multiplyMat(double **matrix0, double **matrix1, int n, int m, int t, int isItTranspose);
void getUpdatedH(double **matrix, double **W, double **updatedH, int r, int k);
void printMat(double **matrix, int n, int m);
double frobeniusNorm(double **matrix0, double **matrix1, int n, int m);
double euclidianDist(double *x, double *y, int c);

// symnmfmodule.c functions
static PyObject *fitSym(PyObject *self, PyObject *args);
static PyObject *fitNorm(PyObject *self, PyObject *args);
static PyObject *fitDdg(PyObject *self, PyObject *args);
static PyObject *fitSymnmf(PyObject *self, PyObject *args);
PyMODINIT_FUNC PyInit_geo_capi(void);