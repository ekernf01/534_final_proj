#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

double* allocvector(int n);
void printvector(double* v,int n);
void fprintvector(double* v,int n, char* filename);
void freevector(double* v);
void vectorproduct(int n,double* v1,double* v2,double* v);
void readvector(int n,double* v,char* filename);
