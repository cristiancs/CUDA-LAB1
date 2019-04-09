
#include<iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;


void sumar(float* r1, float* g1, float* b1, float* r, float* g, float* b) {
	//cout << *r << "|" << 1-*r << endl;
	*r1 += *r;
	*g1 += *g;
	*b1 += *b;
}

void promedio(float* r, float* g, float* b, int l) {
	*r = *r/l;
	*g = *g/l;
	*b = *b/l;
}


int main(int argc, char const *argv[]) {
	FILE *pFile;
	int n, m, l;
	float *r, *g, *b; 
	float momentaneo;
	pFile = fopen ("images1.txt","r");
	fscanf(pFile, "%d %d %d", &l, &m, &n);

	r = new float[n*m*l];
	g = new float[n*m*l];
	b = new float[n*m*l];
	for (int j = 0; j < l; ++j){
		for (int i = 0; i < n*m; ++i) {
			fscanf (pFile, "%f", &r[i+(j*n*m)]);
		}

		for (int i = 0; i < n*m; ++i) {
			fscanf (pFile, "%f", &g[i+(j*n*m)]);
		}

		for (int i = 0; i < n*m; ++i) {
			fscanf (pFile, "%f", &b[i+(j*n*m)]);
		}
	}

	fclose (pFile);

	clock_t t1, t2;
	double ms;
	

	t1 = clock();
    for (int i = n*m; i < (n*m*l); ++i) {
		sumar(&r[i%(n*m)],&g[i%(n*m)],&b[i%(n*m)], &r[i], &g[i], &b[i]);
	}

	for (int i = 0; i < n*m; ++i){
		promedio(&r[i], &g[i], &b[i], l);
	}

    t2 = clock();
    ms = 1000.0 * (double)(t2 - t1) / CLOCKS_PER_SEC;
    cout << "Tiempo CPU: " << ms << " [ms]" << endl;

	FILE * pSalida;
	pSalida = fopen ("img_salida.txt","w");
	fprintf(pSalida, "%d %d\n", m, n);
	for (int i = 0; i < n*m; ++i) {
		if(i == n*m - 1) {
			fprintf(pSalida, "%f", r[i]);
		} else {
			fprintf(pSalida, "%f ", r[i]);
		}
		
	}
	fprintf(pSalida, "\n");
	for (int i = 0; i < n*m; ++i) {
		if(i == n*m - 1) {
			fprintf(pSalida, "%f", g[i]);
		} else {
			fprintf(pSalida, "%f ", g[i]);
		}
	}
	fprintf(pSalida, "\n");
	for (int i = 0; i < n*m; ++i) {
		if(i == n*m - 1) {
			fprintf(pSalida, "%f", b[i]);
		} else {
			fprintf(pSalida, "%f ", b[i]);
		}
	}
	delete r;
	delete g;
	delete b;

	//cin.get();
	return 0;
}