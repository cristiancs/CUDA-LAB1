
#include<iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

__global__ void kernel( float* r_gpu,  float* g_gpu,   float* b_gpu, int N, int n_m) {
	int tId = threadIdx.x + blockIdx.x * blockDim.x;
	int i=n_m;
	while(i < N && tId < n_m) {
			r_gpu[tId] += r_gpu[tId+i];
			g_gpu[tId] += g_gpu[tId+i];
			b_gpu[tId] += b_gpu[tId+i];
			i+=n_m;
	}
}

__global__ void kernel2( float* r_gpu,  float* g_gpu,   float* b_gpu, int N, int l) {
	int tId = threadIdx.x + blockIdx.x * blockDim.x;
	if(tId < N) {
		r_gpu[tId] = r_gpu[tId]/l;
		g_gpu[tId] = g_gpu[tId]/l;
		b_gpu[tId] = b_gpu[tId]/l;
	}
}


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
	pFile = fopen ("images6.txt","r");
	fscanf(pFile, "%d %d %d", &l, &m, &n);

	r = new float[n*m*l];
	g = new float[n*m*l];
	b = new float[n*m*l];

	int block_size = 256;
	int grid_size = (int) ceil((float) n*m / block_size);



	float* r_gpu, *g_gpu, *b_gpu;

	cudaMalloc(&r_gpu, sizeof(float) * n * m * l);
	cudaMalloc(&g_gpu, sizeof(float) * n * m * l);
	cudaMalloc(&b_gpu, sizeof(float) * n * m * l);

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

	cudaMemcpy(r_gpu, r, sizeof(float) * n * m * l, cudaMemcpyHostToDevice);
	cudaMemcpy(g_gpu, g, sizeof(float) * n * m * l, cudaMemcpyHostToDevice);
	cudaMemcpy(b_gpu, b, sizeof(float) * n * m * l, cudaMemcpyHostToDevice);

	

	int tamanio = n * m * l;
	int nm=n * m;

	cudaEvent_t ct1, ct2;
	float dt;
	cudaEventCreate(&ct1);
	cudaEventCreate(&ct2);
	cudaEventRecord(ct1);

	kernel<<<grid_size, block_size>>>(r_gpu, g_gpu, b_gpu, tamanio, nm);
	kernel2<<<grid_size, block_size>>>(r_gpu, g_gpu, b_gpu, nm, l);

	cudaEventRecord(ct2);
	cudaEventSynchronize(ct2);
	cudaEventElapsedTime(&dt, ct1, ct2);

	cout << "Tiempo GPU: " << dt << " [ms]" << endl; 

	cudaMemcpy(r, r_gpu, sizeof(float) * n * m * l, cudaMemcpyDeviceToHost);
	cudaMemcpy(g, g_gpu, sizeof(float) * n * m * l, cudaMemcpyDeviceToHost);
	cudaMemcpy(b, b_gpu, sizeof(float) * n * m * l, cudaMemcpyDeviceToHost);

	FILE * pSalida;
	pSalida = fopen ("gpu_img_salida.txt","w");
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

	cudaFree(r_gpu);
	cudaFree(g_gpu);
	cudaFree(b_gpu);

	//cin.get();
	return 0;
}