#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "cuda_runtime.h"
#include <time.h>

#define NONE -1
#define CLOSED -2

#define N 3000

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void shortest_path_initialize(int32_t numOfNodes,
                              float **minValues,
                              int32_t** parents,
                              int32_t** heap,
                              int32_t** qpos,
                              float **d_minValues,
                              int32_t **d_parents,
                              int32_t **d_heap,
                              int32_t **d_qpos) {
  printf("begin: shortest_path_initialize\n");
  cudaMallocHost(minValues, N*numOfNodes * sizeof(float));
  cudaMallocHost(parents, N*numOfNodes *sizeof(int32_t));
  cudaMallocHost(heap, N*numOfNodes * sizeof(int32_t));
  cudaMallocHost(qpos, N*numOfNodes *sizeof(int32_t));

   for (int i=0; i < N * numOfNodes; ++i) {
     (*minValues)[i] = FLT_MAX;
    (*parents)[i] = NONE;
    (*qpos)[i] = NONE;
  }

  cudaMalloc(d_minValues, N*numOfNodes*sizeof(float));
  cudaMemcpy(*d_minValues, *minValues, N*numOfNodes*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc(d_parents, N*numOfNodes*sizeof(int32_t));
  cudaMemcpy(*d_parents, *parents, N*numOfNodes*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_heap, N*numOfNodes*sizeof(int32_t));
  cudaMemcpy(*d_heap, *heap, N*numOfNodes*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_qpos, N*numOfNodes*sizeof(int32_t));
  cudaMemcpy(*d_qpos, *qpos, N*numOfNodes*sizeof(int32_t), cudaMemcpyHostToDevice );

  printf("end: shortest_path_initialize\n");
}

__device__ void shortest_path_execute(int aCarID,
                                     int aStartNode,
                                     int32_t numOfNodes,
                                     int32_t *d_parents,
                                     int32_t *d_heap,
                                     int32_t *d_qpos,
                                     float *d_minValues,
                                     int32_t *d_aNodes,
                                     int32_t *d_bNodes,
                                     float *d_impedances) {
  int32_t current;
  int32_t hp1;
  int32_t hp2;
  int32_t hp3;
  int32_t j;
  int32_t k;
  int32_t k2;
  int32_t next;
  float nextValue;
  int32_t nhp = 0;

  d_minValues[aCarID*numOfNodes + aStartNode] = 0.0;
  d_parents[aCarID*numOfNodes + aStartNode] = NONE;
  d_heap[aCarID*numOfNodes + (nhp++)] = aStartNode;
  current = NONE;
  while (nhp > 0) {
    current = d_heap[aCarID*numOfNodes + 0];
    d_qpos[aCarID*numOfNodes + current] = CLOSED;
    if (nhp > 1) {
      hp1 = d_heap[aCarID*numOfNodes + (--nhp)];
      nextValue = d_minValues[aCarID*numOfNodes + hp1];
      k = 0;
      k2 = k * 2 + 1;
      while (k2 < nhp) {
        hp2 = d_heap[aCarID*numOfNodes + k2];
        if (k2 < nhp - 1) {
          hp3 = d_heap[aCarID*numOfNodes + k2 + 1];
          if (d_minValues[aCarID*numOfNodes + hp3] < d_minValues[aCarID*numOfNodes + hp2]) {
            hp2 = hp3;
            k2++;
          }
        }
        if (nextValue > d_minValues[aCarID*numOfNodes + hp2]) {
          d_heap[aCarID*numOfNodes + k] = hp2;
          d_qpos[aCarID*numOfNodes + hp2] = k;
          k = k2;
          k2 = k * 2 + 1;
        } else
          break;
      }
      d_heap[aCarID*numOfNodes + k] = hp1;
      d_qpos[aCarID*numOfNodes + hp1] = k;
    } else
      nhp--;

    for (j = d_aNodes[current]; j < d_aNodes[current + 1]; ++j) {
      next = d_bNodes[j];
      if (d_qpos[aCarID*numOfNodes + next] == CLOSED)
        continue;
      nextValue = d_minValues[aCarID*numOfNodes + current] + d_impedances[j];
      if (nextValue < d_minValues[aCarID*numOfNodes + next]) {
        d_minValues[aCarID*numOfNodes + next] = nextValue;
        d_parents[aCarID*numOfNodes + next] = current;
        if (d_qpos[aCarID*numOfNodes + next] == NONE)
          d_qpos[aCarID*numOfNodes + next] = nhp++;
        k = d_qpos[aCarID*numOfNodes + next];
        while (k > 0) {
          k2 = (k - 1) / 2;
          hp2 = d_heap[aCarID*numOfNodes +k2];
          if (nextValue < d_minValues[aCarID*numOfNodes +hp2]) {
            d_heap[aCarID*numOfNodes +k] = hp2;
            d_qpos[aCarID*numOfNodes +hp2] = k;
            k = k2;
          } else
            break;
        }
        d_heap[aCarID*numOfNodes + k] = next;
        d_qpos[aCarID*numOfNodes + next] = k;
      }
    }
  }
}

__global__ void shortest_path_kernel(int32_t numOfNodes,
                                     int32_t *d_parents,
                                     int32_t *d_heap,
                                     int32_t *d_qpos,
                                     float *d_minValues,
                                     int32_t *d_aNodes,
                                     int32_t *d_bNodes,
                                     float *d_impedances) {
  int tid = blockIdx.x *blockDim.x + threadIdx.x;

  if(tid < N) {
    shortest_path_execute(tid, tid, numOfNodes, d_parents,
                                                    d_heap,
                                                    d_qpos,
                                                    d_minValues,
                                                    d_aNodes,
                                                    d_bNodes,
                                                    d_impedances);
  }
}

void shortest_path_gpu(int32_t numOfNodes,
                          int32_t *d_parents,
                          int32_t *parents,
                          int32_t *d_heap,
                          int32_t *d_qpos,
                          float *d_minValues,
                          int32_t *d_aNodes,
                          int32_t *d_bNodes,
                          float *d_impedances) {

  int block_size = 512;
  int grid_size = (int)ceil((float)N/block_size);

  shortest_path_kernel<<<grid_size, block_size>>>(numOfNodes,
                                                  d_parents,
                                                  d_heap,
                                                  d_qpos,
                                                  d_minValues,
                                                  d_aNodes,
                                                  d_bNodes,
                                                  d_impedances);

//  gpuErrchk( cudaPeekAtLastError() );
  cudaMemcpy(parents, d_parents, N*numOfNodes*sizeof(int32_t), cudaMemcpyDeviceToHost);

//  gpuErrchk( cudaPeekAtLastError() );
}

void read_forward_star_from_bin_file(int32_t *numOfNodes,
                                     int32_t *numOfLinks,
                                     int32_t **aNodes,
                                     int32_t **bNodes,
                                     float **impedances,
                                     int32_t **d_aNodes,
                                     int32_t **d_bNodes,
                                     float **d_impedances
                                    ) {
  printf("begin:read forward star data...\n");
  FILE *ifp;
  if ((ifp = fopen("GDT_100KFS.bin", "r")) == NULL) {
    printf("open file failed\n");
  }
  size_t bytes;
  int64_t numOfNodes64;
  bytes = fread(&numOfNodes64, sizeof(int64_t), 1, ifp);
  *numOfNodes = (int32_t)numOfNodes64;

  int theRecordSize;
  bytes = fread(&theRecordSize, sizeof(int), 1, ifp);

  int64_t *tmp_buff = (int64_t*)malloc((*numOfNodes+1)*sizeof(int64_t));

  cudaMallocHost(aNodes, sizeof(int32_t) * (*numOfNodes + 1));
  bytes = fread(tmp_buff, theRecordSize, (*numOfNodes) + 1, ifp);
  for(int i=0; i<(*numOfNodes+1); i++)
    (*aNodes)[i] = (int32_t)tmp_buff[i];

  int64_t numOfLinks64;
  bytes = fread(&numOfLinks64, sizeof(int64_t), 1, ifp);
  *numOfLinks = (int32_t)numOfLinks64;

  free(tmp_buff);
  tmp_buff = (int64_t*)malloc(*numOfLinks*sizeof(int64_t));

  bytes = fread(&theRecordSize, sizeof(int), 1, ifp);
  cudaMallocHost(bNodes, sizeof(int32_t) * (*numOfLinks));
  bytes = fread(tmp_buff, theRecordSize, *numOfLinks, ifp);
  for(int i=0; i<*numOfLinks; i++)
    (*bNodes)[i] = (int32_t)tmp_buff[i];

  bytes = fread(&theRecordSize, sizeof(int), 1, ifp);
  cudaMallocHost(impedances, sizeof(float) * (*numOfLinks));
  bytes = fread(*impedances, sizeof(float), *numOfLinks, ifp);

  free(tmp_buff);
  fclose(ifp);

  cudaMalloc(d_aNodes, (*numOfNodes+1)*sizeof(int32_t));
  cudaMemcpy(*d_aNodes, *aNodes, (*numOfNodes+1)*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_bNodes, *numOfLinks*sizeof(int32_t));
  cudaMemcpy(*d_bNodes, *bNodes, *numOfLinks*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_impedances, *numOfLinks*sizeof(int32_t));
  cudaMemcpy(*d_impedances, *impedances, *numOfLinks*sizeof(float), cudaMemcpyHostToDevice );

  printf("end:read forward star data...\n");
}

int main(int argc, char **argv) {

  int32_t numOfNodes;
  int32_t numOfLinks;
  int32_t *aNodes = NULL;
  int32_t *bNodes = NULL;
  float *impedances = NULL;

  int32_t *parents = NULL;
  int32_t *heap = NULL;
  int32_t *qpos = NULL;
  float *minValues = NULL;

  int32_t *d_parents = NULL;
  int32_t *d_heap = NULL;
  int32_t *d_qpos = NULL;
  float *d_minValues = NULL;

  int32_t *d_aNodes = NULL;
  int32_t *d_bNodes = NULL;
  float *d_impedances = NULL;

  read_forward_star_from_bin_file(&numOfNodes,
                                  &numOfLinks,
                                  &aNodes,
                                  &bNodes,
                                  &impedances,
                                  &d_aNodes,
                                  &d_bNodes,
                                  &d_impedances);

  clock_t start = clock(), diff;

   shortest_path_initialize(numOfNodes,
                           &minValues,
                           &parents,
                           &heap,
                           &qpos,
                           &d_minValues,
                           &d_parents,
                           &d_heap,
                           &d_qpos);

  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
  start = clock();

  shortest_path_gpu(numOfNodes,
                    d_parents,
                    parents,
                    d_heap,
                    d_qpos,
                    d_minValues,
                    d_aNodes,
                    d_bNodes,
                    d_impedances);

  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);

  cudaMemcpy(minValues, d_minValues, N*numOfNodes*sizeof(float), cudaMemcpyDeviceToHost);

  printf("mineValues[10] = %f parents[10]: %d\n", minValues[10], parents[10]);

  return 0;
}
