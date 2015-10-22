#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "cuda_runtime.h"
#include <time.h>

#define NONE -1
#define CLOSED -2

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void shortest_path_initialize(int32_t N,
                              int32_t N_block,
                              int32_t numOfNodes,
                              int32_t** parents,
                              float **d_minValues,
                              int32_t **d_parents,
                              int32_t **d_heap,
                              int32_t **d_qpos) {
  printf("begin: shortest_path_initialize\n");
  cudaMallocHost(parents, N*numOfNodes *sizeof(int32_t));

  printf("bytes allocated: %lu \n", 4*N_block*numOfNodes*4);

  cudaMalloc(d_minValues, N_block*numOfNodes*sizeof(float));
  cudaMalloc(d_parents, N_block*numOfNodes*sizeof(int32_t));
  cudaMalloc(d_heap, N_block*numOfNodes*sizeof(int32_t));
  cudaMalloc(d_qpos, N_block*numOfNodes*sizeof(int32_t));

  printf("end: shortest_path_initialize\n");
}

__global__ void init_data_kernel(int32_t length, float *d_minValues,
                              int32_t *d_parents, int32_t *d_heap, int32_t *d_qpos) {

  int tid = threadIdx.x + (blockDim.x * ((gridDim.x * blockIdx.y) + blockIdx.x));

  if(tid < length) {
    d_minValues[tid] = FLT_MAX; 
    d_parents[tid] = NONE;
    d_heap[tid] = 0;
    d_qpos[tid] = NONE;
  }

}

void init_data_gpu(int32_t length, float *d_minValues,
                   int32_t *d_parents, int32_t *d_heap, int32_t *d_qpos) {
  int block_size = 1024;
  int size_y = 10;
  dim3 grid_size((int)ceil((float)length/block_size/size_y), size_y, 1);

//  printf("length: %d, block_size: %d, grid_size %d\n", length, block_size, grid_size);

  init_data_kernel<<<grid_size, block_size>>>(length, d_minValues, d_parents, d_heap, d_qpos);

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

  d_minValues[aStartNode] = 0.0;
  d_parents[aStartNode] = NONE;
  d_heap[nhp++] = aStartNode;
  current = NONE;
  while (nhp > 0) {
    current = d_heap[0];
    d_qpos[current] = CLOSED;
    if (nhp > 1) {
      hp1 = d_heap[--nhp];
      nextValue = d_minValues[hp1];
      k = 0;
      k2 = k * 2 + 1;
      while (k2 < nhp) {
        hp2 = d_heap[k2];
        if (k2 < nhp - 1) {
          hp3 = d_heap[k2 + 1];
          if (d_minValues[hp3] < d_minValues[hp2]) {
            hp2 = hp3;
            k2++;
          }
        }
        if (nextValue > d_minValues[hp2]) {
          d_heap[k] = hp2;
          d_qpos[hp2] = k;
          k = k2;
          k2 = k * 2 + 1;
        } else
          break;
      }
      d_heap[k] = hp1;
      d_qpos[hp1] = k;
    } else
      nhp--;

    for (j = d_aNodes[current]; j < d_aNodes[current + 1]; ++j) {
      next = d_bNodes[j];
      if (d_qpos[next] == CLOSED)
        continue;
      nextValue = d_minValues[current] + d_impedances[j];
      if (nextValue < d_minValues[next]) {
        d_minValues[next] = nextValue;
        d_parents[next] = current;
        if (d_qpos[next] == NONE)
          d_qpos[next] = nhp++;
        k = d_qpos[next];
        while (k > 0) {
          k2 = (k - 1) / 2;
          hp2 = d_heap[k2];
          if (nextValue < d_minValues[hp2]) {
            d_heap[k] = hp2;
            d_qpos[hp2] = k;
            k = k2;
          } else
            break;
        }
        d_heap[k] = next;
        d_qpos[next] = k;
      }
    }
  }
}

__global__ void shortest_path_kernel(int32_t N_block,
                                     int32_t numOfNodes,
                                     int32_t start_offset,
                                     int32_t *gd_parents,
                                     int32_t *gd_heap,
                                     int32_t *gd_qpos,
                                     float *gd_minValues,
                                     int32_t *d_aNodes,
                                     int32_t *d_bNodes,
                                     float *d_impedances) {
  int tid = blockIdx.x *blockDim.x + threadIdx.x;

  if(tid < N_block) {
    int aCarID = tid;
    int aStartNode = tid + start_offset;
    int32_t *d_parents = gd_parents + aCarID*numOfNodes;
    int32_t *d_heap = gd_heap + aCarID*numOfNodes;
    int32_t *d_qpos = gd_qpos + aCarID*numOfNodes;
    float *d_minValues = gd_minValues + aCarID*numOfNodes;

    shortest_path_execute(aCarID, aStartNode, numOfNodes, d_parents,
                                                    d_heap,
                                                    d_qpos,
                                                    d_minValues,
                                                    d_aNodes,
                                                    d_bNodes,
                                                    d_impedances);
  }
}

void shortest_path_gpu(int32_t N_block,
                          int32_t numOfNodes,
                          int32_t start_offset,
                          int32_t *d_parents,
                          int32_t *parents,
                          int32_t *d_heap,
                          int32_t *d_qpos,
                          float *d_minValues,
                          int32_t *d_aNodes,
                          int32_t *d_bNodes,
                          float *d_impedances) {

  int block_size = 8;
  int grid_size = (int)ceil((float)N_block/block_size);

  shortest_path_kernel<<<grid_size, block_size>>>(N_block,numOfNodes,start_offset,
                                                  d_parents,
                                                  d_heap,
                                                  d_qpos,
                                                  d_minValues,
                                                  d_aNodes,
                                                  d_bNodes,
                                                  d_impedances);

    cudaMemcpy(parents, d_parents, N_block*numOfNodes*sizeof(int32_t), cudaMemcpyDeviceToHost);
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

  printf("end:read forward star data...\n");
}

void copy_to_device(int32_t *numOfNodes,
                                     int32_t *numOfLinks,
                                     int32_t **aNodes,
                                     int32_t **bNodes,
                                     float **impedances,
                                     int32_t **d_aNodes,
                                     int32_t **d_bNodes,
                                     float **d_impedances) {
  cudaMalloc(d_aNodes, (*numOfNodes+1)*sizeof(int32_t));
  cudaMemcpy(*d_aNodes, *aNodes, (*numOfNodes+1)*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_bNodes, *numOfLinks*sizeof(int32_t));
  cudaMemcpy(*d_bNodes, *bNodes, *numOfLinks*sizeof(int32_t), cudaMemcpyHostToDevice );
  cudaMalloc(d_impedances, *numOfLinks*sizeof(int32_t));
  cudaMemcpy(*d_impedances, *impedances, *numOfLinks*sizeof(float), cudaMemcpyHostToDevice );

}

int main(int argc, char **argv) {

  int32_t numOfNodes;
  int32_t numOfLinks;
  int32_t *aNodes = NULL;
  int32_t *bNodes = NULL;
  float *impedances = NULL;

  int32_t *d_aNodes = NULL;
  int32_t *d_bNodes = NULL;
  float *d_impedances = NULL;

  int32_t *parents = NULL;
  int32_t *d_parents = NULL;
  int32_t *d_heap = NULL;
  int32_t *d_qpos = NULL;
  float *d_minValues = NULL;

  read_forward_star_from_bin_file(&numOfNodes,
                                  &numOfLinks,
                                  &aNodes,
                                  &bNodes,
                                  &impedances,
                                  &d_aNodes,
                                  &d_bNodes,
                                  &d_impedances);

  if(argc < 2) {
    printf("Enter number of cars and block size");
    return 1;
  }

  int N = atoi(argv[1]);
  int N_block = atoi(argv[2]);

   printf("num cars: %d, block size: %d\n",N, N_block);

   clock_t start = clock(), diff;

   shortest_path_initialize(N,
                             N_block,
                           numOfNodes,
                           &parents,
                           &d_minValues,
                           &d_parents,
                           &d_heap,
                           &d_qpos);


  copy_to_device(&numOfNodes,
                                  &numOfLinks,
                                  &aNodes,
                                  &bNodes,
                                  &impedances,
                                  &d_aNodes,
                                  &d_bNodes,
                                  &d_impedances);

  for(int block=0; block< ceil(N/(float)N_block); block++) {
    int32_t num_block = min(N_block, N - block*N_block);
    int32_t start_offset = N_block*block;

    init_data_gpu(num_block*numOfNodes, d_minValues, d_parents, d_heap, d_qpos);
//    gpuErrchk( cudaPeekAtLastError() );

    shortest_path_gpu(num_block,
                     numOfNodes,
                     start_offset,
                      d_parents,
                      parents + N_block*block*numOfNodes,
                      d_heap,
                      d_qpos,
                      d_minValues,
                      d_aNodes,
                      d_bNodes,
                      d_impedances);

//  printf("parents[4000]: %d\n", parents[4000]);
  }

  diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Time taken to compute and copy %d seconds %d milliseconds\n", msec/1000, msec%1000);
  
  printf("parents[4000]: %d\n", parents[4000]);

  return 0;
}
