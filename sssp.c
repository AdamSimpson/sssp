// Implimentation of "A New GPU-based Approach to the Shortest Path Problem"

#include <stdio.h>
#include <stdlib.h>
//#include <stdbool.h>
#include <float.h>
#include "safe_alloc.h"
#include <math.h>

typedef int bool;
#define true  1
#define false 0

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct GRAPH {
  int64_t arc_count;
  int64_t node_count;

  // Using forward star representation
  // edge_index[i] is first index in edge_tail and weights for i'th node
  // edge_index[i+1] is the last index + 1 in edge_tail and weights for i'th node
  int64_t *restrict edge_index;
  int64_t *restrict edge_tail;
  float   *restrict weights;
};

struct SSSP {
  bool    *restrict unsettled;
  bool    *restrict frontier;
  float   *restrict distance;
  float   *restrict delta_node;
  int64_t *restrict parents; 
  int64_t start_node;
};

void Initialize(struct SSSP *sssp, struct GRAPH *graph, float *delta_dist) {
  sssp->unsettled  = SAFE_ALLOC(graph->node_count, sizeof(bool));
  sssp->frontier   = SAFE_ALLOC(graph->node_count, sizeof(bool));
  sssp->distance   = SAFE_ALLOC(graph->node_count, sizeof(float));
  sssp->delta_node = SAFE_ALLOC(graph->node_count, sizeof(float));
  sssp->parents    = SAFE_ALLOC(graph->node_count, sizeof(int64_t));

  *delta_dist = 0.0f;
  sssp->start_node = 0;

  for(int64_t i=0; i<graph->node_count; i++) {
    sssp->distance[i]  = FLT_MAX;
    sssp->frontier[i]  = false;
    sssp->unsettled[i] = true;
    sssp->parents[i]   = -1;
  }

  sssp->distance[sssp->start_node]  = 0.0f;
  sssp->frontier[sssp->start_node]  = true;
  sssp->unsettled[sssp->start_node] = false;
  sssp->parents[sssp->start_node]   = sssp->start_node;

  #pragma acc enter data copyin(sssp[:1],                            \
                                sssp->unsettled[0:graph->node_count],  \
                                sssp->frontier[0:graph->node_count],   \
                                sssp->distance[0:graph->node_count],   \
                                sssp->delta_node[0:graph->node_count], \
                                sssp->parents[0:graph->node_count])
}

void Relaxation(struct SSSP *sssp, struct GRAPH *graph) {

  bool* frontier = sssp->frontier;
  int64_t *edge_index = graph->edge_index;
  int64_t *edge_tail = graph->edge_tail;
  float *distance = sssp->distance;
  int64_t *parents = sssp->parents;
  bool *unsettled = sssp->unsettled;
  float *weights = graph->weights;

  int64_t node_count = graph->node_count;

  #pragma acc parallel loop default(present)
  for(int64_t i=0; i<node_count; i++) {
    if(frontier[i] == true) {
      int64_t begin = edge_index[i];
      int64_t end = edge_index[i+1];

      for(int64_t j=begin; j<end; j++) {
        int64_t tail_node = edge_tail[j];
        if(unsettled[tail_node] == true) {
          float tmp_dist = distance[i] + weights[j];
         if(tmp_dist < distance[tail_node]) {
            distance[tail_node] = tmp_dist;
            parents[tail_node] = i;
          }

        }
      }
    }
  }
}

void SettlementMin(struct SSSP *sssp, struct GRAPH *graph, float *g_delta_dist) {
//  sssp->delta_dist = FLT_MAX;

  float delta_dist = FLT_MAX;
  float *distance = sssp->distance;
  float *delta_node = sssp->delta_node;
  bool *unsettled = sssp->unsettled;

  #pragma acc parallel loop reduction(min:delta_dist) default(present)
  for(int64_t i=0; i<graph->node_count; i++) {
    if(unsettled[i] == true) {
      float dist_i = distance[i] + delta_node[i];

      delta_dist = fminf(delta_dist, dist_i);
    }
  }
 
  *g_delta_dist = delta_dist;
}


void SettlementUpdate(struct SSSP* sssp, struct GRAPH* graph, float delta_dist) {
  bool *frontier = sssp->frontier;
  bool *unsettled = sssp->unsettled;
  float *distance = sssp->distance;


  #pragma acc parallel loop default(present)
  for(int64_t i=0; i<graph->node_count; i++) {
    frontier[i] = false;
    if(unsettled[i] == true && distance[i] <= delta_dist) {
      unsettled[i] = false;
      frontier[i] = true;
    }
  }

}

// compute the minimum arch weight for each node
void PrecomputeNodeMinimum(struct SSSP* sssp, struct GRAPH* graph) {

  int64_t *edge_index = graph->edge_index;
  float   *weights = graph->weights;
  float *delta_node = sssp->delta_node;

  #pragma acc parallel loop default(present)
  for(int64_t i=0; i<graph->node_count; i++) {
    int64_t begin = edge_index[i];
    int64_t end = edge_index[i+1];

    float minimum_weight = FLT_MAX;
    #pragma acc loop seq
    for(int64_t j=begin; j<end; j++) {
      float weight = weights[j];
      if(weight < minimum_weight)
        minimum_weight = weight;
    }
    delta_node[i] = minimum_weight;
  }
}

void ReadInput(struct GRAPH *graph) {
  size_t bytes_read;
  printf("read forward star data...\n");
  FILE *ifp;
  if ((ifp = fopen("/home/atj/Documents/GISData/GDT/GDT_10MFS.bin", "r")) == NULL) {
    printf("open file failed\n");
  }
  bytes_read = fread(&graph->node_count, sizeof(int64_t), 1, ifp);
  printf("node count: %lu\n", graph->node_count);
  int record_size;
  bytes_read = fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->edge_index = SAFE_ALLOC(sizeof(int64_t), (graph->node_count + 1));
  bytes_read = fread(graph->edge_index, record_size, graph->node_count + 1, ifp);
  bytes_read = fread(&graph->arc_count, sizeof(int64_t), 1, ifp);
  printf("arc count: %lu\n", graph->arc_count);
  bytes_read = fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->edge_tail = SAFE_ALLOC(sizeof(int64_t), graph->arc_count);
  bytes_read = fread(graph->edge_tail, record_size, graph->arc_count, ifp);

  bytes_read = fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(float)) {
    printf("record size = %d != %lu", record_size, sizeof(float));
  }
  graph->weights = SAFE_ALLOC(sizeof(float), graph->arc_count);
  bytes_read = fread(graph->weights, sizeof(float), graph->arc_count, ifp);

  fclose(ifp);

  #pragma acc enter data copyin(graph[:1], \
                                graph->edge_index[0:(graph->node_count + 1)], \
                                graph->edge_tail[0:graph->arc_count], \
                                graph->weights[0:graph->arc_count])
}

void PrintSettledCount(struct SSSP sssp, struct GRAPH graph) {
  int64_t settled = 0;
  for(int64_t i=0; i<graph.node_count; i++) {
    if(sssp.unsettled[i] == false)
      settled++;
  }
  printf("Settled count: %lu \n", settled);
}
void CheckResults(struct SSSP sssp){
  float *distance = sssp.distance;
  int64_t *parents = sssp.parents;
  #pragma acc update host(distance[10000:1], parents[10000:1])

  printf("node %d distance: %f parent %lu\n", 10000, distance[10000], parents[10000]);
}

int main(int argc, char **argv) {
  struct GRAPH graph;
  struct SSSP sssp;

  float delta_dist;

  ReadInput(&graph);
  Initialize(&sssp, &graph, &delta_dist);

  PrecomputeNodeMinimum(&sssp, &graph);

  int64_t loop_count = 0;
  while(delta_dist < FLT_MAX) {
    printf("loop count: %lu\n", loop_count);
    Relaxation(&sssp, &graph);
    SettlementMin(&sssp, &graph, &delta_dist);
    SettlementUpdate(&sssp, &graph, delta_dist);

//    PrintSettledCount(sssp, graph);
    loop_count++;
  }

  CheckResults(sssp);

  return 0;
}
