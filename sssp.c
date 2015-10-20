// Implimentation of "A New GPU-based Approach to the Shortest Path Problem"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include "safe_alloc.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct GRAPH {
  int64_t arc_count;
  int64_t node_count;

  // Using forward star representation
  // edge_index[i] is first index in edge_tail and weights for i'th node
  // edge_index[i+1] is the last index + 1 in edge_tail and weights for i'th node
  int64_t *restrict edge_index;
  int64_t *restrict edge_tail;
  float  *restrict weights;
};

struct SSSP {
  bool *restrict unsettled;
  bool *restrict frontier;
  float *restrict distance;
  float *restrict delta_node;
  int64_t *restrict parents; 
  int64_t start_node;
  float delta_dist;
};

void Initialize(struct SSSP *sssp, struct GRAPH graph) {
  sssp->unsettled  = SAFE_ALLOC(graph.node_count, sizeof(int64_t));
  sssp->frontier   = SAFE_ALLOC(graph.node_count, sizeof(int64_t));
  sssp->distance   = SAFE_ALLOC(graph.node_count, sizeof(float));
  sssp->delta_node = SAFE_ALLOC(graph.node_count, sizeof(float));
  sssp->parents    = SAFE_ALLOC(graph.node_count, sizeof(int64_t));

  sssp->delta_dist = 0.0f;
  sssp->start_node = 0;

  for(int64_t i=0; i<graph.node_count; i++) {
    sssp->distance[i] = FLT_MAX;
    sssp->frontier[i] = false;
    sssp->unsettled[i] = true;
    sssp->parents[i] = -1;
  }

  sssp->distance[sssp->start_node] = 0.0f;
  sssp->frontier[sssp->start_node] = true;
  sssp->unsettled[sssp->start_node] = false;
  sssp->parents[sssp->start_node] = sssp->start_node;
}

void Relaxation(struct SSSP sssp, struct GRAPH graph) {
  for(int64_t i=0; i<graph.node_count; i++) {
    if(sssp.frontier[i] == true) {
      int64_t begin = graph.edge_index[i];
      int64_t end = graph.edge_index[i+1];

      for(int64_t j=begin; j<end; j++) {
        int64_t tail_node = graph.edge_tail[j];
        if(sssp.unsettled[tail_node] == true) {
          float tmp_dist = sssp.distance[i] + graph.weights[j];
          if(tmp_dist < sssp.distance[tail_node]) {
            sssp.distance[tail_node] = tmp_dist;
            sssp.parents[tail_node] = i;
 
          }
            
        }
      }
    }
  }
}

void SettlementMin(struct SSSP *sssp, struct GRAPH graph) {
  sssp->delta_dist = FLT_MAX;
  for(int64_t i=0; i<graph.node_count; i++) {
    if(sssp->unsettled[i] == true) {
      float dist_i = sssp->distance[i] + sssp->delta_node[i];
      if(dist_i < sssp->delta_dist)
        sssp->delta_dist = dist_i;
    }
  }
}

void SettlementUpdate(struct SSSP sssp, struct GRAPH graph) {
  for(int64_t i=0; i<graph.node_count; i++) {
    sssp.frontier[i] = false;
    if(sssp.unsettled[i] == true && sssp.distance[i] <= sssp.delta_dist) {
      sssp.unsettled[i] = false;
      sssp.frontier[i] = true;
    }
  }
}

// compute the minimum arch weight for each node
void PrecomputeNodeMinimum(struct SSSP sssp, struct GRAPH graph) {
  for(int64_t i=0; i<graph.node_count; i++) {
    int64_t begin = graph.edge_index[i];
    int64_t end = graph.edge_index[i+1];

    float minimum_weight = FLT_MAX;
    for(int64_t j=begin; j<end; j++) {
      float weight = graph.weights[j];
      if(weight < minimum_weight)
        minimum_weight = weight;
    }
    sssp.delta_node[i] = minimum_weight;
  }
}

void ReadInput(struct GRAPH *graph) {
  printf("read forward star data...\n");
  FILE *ifp;
  if ((ifp = fopen("/home/atj/Documents/GISData/GDT/GDT_10MFS.bin", "r")) == NULL) {
    printf("open file failed\n");
  }
  fread(&graph->node_count, sizeof(int64_t), 1, ifp);
  printf("node count: %lu\n", graph->node_count);
  int record_size;
  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->edge_index = SAFE_ALLOC(sizeof(int64_t), (graph->node_count + 1));
  fread(graph->edge_index, record_size, graph->node_count + 1, ifp);
  fread(&graph->arc_count, sizeof(int64_t), 1, ifp);
  printf("arc count: %lu\n", graph->arc_count);
  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->edge_tail = SAFE_ALLOC(sizeof(int64_t), graph->arc_count);
  fread(graph->edge_tail, record_size, graph->arc_count, ifp);

  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(float)) {
    printf("record size = %d != %lu", record_size, sizeof(float));
  }
  graph->weights = SAFE_ALLOC(sizeof(float), graph->arc_count);
  fread(graph->weights, sizeof(float), graph->arc_count, ifp);

  fclose(ifp);
}

void printSettledCount(struct SSSP sssp, struct GRAPH graph) {
  int64_t settled = 0;
  for(int64_t i=0; i<graph.node_count; i++) {
    if(sssp.unsettled[i] == false)
      settled++;
  }
  printf("Settled count: %lu \n", settled);
}

int main(int argc, char **argv) {
  struct GRAPH graph;
  struct SSSP sssp;

  ReadInput(&graph);
  Initialize(&sssp, graph);
  PrecomputeNodeMinimum(sssp, graph);

  int64_t loop_count = 0;
  while(sssp.delta_dist < FLT_MAX) {
//    printf("loop count: %lu\n", loop_count);
    Relaxation(sssp, graph);
    SettlementMin(&sssp, graph);
    SettlementUpdate(sssp, graph);

//    printSettledCount(sssp, graph);
    loop_count++;
  }

  printf("node %d distance: %f parent %lu\n", 1000, sssp.distance[1000], sssp.parents[1000]);

  return 0;
}
