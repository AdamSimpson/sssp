// Implimentation of "A New GPU-based Approach to the Shortest Path Problem"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct GRAPH {
  int64_t arc_count;
  int64_t node_count;

  // Using forward star representation
  // arcs[i] is first index in arc_tail and weights for i'th node
  // arcs[i+1] is the last index + 1 in arc_tail and weights for i'th node
  int64_t *arcs;
  int64_t *arc_tail;
  float  *weights;
};

struct SSSP {
  bool *unsettled;
  bool *frontier;
  float *distance;
  float *delta_node;
  int64_t start_node;
  float delta_dist;
};

void Initialize(struct SSSP *sssp, struct GRAPH graph) {
  sssp->unsettled = malloc(graph.node_count * sizeof(int64_t));
  sssp->frontier = malloc(graph.node_count * sizeof(int64_t));
  sssp->distance = malloc(graph.node_count * sizeof(float));
  sssp->delta_node = malloc(graph.node_count * sizeof(float));
  sssp->delta_dist = FLT_MAX;
  sssp->start_node = 0; // Change this to whatever

  for(int i=0; i<graph.node_count; i++) {
    sssp->distance[i] = FLT_MAX;
    sssp->frontier[i] = false;
    sssp->unsettled[i] = false;
  }

  sssp->distance[sssp->start_node] = 0.0;
  sssp->frontier[sssp->start_node] = true;
  sssp->unsettled[sssp->start_node] = false;
}

void Relaxation(struct SSSP sssp, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    if(sssp.frontier[i] == true) {
      int64_t begin = graph.arcs[i];
      int64_t end = graph.arcs[i+1];

      for(int64_t j=begin; j<end; j++) {
        int64_t tail_node = graph.arc_tail[j];
        if(sssp.unsettled[tail_node] == true) {
          sssp.distance[tail_node] = MIN(sssp.distance[tail_node], sssp.distance[i] + graph.weights[j]);
        }
      }
    }
  }

}

void SettlementMin(struct SSSP* sssp, struct GRAPH graph) {
  sssp->delta_dist = FLT_MAX;
  for(int i=0; i<graph.node_count; i++) {
    if(sssp->unsettled[i] == true) {
      float dist_i = sssp->distance[i] + sssp->delta_node[i];
      if(dist_i < sssp->delta_dist)
        sssp->delta_dist = dist_i;
    }
  }
}

void SettlementUpdate(struct SSSP sssp, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    sssp.frontier[i] = false;
    if(sssp.unsettled[i] == true && sssp.distance[i] <= sssp.delta_dist) {
      sssp.unsettled[i] = false;
      sssp.frontier[i] = true;
    }
  }
}

// compute the minimum arch weight for each node
void PrecomputeNodeMinimum(struct SSSP sssp, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    int64_t begin = graph.arcs[i];
    int64_t end = graph.arcs[i+1];

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
  if ((ifp = fopen("put/file/name/here", "r")) == NULL) {
    printf("open file failed\n");
  }
  fread(&graph->node_count, sizeof(int64_t), 1, ifp);
  int record_size;
  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->arcs = (int64_t *) malloc (sizeof(int64_t) * (graph->node_count + 1));
  fread(graph->arcs, record_size, graph->node_count + 1, ifp);
  fread(&graph->arc_count, sizeof(int64_t), 1, ifp);
  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(int64_t)) {
    printf("record size = %d != %lu", record_size, sizeof(int64_t));
  }
  graph->arc_tail = (int64_t *) malloc(sizeof(int64_t) * graph->arc_count);
  fread(graph->arc_tail, record_size, graph->arc_count, ifp);

  fread(&record_size, sizeof(int), 1, ifp);
  if (record_size != sizeof(float)) {
    printf("record size = %d != %lu", record_size, sizeof(float));
  }
  graph->weights = (float *) malloc (sizeof(float) * graph->arc_count);
  fread(graph->weights, sizeof(float), graph->arc_count, ifp);

  fclose(ifp);
}

int main(int argc, char **argv) {
  struct GRAPH graph;
  struct SSSP sssp;

  ReadInput(&graph);
  PrecomputeNodeMinimum(sssp, graph);
  Initialize(&sssp, graph);

  while(sssp.delta_dist < FLT_MAX) {
    Relaxation(sssp, graph);
    SettlementMin(&sssp, graph);
    SettlementUpdate(sssp, graph);
  }

  return 0;
}
