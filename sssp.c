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
  int64_t *arcs; // arcs[i] is index in arc_begin and weights for i'th node
  int64_t *arc_begin;
  float  *weights;
};

void Initialize(bool *unsettled, bool *frontier, float *distance,
                int64_t start_node, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    distance[i] = FLT_MAX;
    frontier[i] = false;
    unsettled[i] = false;
  }

  distance[start_node] = 0.0;
  frontier[start_node] = true;
  unsettled[start_node] = false;
}

void Relaxation(bool *unsettled, bool *frontier, float *distance, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    if(frontier[i] == true) {
      int64_t arc = graph.arcs[i];
      int64_t arc_next = graph.arcs[i+1];
      int64_t begin = graph.arc_begin[arc];
      int64_t end = graph.arc_begin[arc_next];
      for(int64_t j=begin; j<end; j++) {
        if(unsettled[j] == true) {
          distance[j] = MIN(distance[j], distance[i] + graph.weights[j]);
        }
      }
    }
  }

}

void SettlementMin(bool *unsettled, float *distance,
                   float *delta_node, float *delta_dist, struct GRAPH graph) {
  *delta_dist = FLT_MAX;
  for(int i=0; i<graph.node_count; i++) {
    if(unsettled[i] == true) {
      float dist_i = distance[i] + delta_node[i];
      if(dist_i < *delta_dist)
        *delta_dist = dist_i;
    }
  }
}

void SettlementUpdate(bool *unsettled, bool *frontier,
                      float *distance, float delta_dist, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    frontier[i] = false;
    if(unsettled[i] == true && distance[i] <= delta_dist) {
     unsettled[i] = false;
     frontier[i] = true;
    }
  }
}

// compute the minimum arch weight for each node
void PrecomputeNodeMinimum(float *delta_node, struct GRAPH graph) {
  for(int i=0; i<graph.node_count; i++) {
    int64_t arc = graph.arcs[i];
    int64_t arc_next = graph.arcs[i+1];
    int64_t begin = graph.arc_begin[arc];
    int64_t end = graph.arc_begin[arc_next];

    float minimum_weight = FLT_MAX;
    for(int64_t j=begin; j<end; j++) {
      float weight = graph.weights[j];
      if(weight < minimum_weight)
        minimum_weight = weight;
    }
    delta_node[i] = minimum_weight;
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
  graph->arc_begin = (int64_t *) malloc(sizeof(int64_t) * graph->arc_count);
  fread(graph->arc_begin, record_size, graph->arc_count, ifp);

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

  ReadInput(&graph);
  bool *unsettled = malloc(graph.node_count * sizeof(int64_t));
  bool *frontier = malloc(graph.node_count * sizeof(int64_t));
  float *distance = malloc(graph.node_count * sizeof(float));
  float *delta_node = malloc(graph.node_count * sizeof(float));
  float delta_dist = FLT_MAX;

  PrecomputeNodeMinimum(delta_node,graph);

  int64_t start_node = 0;
  Initialize(unsettled, frontier, distance, start_node, graph);
  while(delta_dist < FLT_MAX) {
    Relaxation(unsettled, frontier, distance, graph);
    SettlementMin(unsettled, distance, delta_node, &delta_dist, graph);
    SettlementUpdate(unsettled, frontier, distance, delta_dist, graph);
  }

  return 0;
}
