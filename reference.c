#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define NONE -1
#define CLOSED -2

long long numOfNodes;
long long numOfLinks;
long long *aNodes;
long long *bNodes;
float *impedances;
char forwardStarBinFileName[2048];

long long **parents;
long long **heap;
long long **qpos;
float **minValues;

#define N 1

void shortest_path_initialize() {
	printf("begin: shortest_path_initialize\n");
	heap = (long long **) calloc(N, sizeof(long long *));
 	minValues = (float **)calloc(N, sizeof(float *));
 	parents = (long long **)calloc(N, sizeof(long long *));
 	qpos = (long long **)calloc(N, sizeof(long long *));
	int i;
	for (i = 0; i < N; ++i) {
  		heap[i] = (long long *)calloc(numOfNodes, sizeof(long long));
  		minValues[i] = (float *)calloc(numOfNodes, sizeof(float));
  		parents[i] = (long long *)calloc(numOfNodes, sizeof(long long));
  		qpos[i] = (long long *)calloc(numOfNodes, sizeof(long long));
	}
	int j;
 	for (i = 0; i < N; ++i) {
 		for (j = 0; j < numOfNodes; ++j) {
 			minValues[i][j] = FLT_MAX;
    	parents[i][j] = NONE;
    	qpos[i][j] = NONE;
		}
  }
	printf("end: shortest_path_initialize\n");
}

void shortest_path_execute(long long aCarID, long long aStartNode) {
  printf("shortest path...\n");
  long long current;
  long long hp1;
  long long hp2;
  long long hp3;
  long long j;
  long long k;
  long long k2;
  long long next;
  float nextValue;
  long long nhp = 0;;
  long long theTotalNodesVisited = 0;

  minValues[aCarID][aStartNode] = 0.0;
  parents[aCarID][aStartNode] = NONE;
  heap[aCarID][nhp++] = aStartNode;
  current = NONE;
  while (nhp > 0) {
    current = heap[aCarID][0];
    qpos[aCarID][current] = CLOSED;
    if (nhp > 1) {
      hp1 = heap[aCarID][--nhp];
      nextValue = minValues[aCarID][hp1];
      ++theTotalNodesVisited;
      k = 0;
      k2 = k * 2 + 1;
      while (k2 < nhp) {
        hp2 = heap[aCarID][k2];
        if (k2 < nhp - 1) {
          hp3 = heap[aCarID][k2 + 1];
          if (minValues[aCarID][hp3] < minValues[aCarID][hp2]) {
            hp2 = hp3;
            k2++;
          }
        }
        if (nextValue > minValues[aCarID][hp2]) {
          heap[aCarID][k] = hp2;
          qpos[aCarID][hp2] = k;
          k = k2;
          k2 = k * 2 + 1;
        } else
          break;
      }
      heap[aCarID][k] = hp1;
      qpos[aCarID][hp1] = k;
    } else
      nhp--;

    for (j = aNodes[current]; j < aNodes[current + 1]; ++j) {
      next = bNodes[j];
      if (qpos[aCarID][next] == CLOSED)
        continue;
      nextValue = minValues[aCarID][current] + impedances[j];
      if (nextValue < minValues[aCarID][next]) {
        minValues[aCarID][next] = nextValue;
        parents[aCarID][next] = current;
        if (qpos[aCarID][next] == NONE)
          qpos[aCarID][next] = nhp++;
        k = qpos[aCarID][next];
        while (k > 0) {
          k2 = (k - 1) / 2;
          hp2 = heap[aCarID][k2];
          if (nextValue < minValues[aCarID][hp2]) {
            heap[aCarID][k] = hp2;
            qpos[aCarID][hp2] = k;
            k = k2;
          } else
            break;
        }
        heap[aCarID][k] = next;
        qpos[aCarID][next] = k;
      }
    }
  }
  printf("total nodes visited = %ld\n", theTotalNodesVisited);
}

void read_forward_star_from_bin_file() {
  printf("begin:read forward star data...\n");
  FILE *ifp;
  if ((ifp = fopen(forwardStarBinFileName, "r")) == NULL) {
    printf("open %s failed", forwardStarBinFileName);
  }
  fread(&numOfNodes, sizeof(long long), 1, ifp);
  int theRecordSize;
  fread(&theRecordSize, sizeof(int), 1, ifp);
  if (theRecordSize != sizeof(long long)) {
    printf("record size = %d != %d", theRecordSize, sizeof(long long));
  }
  aNodes = (long long *) malloc (sizeof(long long) * (numOfNodes + 1));
  fread(aNodes, theRecordSize, numOfNodes + 1, ifp);
  fread(&numOfLinks, sizeof(long long), 1, ifp);
  fread(&theRecordSize, sizeof(int), 1, ifp);
  if (theRecordSize != sizeof(long long)) {
    printf("record size = %d != %d", theRecordSize, sizeof(long long));
  }
  bNodes = (long long *) malloc(sizeof(long long) * numOfLinks);
  fread(bNodes, theRecordSize, numOfLinks, ifp);

  fread(&theRecordSize, sizeof(int), 1, ifp);
  if (theRecordSize != sizeof(float)) {
    printf("record size = %d != %d", theRecordSize, sizeof(float));
  }
  impedances = (float *) malloc (sizeof(float) * numOfLinks);
  fread(impedances, sizeof(float), numOfLinks, ifp);

  fclose(ifp);
  printf("end:read forward star data...\n");

}

void dump_forward_star() {
  printf("number of nodes = %d", numOfNodes);
  int i, j;
  for (i = 0; i < numOfNodes; ++i) {
    for (j = aNodes[i]; j < aNodes[i + 1]; ++j) {
      printf("%d - %ld", i, bNodes[j]);
    }
  }
}

void check() {
//	if (fabs(minValues[0][10000] - 59.544136) < 0.00001) 
//		printf ("check OK\n");
//	else
//		printf("check failed\n");
	printf ("min_value %f, parent: %lu\n", minValues[0][10000], parents[0][10000]);
}

void usage(int argc, char *argv[]) {
  if (argc < 4) {
    printf ("usage: %s <project> <input> <root> ", argv[0]);
    exit(1);
  }
}

int main(int argc, char **argv) {
/*
  usage(argc, argv);

  char *theEnvironment = getenv("GISDATAHOME");
  if (theEnvironment == NULL) {
    printf("Environment variable GISDATAHOME not set");
    exit(1);
  }
  char * theGISDataHome = theEnvironment;
  char * theProject = argv[1];
*/
//  char theInputFileName[128] = "/home/atj/Documents/GISData/GDT/GDT_10MFS.bin";//argv[2];
  long long theRoot = 0;//atol(argv[3]);
  sprintf (forwardStarBinFileName, "%s", "/home/atj/Documents/GISData/GDT/GDT_10MFS.bin");

  printf("input       = %s\n", forwardStarBinFileName);


	read_forward_star_from_bin_file();
 	shortest_path_initialize();
	int i;
	for (i = 0; i < N; ++i) {
  		shortest_path_execute(i, i);
	}


	check();

  exit(0);
}
