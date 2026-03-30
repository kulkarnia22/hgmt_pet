#ifndef imager_h
#define imager_h

#include <stdbool.h>
// including custom files
#include "../src/compton_chain_ordering.h"
#include "../src/helper_functions.h"
#include "../src/hgmt_structs.h"
#include "../src/linear_algebra.h"
#define X_LENGTH 10.0 //cm
#define Y_LENGTH 10.0 //cm
#define Z_LENGTH 4.0 //cm
#define X_RES 100 //pixels
#define Y_RES 100 //pixels
#define Z_RES 1 //pixels
#define RES (X_RES * Y_RES * Z_RES)
typedef struct grid {
  uint counts[RES];
} grid;
typedef struct graph {
  double values[RES];
  double background;
} graph;
typedef struct comp_lor {
  vec3d center;
  sym_matrix new_dot;
  double constant;
  int i_left;
  int i_right;
  int j_left;
  int j_right;
  int k_left;
  int k_right;

} comp_lor;
typedef struct thread_data {
  uint num_lors;
  comp_lor *data;
  graph *current_image;
} thread_data;
#endif
