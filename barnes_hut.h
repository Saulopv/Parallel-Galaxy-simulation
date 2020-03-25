#ifndef BARNES_HUT_H
#define BARNES_HUT_H

typedef struct vector {
    double x;
    double y;
} vector_t;

typedef struct particle {
    double xpos;
    double ypos;
    double mass;
    double xvel;
    double yvel;
    double brightness;
} particle_t;

typedef struct coord {
  double x, y;
} coord_t;

typedef struct node {
  double xpos, ypos, mass;
} node_t;

typedef struct quad_tree {
  struct quad_tree *TL;
  struct quad_tree *TR;
  struct quad_tree *BL;
  struct quad_tree *BR;
  node_t *node;
  coord_t TL_point;
  coord_t BR_point;
  int internal;
} quad_tree_t;

#endif


