#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "barnes_hut.h"
#include <pthread.h>

#define TOPLEFT 1
#define TOPRIGHT 2
#define BOTTOMLEFT 3
#define BOTTOMRIGHT 4

#ifdef _MSC_VER
#define forceinline __forceinline
#elif defined(__GNUC__)
#define forceinline inline __attribute__((__always_inline__))
#elif defined(__CLANG__)
#if __has_attribute(__always_inline__)
#define forceinline inline __attribute__((__always_inline__))
#else
#define forceinline inline
#endif
#else
#define forceinline inline
#endif

static inline void insert_node(quad_tree_t *root, node_t *node);

static inline quad_tree_t *init_tree() {
  quad_tree_t *tree = calloc(1, sizeof(quad_tree_t));
  tree->TL_point = (coord_t) {0, 1};
  tree->BR_point = (coord_t) {1, 0};
  tree->node = NULL;
  return tree;
}

static inline void destroy_tree(quad_tree_t *root) {
  if (!root) return;
  destroy_tree(root->TL);
  destroy_tree(root->TR);
  destroy_tree(root->BL);
  destroy_tree(root->BR);
  free(root->node);
  free(root);
}
static inline void clear_tree(quad_tree_t *root) {
  destroy_tree(root->TL);
  root->TL = NULL;
  destroy_tree(root->TR);
  root->TR = NULL;
  destroy_tree(root->BL);
  root->BL = NULL;
  destroy_tree(root->BR);
  root->BR = NULL;
  free(root->node);
  root->TL_point = (coord_t) {0, 1};
  root->BR_point = (coord_t) {1, 0};
  root->internal = 0;
  root->node = NULL;

}

static inline int get_direction(quad_tree_t *root, node_t *node) {
  int left = node->xpos <= (root->TL_point.x + root->BR_point.x) / 2;
  int bottom = node->ypos <= (root->TL_point.y + root->BR_point.y) / 2;
  if (left) {
    return bottom ? BOTTOMLEFT : TOPLEFT;
  }
  return bottom ? BOTTOMRIGHT : TOPRIGHT;
}

static inline void update_tree(quad_tree_t *root, node_t *node) {
  double new_mass = root->node->mass + node->mass;
  root->node->xpos = (root->node->xpos * root->node->mass + node->xpos * node->mass) / new_mass;
  root->node->ypos = (root->node->ypos * root->node->mass + node->ypos * node->mass) / new_mass;
  root->node->mass = new_mass;
}

static inline void find_and_insert(quad_tree_t *root, node_t *node) {
  int direction = get_direction(root, node);
  switch (direction) {
    case BOTTOMLEFT:
      if (!root->BL) {
        quad_tree_t *BL = calloc(1, sizeof(quad_tree_t));
        BL->TL_point = (coord_t) {root->TL_point.x, (root->TL_point.y + root->BR_point.y) / 2};
        BL->BR_point = (coord_t) {(root->BR_point.x + root->TL_point.x) / 2, root->BR_point.y};
        BL->internal = 0;
        root->BL = BL;
      }
      insert_node(root->BL, node);
      return;
    case TOPLEFT:
      if (!root->TL) {
        quad_tree_t *TL = calloc(1, sizeof(quad_tree_t));
        TL->TL_point = root->TL_point;
        TL->BR_point = (coord_t) {(root->BR_point.x + root->TL_point.x) / 2, (root->TL_point.y + root->BR_point.y) / 2};
        TL->internal = 0;
        root->TL = TL;
      }
      insert_node(root->TL, node);
      return;
    case BOTTOMRIGHT:
      if (!root->BR) {
        quad_tree_t *BR = calloc(1, sizeof(quad_tree_t));
        BR->TL_point = (coord_t) {(root->BR_point.x + root->TL_point.x) / 2, (root->TL_point.y + root->BR_point.y) / 2};
        BR->BR_point = root->BR_point;
        BR->internal = 0;
        root->BR = BR;
      }
      insert_node(root->BR, node);
      return;
    case TOPRIGHT:
      if (!root->TR) {
        quad_tree_t *TR = calloc(1, sizeof(quad_tree_t));
        TR->TL_point = (coord_t) {(root->BR_point.x + root->TL_point.x) / 2, root->TL_point.y};
        TR->BR_point = (coord_t) {root->BR_point.x, (root->TL_point.y + root->BR_point.y) / 2};
        TR->internal = 0;
        root->TR = TR;
      }
      insert_node(root->TR, node);
      return;
    default:
      exit(-1);
  }
}

static inline void insert_node(quad_tree_t *root, node_t *node) {
  if (root->node == NULL) {
    root->node = node;
    return;
  }
  if (root->internal) {
    update_tree(root, node);
    find_and_insert(root, node);
    return;
  }
  find_and_insert(root, node);
  node_t *new_node = calloc(1, sizeof(node_t));
  new_node->xpos = root->node->xpos;
  new_node->ypos = root->node->ypos;
  new_node->mass = root->node->mass;
  update_tree(root, node);
  find_and_insert(root, new_node);
  root->internal = 1;
}

static inline void barnes_hut(quad_tree_t *tree, particle_t p, const double theta_max, const double e, vector_t *F) {
  if (tree->node != NULL) {
    node_t *node = tree->node;
    double width, x, y, distance;

    int notEqual = p.xpos != node->xpos && p.ypos != node->ypos && p.mass != node->mass;
    width = tree->BR_point.x - tree->TL_point.x;
    x = p.xpos - node->xpos;
    y = p.ypos - node->ypos;
    distance = sqrt(x * x + y * y);
    if ((!tree->internal && notEqual) || width / distance <= theta_max) {
      F->x += (node->mass / pow(distance + e, 3)) * x;
      F->y += (node->mass / pow(distance + e, 3)) * y;
    }
    else {
      if (tree->TL != NULL) barnes_hut(tree->TL, p, theta_max, e, F);
      if (tree->TR != NULL) barnes_hut(tree->TR, p, theta_max, e, F);
      if (tree->BL != NULL) barnes_hut(tree->BL, p, theta_max, e, F);
      if (tree->BR != NULL) barnes_hut(tree->BR, p, theta_max, e, F);
    }
  }
}

pthread_mutex_t stepLock;
pthread_cond_t stepReady;

//semaphores for synchronization
int waiting = 0, waitState = 0;

void barrier(int participants) {
  int myWaitState;
  pthread_mutex_lock(&stepLock);
  myWaitState = waitState;
  waiting++;
  if (waiting == participants) {
    waiting = 0;
    waitState = 1 - myWaitState;
    pthread_cond_broadcast(&stepReady);
  }
  while (myWaitState == waitState) pthread_cond_wait(&stepReady, &stepLock);
  pthread_mutex_unlock(&stepLock);
}

typedef struct payload {
    particle_t *particles;
    quad_tree_t *tree;
    int start;
    int stop;
    int participants;
    int nsteps;
    double G;
    double theta;
    double e;
    double delta_t;

} payload_t;

void *update_particles(void *args) {
  payload_t *payload = (payload_t *) args;
  vector_t F;
  for(int i = 0; i < payload->nsteps; i++){
    //wait for clear and insert
    barrier(payload->participants);
    for (int j = payload->start; j < payload->stop; j++) {
      F.x = 0;
      F.y = 0;
      barnes_hut(payload->tree, payload->particles[j], payload->theta, payload->e, &F);
      F.x *= -(payload->G);
      F.y *= -(payload->G);
      particle_t p = payload->particles[j];
      p.xvel = p.xvel + payload->delta_t * F.x;
      p.yvel = p.yvel + payload->delta_t * F.y;
      p.xpos = p.xpos + payload->delta_t * p.xvel;
      p.ypos = p.ypos + payload->delta_t * p.yvel;
      payload->particles[j] = p;
    }
    //broadcast done with updates/wait for others
    barrier(payload->participants);
  }


  pthread_exit(NULL);
}

int main(const int argc, const char *argv[]) {
  if (argc != 7) {
    printf("Program takes 6 arguments!\n");
    return 1;
  }
  const double N = atof(argv[1]);
  const unsigned int nsteps = atoi(argv[3]);
  const double delta_t = atof(argv[4]);
  const double theta = atof(argv[5]);
  const double G = (double) 100 / N;
  const double e = 0.001;

  int noThreads = 10, participants = noThreads + 1;

  particle_t particles[(int) N];
  FILE *f = fopen(argv[2], "r");
  fread(particles, sizeof(particle_t), N, f);
  fclose(f);


  quad_tree_t *tree = init_tree();
  //initializes threads, locks and conditions
  pthread_t threads[noThreads];
  pthread_attr_t attr;
  pthread_cond_init(&stepReady, NULL);
  pthread_mutex_init(&stepLock, NULL);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  //call threads with payloads
  for (int i = 0; i < noThreads; i++) {
    payload_t *payload = calloc(1, sizeof(payload_t));
    payload->start = i * (int)N / noThreads;
    payload->stop = (i + 1) * (int)N / noThreads;
    payload->G = G;
    payload->delta_t = delta_t;
    payload->e = e;
    payload->particles = particles;
    payload->tree = tree;
    payload->theta = theta;
    payload->participants = participants;
    payload->nsteps = nsteps;
    pthread_create(&threads[i], NULL, update_particles, (void *) payload);
  }

  for (int i = 0; i < nsteps; i++) {
    //insert particles into the tree
    for (int k = 0; k < N; k++) {
      node_t *node = calloc(1, sizeof(node_t));
      node->xpos = particles[k].xpos;
      node->ypos = particles[k].ypos;
      node->mass = particles[k].mass;
      insert_node(tree, node);
    }
    //broadcast done
    barrier(participants);
    //wait for particle update
    barrier(participants);
    //clear the tree
    clear_tree(tree);
  }
  void *status;
  for (int i = 0; i < noThreads; i++)
    pthread_join(threads[i], &status);

  FILE *file = fopen("result.gal", "w");
  fwrite(particles, sizeof(particle_t), N, file);
  fclose(file);
  return 0;
}
