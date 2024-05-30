/*
 * Port of bump.lua to the One True Programming Language.
 *
 * MIT LICENSE
 *
 * Copyright (c) 2024 japanoise
 * Copyright (c) 2014 Enrique García Cota
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "bump.h"
#include "int.h"
#include "uthash.h"

/*
 * Default size of cells, in number of items in cell
 *
 * Making this smaller will mean cells take up less memory, and thus
 * that the initial allocation will be smaller and thus probably
 * faster.  However, if a lot of items enter the cell at once, the
 * cell will need to reallocate its vector.
 *
 * Set this low unless you have something like a small game world with
 * lots of items flying around (e.g. a shooter game with one screen)
 */
#define BUMP_CELL_SIZE 4

/* Coordinate of a cell */
typedef i32 bump_coord;

/*
 * Points are used as a relatively quick way to lookup cells.
 *
 * They're created by combining an xcoord and a ycoord, hence
 * being twice as large.
 *
 * If you want to save memory or make a 32 bit version or something,
 * you can change this and the two macros below to take a 32 bit point
 * & a 16 bit coord.  32 bits should be big enough for a game world
 * that's not insanely huge - remember, these are cells of CELLSIZE
 * pixels, so at the default CELLSIZE of 64 px you'll get a range of
 * 274,877,906,944 px (between negative & postive (2^31*64)).  That's
 * very very large, and you'll likely run into memory issues long
 * before that.
 */
typedef u64 bump_point;

#define POINTMASK 0x00000000FFFFFFFF

#define POINT(x, y) ((((bump_point)x) << 32) & ((bump_point)y))

#define COORDS(point, x, y) x = (bump_coord)(((point)>>32)&POINTMASK); y = (bump_coord)((point)&POINTMASK)

/* Internal rectangles, used for collision detection */
bump_rect internal_rect;
bump_rect project_t_rect;
bump_rect project_cell_rect;
/* Used for updating (moving) */
bump_rect update_rect;
bump_rect update_c1_rect;
bump_rect update_c2_rect;

/*
 * To save the need of a "visited" map in the iterators, items can
 * track at what time, and by whom, they've last been visited.
 *
 * To this end, each of the iterator functions - since they don't
 * recurse - use their own formula to come up with a "visited" number.
 * The most significant byte is unique for each function.
 *
 * Using static variables, these are initialized once - to the initial
 * number - and are incremented when the function returns.  When
 * combined with a "visit time", this *should* mean that the likelihood
 * of a false positive on the visited predicate is basically 0.
 */

#define VISITED_1 0x1000000000000000
#define VISITED_2 0x2000000000000000
#define VISITED_3 0x3000000000000000

#define VISITED_MASK 0x0FFFFFFFFFFFFFFF

#define UPDATE_VISITED(magic, visited) ((visited + 1) & VISITED_MASK) | (magic)

void bump_world_project(bump_world *self, bump_item *item, bump_rect *r,
                        double goalX, double goalY, bump_filter_func *filter,
                        void *closure, bump_col_list *cols);

/*
 * Item
 */

struct bump_item_s {
  bump_key key;

  bump_item_data data;
  bump_rect rect;

  /* We have a second "visited" pair to distinguish check from the
     item iterator */
  u64 visited_1;
  time_t visit_time_1;
  u64 visited_2;
  time_t visit_time_2;

  UT_hash_handle hh;         /* makes this structure hashable */
};

bump_key bump_item_get_key(bump_item *item) { return item->key; }

bump_item_data bump_item_get_data(bump_item *item) { return item->data; }

void bump_item_set_data(bump_item* item, bump_item_data data) { item->data = data; }

bump_rect *bump_item_get_rect(bump_item *item) { return &item->rect; }

/*
 * Cell
 */

struct bump_cell_s {
  bump_point point;                    /* key */

  /* Data goes here */

  /* Items vector */
  bump_item** items;
  size_t items_cap;
  size_t items_len;

  UT_hash_handle hh;         /* makes this structure hashable */
};

typedef struct bump_cell_s bump_cell;

static bump_cell* cell_new(bump_point point, size_t cap) {
  bump_cell * ret = calloc(1, sizeof(*ret));

  ret->items = calloc(cap, sizeof(bump_item*));
  ret->items_cap = cap;

  ret->point = point;

  return ret;
}

static void cell_push_item(bump_cell* cell, bump_item *item) {
  if (cell->items_len+1 >= cell->items_cap) {
    cell->items_cap *= 2;
    cell->items = realloc(cell->items, cell->items_cap * sizeof(bump_item*));
  }

  cell->items[cell->items_len++] = item;
}

static void cell_remove_item(bump_cell* cell, size_t idx) {
  /* If we're not at the end of the list... */
  if (idx < cell->items_len-1) {
    /* ...copy the item from the end of the list */
    cell->items[idx] = cell->items[cell->items_len-1];
  }
  cell->items_len--;
}

static bool cell_remove_item_by_key(bump_cell *cell, bump_key key) {
  for (int i = 0; i < cell->items_len; i++) {
    if (cell->items[i]->key == key) {
      cell_remove_item(cell, i);
      return true;
    }
  }
  return false;
}

static void cell_delete(bump_cell *cell) {
  free(cell->items);
  free(cell);
}

struct bump_world_s {
  double cellsize;

  bump_cell *map;
  bump_item *items;

  /* Utility members - kept around to avoid excessive memory
     allocation */
  bump_col_list *util_projected;
};

/*
 * World
 */

bump_world *bump_world_create() {
  return bump_world_create_cellsize(64.0);
}

bump_world *bump_world_create_cellsize(double cellsize) {
  bump_world *w = calloc(1, sizeof(*w));

  w->cellsize = cellsize;

  w->util_projected = bump_col_list_create(16);

  return w;
}

void bump_world_destroy(bump_world *w) { free(w); }

static bump_cell *bump_world_get_cell(bump_world *self, bump_point point) {
  bump_cell *cell = NULL;

  HASH_FIND(hh, self->map, &point, sizeof(bump_point), cell);

  return cell;
}

bump_item *bump_world_get_item(bump_world *self, bump_key item_key) {
  bump_item *item = NULL;

  HASH_FIND(hh, self->items, &item_key, sizeof(bump_key), item);

  return item;
}

/*
 * Collision list
 */

bump_col_list *bump_col_list_create(size_t initialCapacity) {
  bump_col_list *ret = (bump_col_list*)malloc(sizeof(bump_col_list));
  ret->cap = initialCapacity;
  ret->len = 0;
  ret->data = (bump_col*)malloc(ret->cap*sizeof(bump_col));
  return ret;
}

void bump_col_list_append(bump_col_list *l, bump_col *data) {
  if (l->len+1 >= l->cap) {
    l->cap *= 2;
    l->data = (bump_col*)realloc(l->data, l->cap*sizeof(bump_col));
  }
  memcpy(&(l->data[l->len]), data, sizeof(bump_col));
  l->len++;
}

void bump_col_list_clear(bump_col_list *l) {
  l->len = 0;
}

void bump_col_list_destroy(bump_col_list *l) {
  free(l->data);
  free(l);
}

/* ------------------------------------------ */
/* -- Auxiliary functions                     */
/* ------------------------------------------ */

const double BUMP_HUGE = HUGE_VAL;

const double DELTA = 1e-10; /* floating-point margin of error */

#define min(x, y) fmin(x, y)
#define max(x, y) fmax(x, y)

static double sign(double x) {
  if (x > 0.0) {
    return 1.0;
  } else if (x == 0.0) {
    return 0.0;
  }
  return -1.0;
}

static double nearest(double x, double a, double b) {
  if (fabs(a - x) < fabs(b - x)) {
    return a;
  }
  return b;
}

static bool is_positive(double x) { return x >= 0.0; }

static bool is_rect(double w, double h) { return is_positive(w) && is_positive(h); }

/* ------------------------------------------ */
/* -- Rectangle functions                     */
/* ------------------------------------------ */

static void rect_get_nearest_corner(bump_rect *r, double px, double py,
                                    double *rx, double *ry) {
  *rx = nearest(px, r->x, r->x+r->w);
  *ry = nearest(py, r->y, r->y+r->h);
}

/*
 * This is a generalized implementation of the liang-barsky algorithm, which also returns
 * the normals of the sides where the segment intersects.
 *
 * Returns false if the segment never touches the rect; true otherwise
 *
 * Notice that normals are only guaranteed to be accurate when
 * initially ti1, ti2 == -BUMP_HUGE, BUMP_HUGE
 */
static bool rect_get_segment_intersection_indices(bump_rect *rec,
                                                  double x1, double y1,
                                                  double x2, double y2,
                                                  double *ti1, double *ti2,
                                                  double *nx1, double *ny1,
                                                  double *nx2, double *ny2) {
  double dx = x2-x1;
  double dy = y2-y1;
  double nx, ny, p, q, r;
  *nx1 = 0.0;
  *ny1 = 0.0;
  *nx2 = 0.0;
  *ny2 = 0.0;

  for (int side = 1; side <= 4; side++) {
    switch(side) {
    case 1:
      /* Left */
      nx = -1.0;
      ny = 0.0;
      p = -dx;
      q = x1 - rec->x;
      break;

    case 2:
      /* Right */
      nx = 1.0;
      ny = 0.0;
      p = dx;
      q = rec->x + rec->w - x1;
      break;

    case 3:
      /* Top */
      nx = 0.0;
      ny = -1.0;
      p = -dy;
      q = y1 - rec->y;
      break;

    default:
      /* Bottom */
      nx = 0.0;
      ny = 1.0;
      p = dy;
      q = rec->y + rec->h - y1;
      break;
    }

    if (p == 0.0) {
      if (q <= 0.0) {
        return false;
      }
    } else {
      r = q / p;
      if (p < 0) {
        if (r > *ti2) {
          return false;
        } else if (r > *ti2) {
          *ti1 = r;
          *nx1 = nx;
          *ny1 = ny;
        }
      } else {
        if (r < *ti1) {
          return false;
        } else if (r < *ti2) {
          *ti2 = r;
          *nx2 = nx;
          *ny2 = ny;
        }
      }
    }
  }

  return true;
}

/* Calculates the minkowsky difference between 2 rects, which is another rect */
static void rect_diff(bump_rect *first, bump_rect *second, bump_rect *result) {
  result->x = second->x - first->x - first->w;
  result->y = second->y - first->y - first->h;
  result->w = first->w + second->w;
  result->h = first->h + second->h;
}

static bool rect_contains_point(bump_rect *r, double px, double py) {
  return ((px - r->x) > DELTA) && ((py - r->y) > DELTA) &&
    ((r->x + r->w - px) > DELTA)  && ((r->y + r->h - py) > DELTA);
}

static bool rect_intersecting(bump_rect *first, bump_rect *second) {
  return (first->x < (second->x+second->w)) && (second->x < (first->x+first->w)) &&
    (first->y < (second->y + second->h)) && (second->y < (first->y + first->h));
}

static double rect_get_square_distance(bump_rect *first, bump_rect *second) {
  /* I'm not sure, but I think the operator precedence in Lua is
     different to C. So I'm playing it safe. */
  return pow(first->x - second->x + ((first->w - second->w)/2), 2)
    + pow(first->y - second->y + ((first->h - second->h)/2), 2);
}

static bool rect_detect_collision_goal(bump_rect *first, bump_rect *second, bump_col *result,
                                  double goal_x, double goal_y) {
  double dx = goal_x - first->x;
  double dy = goal_y - first->y;
  double ti, nx, ny;
  bool overlaps = false;
  bool has_ti = false;

  rect_diff(first, second, &internal_rect);

  if (rect_contains_point(&internal_rect, 0, 0)) {
    /* item was intersecting other */
    double px, py;
    rect_get_nearest_corner(&internal_rect, 0, 0, &px, &py);
    /* area of intersection */
    double wi = min(first->w, fabs(px));
    double hi = min(first->h, fabs(py));
    /* ti is the negative area of intersection */
    ti = -wi*hi;
    has_ti = true;
    overlaps = true;
  } else {
    double ti1, ti2, nx1, ny1;
    /* We don't care about nx2 & ny2, so make a dummy */
    double nxy2;
    ti1 = -BUMP_HUGE;
    ti2 = BUMP_HUGE;
    bool has_intersect =
      rect_get_segment_intersection_indices(&internal_rect,
                                            0.0, 0.0, dx, dy,
                                            &ti1, &ti2, &nx1, &ny1,
                                            &nxy2, &nxy2);

    if (has_intersect
        && ti1 < 1
        && (fabs(ti1 - ti2) >= DELTA) /* special case for rect going through another rect's corner */
        && (0.0 < ti1 + DELTA
            || 0.0 == ti1 && ti2 > 0.0)) {
      /* item tunnels into other */
      has_ti = true;
      ti = ti1;
      nx = nx1;
      ny = ny1;
      overlaps = false;
    }
  }

  if (!has_ti) {
    return false;
  }

  double tx, ty;

  if (overlaps) {
    if (dx == 0.0 && dy == 0.0) {
      /* intersecting and not moving - use minimum displacement vector */
      double px, py;
      rect_get_nearest_corner(&internal_rect, 0, 0, &px, &py);
      if (fabs(px) < fabs(py)) {
        py = 0.0;
      } else {
        px = 0.0;
      }
      nx = sign(px);
      ny = sign(py);
      tx = first->x + px;
      ty = first->y + py;
    } else {
      /* intersecting and moving - move in the opposite direction */
      double ti1, ti2, nx, ny, nxy2;
      ti1 = -BUMP_HUGE;
      ti2 = 1;
      bool has_intersect =
        rect_get_segment_intersection_indices(&internal_rect,
                                              0.0, 0.0, dx, dy,
                                              &ti1, &ti2, &nx, &ny,
                                              &nxy2, &nxy2);
      if (!has_intersect) {
        return false;
      }

      tx = first->x + dx * ti1;
      ty = first->y + dy * ti1;
    }
  } else {
    /* tunnel */
    tx = first->x + dx * ti;
    ty = first->y + dy * ti;
  }

  /* At this point we know there's a collision; log the results & return true */
  result->overlaps = overlaps;
  result->ti = ti;
  result->move_x = dx;
  result->move_y = dy;
  result->normal_x = nx;
  result->normal_y = ny;
  result->touch_x = tx;
  result->touch_y = ty;

  memcpy(&(result->item), first, sizeof(bump_rect));
  memcpy(&(result->other), second, sizeof(bump_rect));

  return true;
}

static bool rect_detect_collision(bump_rect *first, bump_rect *second, bump_col *result) {
  return rect_detect_collision_goal(first, second, result, first->x, first->y);
}

static bool rect_is_different(bump_rect *first, bump_rect *second) {
  return first->x != second->x || first->y != second->y ||
    first->w != second->w || first->h != second->h;
}

/* ------------------------------------------ */
/* -- Grid functions                          */
/* ------------------------------------------ */

static void grid_to_world(double cellsize, bump_coord cx, bump_coord cy, double *x, double *y) {
  *x = (cx-1)*cellsize;
  *y = (cy-1)*cellsize;
}

static void grid_to_cell(double cellsize, double x, double y, bump_coord *cx, bump_coord *cy) {
  *cx = floor(x / cellsize) + 1;
  *cy = floor(y / cellsize) + 1;
}

/*
 * grid_traverse* functions are based on "A Fast Voxel Traversal
 * Algorithm for Ray Tracing", by John Amanides and Andrew Woo -
 * http://www.cse.yorku.ca/~amana/research/grid.pdf
 *
 * It has been modified to include both cells when the ray "touches a
 * grid corner", and with a different exit condition
 */

static void grid_traverse_init_step(double cellSize, double ct, double t1, double t2,
                                    bump_coord *step, double *d, double *t) {
  double v = t2 - t1;

  if (v > 0.0) {
    *step = 1;
    *d = cellSize / v;
    *t = ((ct + v) * cellSize - t1) / v;
  } else if (v < 0.0) {
    *step = -1;
    *d = -cellSize / v;
    *t = ((ct + v - 1) * cellSize - t1) / v;
  } else {
    *step = 0;
    *d = BUMP_HUGE;
    *t = BUMP_HUGE;
  }
}

/* i.e. callback(closure, cx, cy) */
typedef void (*grid_traverse_callback)(void*, bump_coord, bump_coord);

static void grid_traverse(double cellSize, double x1, double y1, double x2, double y2,
                          grid_traverse_callback cb, void* closure) {
  bump_coord cx1, cy1, cx2, cy2, stepX, stepY;
  double dx, tx, dy, ty;
  grid_to_cell(cellSize, x1, y1, &cx1, &cy1);
  grid_to_cell(cellSize, x2, y2, &cx2, &cy2);
  grid_traverse_init_step(cellSize, cx1, x1, x2, &stepX, &dx, &tx);
  grid_traverse_init_step(cellSize, cy1, y1, y2, &stepY, &dy, &ty);
  bump_coord cx = cx1;
  bump_coord cy = cx2;

  cb(closure, cx, cy);

  /* The default implementation had an infinite loop problem when
   * approaching the last cell in some occassions. We finish iterating
   * when we are *next* to the last cell */
  while (abs(cx - cx2) + abs(cy - cy2) > 1) {
    if (tx < ty) {
      tx += dx;
      cx += stepX;
      cb(closure, cx, cy);
    } else {
      /* Addition: include both cells when going through corners */
      if (tx == ty) {
        cb(closure, cx + stepX, cy);
      }
      ty += dy;
      cy += stepY;
      cb(closure, cx, cy);
    }
  }

  /* If we have not arrived to the last cell, use it */
  if (cx != cx2 || cy != cy2) {
    cb(closure, cx, cy);
  }
}

static void grid_to_cell_rect(double cellSize, bump_rect *input, bump_rect *output) {
  bump_coord cx, cy;
  grid_to_cell(cellSize, input->x, input->y, &cx, &cy);
  double cr = ceil((input->x+input->w) / cellSize);
  double cb = ceil((input->y+input->h) / cellSize);

  output->x = cx;
  output->y = cy;
  output->w = cr - cx + 1;
  output->h = cb - cy + 1;
}

/* ------------------------------------------ */
/* -- Responses                               */
/* ------------------------------------------ */

typedef void(bump_response_func)(bump_world*, bump_col *, bump_rect *,
                                 double, double, bump_filter_func *, void*,
                                 double *, double *, bump_col_list *);

static void bump_nil(bump_world *w, bump_col *col, bump_rect *r, double goalX,
                     double goalY, bump_filter_func *filter, void *closure, double *rx,
                     double *ry, bump_col_list *cols) {
  col->response = r_nil;
}

static void bump_touch(bump_world *w, bump_col *col, bump_rect *r,
                       double goalX, double goalY, bump_filter_func *filter,
                       void *closure, double *rx, double *ry, bump_col_list *cols) {
  col->response = r_touch;
  *rx = col->touch_x;
  *ry = col->touch_y;
}

static void bump_cross(bump_world *w, bump_col *col, bump_rect *r,
                       double goalX, double goalY, bump_filter_func *filter,
                       void *closure, double *rx, double *ry, bump_col_list *cols) {
  col->response = r_cross;
  bump_world_project(w, col->item, r, goalX, goalY, filter, closure, cols);
  *rx = goalX;
  *ry = goalY;
}

static void bump_slide(bump_world *w, bump_col *col, bump_rect *r,
                       double goalX, double goalY, bump_filter_func *filter,
                       void *closure, double *rx, double *ry, bump_col_list *cols) {
  /*
   * GoalX/Y optional in Lua:
   *
   * goalX = goalX or x
   * goalY = goalY or y
   */
  col->response = r_slide;
  double tch_x, tch_y, gx, gy;
  gx = goalX;
  gy = goalY;
  tch_x = col->touch_x;
  tch_y = col->touch_y;

  if (col->move_x != 0.0 || col->move_y != 0.0) {
    if (col->normal_x != 0.0) {
      gx = tch_x;
    } else {
      gy = tch_y;
    }
  }

  col->response_x = gx;
  col->response_y = gy;

  internal_rect.x = tch_x;
  internal_rect.y = tch_y;
  internal_rect.h = r->h;
  internal_rect.w = r->w;

  bump_world_project(w, col->item, &internal_rect, gx, gy, filter, closure, cols);
  *rx = gx;
  *ry = gy;
}

static void bump_bounce(bump_world *w, bump_col *col, bump_rect *r,
                        double goalX, double goalY, bump_filter_func *filter,
                        void *closure, double *rx, double *ry, bump_col_list *cols) {
  /*
   * GoalX/Y optional in Lua:
   *
   * goalX = goalX or x
   * goalY = goalY or y
   */
  double tx, ty, gx, gy, bx, by;
  gx = goalX;
  gy = goalY;
  tx = col->touch_x;
  ty = col->touch_y;
  bx = tx;
  by = ty;

  if (col->move_x != 0.0 || col->move_y != 0.0) {
    double bnx = gx - tx;
    double bny = gy - ty;
    if (col->normal_x == 0.0) {
      bny = -bny;
    } else {
      bnx = -bnx;
    }
    bx = tx + bnx;
    by = ty + bny;
  }

  col->response_x = bx;
  col->response_y = by;

  internal_rect.x = tx;
  internal_rect.y = ty;
  internal_rect.h = r->h;
  internal_rect.w = r->w;
  gx = bx;
  gy = by;

  bump_world_project(w, col->item, &internal_rect, gx, gy, filter, closure, cols);
  *rx = gx;
  *ry = gy;
}

/* ------------------------------------------ */
/* -- World                                   */
/* ------------------------------------------ */

/* Private functions and methods */

typedef struct {
  bump_item *item;
  bump_item *other;
  double ti;
  double ti1;
  double ti2;
  double weight;
} bump_item_info;

bool sort_by_weight(bump_item_info *a, bump_item_info *b) {
  return a->weight < b->weight;
}

bool sort_by_ti_and_distance_ii(bump_item_info *a, bump_item_info *b) {
  if (a->ti == b->ti) {
    double ad = rect_get_square_distance(&a->item->rect, &a->other->rect);
    double bd = rect_get_square_distance(&a->item->rect, &b->other->rect);
    return ad < bd;
  }
  return a->ti < b->ti;
}

int sort_by_ti_and_distance_col(const void *l, const void *r) {
  bump_col *a = *(bump_col**) l;
  bump_col *b = *(bump_col**) r;

  if (a->ti == b->ti) {
    double ad = rect_get_square_distance(&a->item->rect, &a->other->rect);
    double bd = rect_get_square_distance(&a->item->rect, &b->other->rect);
    return ad < bd;
  }
  return a->ti < b->ti;
}

static void world_add_item_to_cell(bump_world *self, bump_item *item, bump_coord cx,
                                   bump_coord cy) {
  bump_point point = POINT(cx, cy);

  bump_cell *cell = bump_world_get_cell(self, point);

  if (cell == NULL) {
    cell = cell_new(point, BUMP_CELL_SIZE);
    HASH_ADD(hh, self->map, point, sizeof(bump_point), cell);
  }

  cell_push_item(cell, item);
}

static bool world_remove_item_from_cell(bump_world *self, bump_item *item, bump_coord cx,
                                        bump_coord cy) {
  bump_point point = POINT(cx, cy);

  bump_cell *cell = bump_world_get_cell(self, point);

  if (cell == NULL) {
    return false;
  }

  return cell_remove_item_by_key(cell, item->key);
}

/*
 * In all cases this replaces getDictItemsInCellRect + subsequent for
 * pairs loops.
 *
 * What we really want to do when we call that function is perform
 * some action on every item in a given cell rect.  It's much better
 * to do this with a callback/lambda/closure/whatever than by
 * allocating a list you're only going to use once.
 */

typedef void(items_in_cell_rect_cb)(bump_world *, bump_item *,
                                    bump_coord, bump_coord, void*);

static void world_for_items_in_cell_rect(bump_world *self, bump_rect *rect,
                                         items_in_cell_rect_cb cb, void *closure) {
  /* The "visited" number should* be different for every call of this
     function. u64 is very big. */
  static u64 visited = VISITED_1;
  /*
   * But just in case, let's use a "visit time" as well. Combined this
   * will mean:
   *
   * 1. subsequent calls in short time should trip visited
   * 2. subsequent calls in long time (i.e. after visited overflows)
   *    should trip visit_time
   * 3. If I'm stupid and malloc items instead of callocing them, it'll
   *    be rare that both of these match by accident.
   */
  time_t visit_time = time(NULL);

  bump_coord cl, ct, cw, ch;
  cl = (bump_coord)rect->x;
  ct = (bump_coord)rect->y;
  cw = (bump_coord)rect->w;
  ch = (bump_coord)rect->h;

  for (bump_coord cy = ct; cy < ct+ch; cy++) {
    for (bump_coord cx = cl; cx < cl+cw; cx++) {
      bump_point point = POINT(cx, cy);
      bump_cell *cell = bump_world_get_cell(self, point);

      if (cell != NULL) {
        for (size_t i = 0; i < cell->items_len; i++) {
          bump_item *item = cell->items[i];

          if (!(item->visited_1 == visited && item->visit_time_1 == visit_time)) {
            /* i.e. don't visit an item more than once */
            item->visited_1 = visited;
            item->visit_time_1 = visit_time;

            cb(self, item, cx, cy, closure);
          }
        }
      }
    }
  }

  /* Since visited is static, it'll be different next time. */
  visited = UPDATE_VISITED(VISITED_1, visited);
}

static bump_response_func *get_response_func(bump_response r) {
  switch (r) {
  case r_touch:
    return bump_touch;
  case r_cross:
    return bump_cross;
  case r_slide:
    return bump_slide;
  case r_bounce:
    return bump_bounce;
  default:
    return bump_nil;
  }
}

/*
 * Custom response logic - not currently supported.

-- misc public methods
function World:addResponse(name, response)
  self.responses[name] = response
end
 */

/* Misc Public Methods */

bump_response default_filter(bump_item *item, bump_item *other, void *closure) {
  return r_slide;
}

struct world_project_closure_s {
  bump_item *item;
  bump_rect *r;
  double goalX;
  double goalY;
  bump_filter_func *filter;
  void *filter_closure;
  bump_col_list *cols;
};

static void project_cb(bump_world *self, bump_item *other,
                       bump_coord goalX, bump_coord goalY, void* closure) {
  struct world_project_closure_s *wpc =
    (struct world_project_closure_s*)closure;

  bump_response response = wpc->filter(wpc->item, other, wpc->filter_closure);

  if (response != r_nil) {
    bump_col col;

    if (rect_detect_collision_goal(wpc->r, &other->rect, &col, goalX, goalY)) {
      col.response = response;
      col.item = wpc->item;
      col.other = other;

      bump_col_list_append(wpc->cols, &col);
    }
  }
}

void bump_world_project(bump_world *self, bump_item *item, bump_rect *r,
                        double goalX, double goalY, bump_filter_func *filter,
                        void *closure, bump_col_list *cols) {
  if (filter == NULL) {
    filter = default_filter;
  }

  struct world_project_closure_s wpc;
  wpc.item = item;
  wpc.r = r;
  wpc.goalX = goalX;
  wpc.goalY = goalY;
  wpc.filter = filter;
  wpc.filter_closure = closure;
  wpc.cols = cols;

  /*
   * This could probably be done with less cells using a polygon
   * raster over the cells instead of a bounding rect of the whole
   * movement.  Conditional to building a queryPolygon method
   */
  project_t_rect.x = min(goalX, r->x);
  project_t_rect.y = min(goalY, r->y);
  project_t_rect.w = max(goalX+r->w, r->x + r->w) - project_t_rect.x;
  project_t_rect.h = max(goalY+r->h, r->y + r->h) - project_t_rect.y;

  grid_to_cell_rect(self->cellsize, &project_t_rect, &project_cell_rect);

  world_for_items_in_cell_rect(self, &project_cell_rect, project_cb, &wpc);

  if (cols->len > 0) {
    qsort(cols->data, cols->len, sizeof(bump_col*),
          sort_by_ti_and_distance_col);
  }
}

size_t bump_world_count_cells(bump_world *self) {
  /*
function World:countCells()
  local count = 0
  for _,row in pairs(self.rows) do
    for _,_ in pairs(row) do
      count = count + 1
    end
  end
  return count
end
  */
  return 0;
}

bool bump_world_has_item(bump_world *self, bump_key item_key) {
  bump_item *item = bump_world_get_item(self, item_key);

  return item != NULL;
}

/* Unsure these are needed.
function World:getItems()
  local items, len = {}, 0
  for item,_ in pairs(self.rects) do
    len = len + 1
    items[len] = item
  end
  return items, len
end

function World:countItems()
  local len = 0
  for _ in pairs(self.rects) do len = len + 1 end
  return len
end
*/

bool bump_world_get_rect(bump_world *self, bump_key item_key, bump_rect *output) {
  bump_item *item = bump_world_get_item(self, item_key);

  if (item == NULL) {
    return false;
  }

  memcpy(output, &item->rect, sizeof(bump_rect));

  return true;
}

void bump_world_to_world(bump_world *self, bump_coord cx, bump_coord cy,
                         double *x, double *y) {
  grid_to_world(self->cellsize, cx, cy, x, y);
}

void bump_world_to_cell(bump_world *self, double x, double y,
                   bump_coord *cx, bump_coord *cy) {
  grid_to_cell(self->cellsize, x, y, cx, cy);
}

/* Query methods */
/* TODO: these methods :^) */

/* Main methods */

bump_item *bump_world_add_rect(bump_world *self, bump_key key,
                               bump_item_data data, bump_rect *r) {
  if (bump_world_has_item(self, key)) {
    /* Already added */
    return NULL;
  }

  bump_item *newItem = calloc(1, sizeof(*newItem));
  newItem->key = key;
  newItem->data = data;
  memcpy(&newItem->rect, r, sizeof(bump_rect));

  HASH_ADD(hh, self->items, key, sizeof(bump_item), newItem);

  grid_to_cell_rect(self->cellsize, r, &project_cell_rect);

  bump_coord cl, ct, cw, ch;
  cl = (bump_coord)project_cell_rect.x;
  ct = (bump_coord)project_cell_rect.y;
  cw = (bump_coord)project_cell_rect.w;
  ch = (bump_coord)project_cell_rect.h;

  for (bump_coord cy = ct; cy < ct+ch; cy++) {
    for (bump_coord cx = cl; cx < cl+cw; cx++) {
      world_add_item_to_cell(self, newItem, cx, cy);
    }
  }

  return newItem;
}

bump_item *bump_world_add(bump_world *self, bump_key key, bump_item_data data,
                          double x, double y, double w, double h) {
  internal_rect.x = x;
  internal_rect.y = y;
  internal_rect.h = h;
  internal_rect.w = w;
  return bump_world_add_rect(self, key, data, &internal_rect);
}

bump_item_data bump_world_remove_item(bump_world *self, bump_item *item) {
  if (item == NULL) {
    return BUMP_ITEM_DATA_NONE;
  }

  grid_to_cell_rect(self->cellsize, &item->rect, &project_cell_rect);

  bump_coord cl, ct, cw, ch;
  cl = (bump_coord)project_cell_rect.x;
  ct = (bump_coord)project_cell_rect.y;
  cw = (bump_coord)project_cell_rect.w;
  ch = (bump_coord)project_cell_rect.h;

  for (bump_coord cy = ct; cy < ct+ch; cy++) {
    for (bump_coord cx = cl; cx < cl+cw; cx++) {
      bump_cell *cell = bump_world_get_cell(self, POINT(cx, cy));
      if (cell != NULL) {
        cell_remove_item_by_key(cell, item->key);
      }
    }
  }

  bump_item_data ret = item->data;
  HASH_DEL(self->items, item);
  free(item);
  return ret;
}

bump_item_data bump_world_remove(bump_world *self, bump_key key) {
  return bump_world_remove_item(self, bump_world_get_item(self, key));
}

void bump_world_update_item_rect(bump_world *self, bump_item *item,
                                 bump_rect *dest) {
  if (item == NULL) {
    return;
  }

  if (rect_is_different(dest, &item->rect)) {
    grid_to_cell_rect(self->cellsize, &item->rect, &update_c1_rect);
    grid_to_cell_rect(self->cellsize, dest, &update_c2_rect);

    /* If the cell rects are different,
       we've moved cells & need to update them */
    if (rect_is_different(&update_c1_rect, &update_c2_rect)) {
      bump_coord cr1 = update_c1_rect.x + update_c1_rect.w - 1;
      bump_coord cb1 = update_c1_rect.y + update_c1_rect.h - 1;
      bump_coord cr2 = update_c2_rect.x + update_c2_rect.w - 1;
      bump_coord cb2 = update_c2_rect.y + update_c2_rect.h - 1;
      bool cyOut;

      for (bump_coord cy = update_c1_rect.y; cy <= cb1; cy++) {
        cyOut = cy < update_c2_rect.y || cy > cb2;
        for (bump_coord cx = update_c1_rect.x; cx <= cr1; cx++) {
          if (cyOut || cx < update_c2_rect.x || cx > cr2) {
            cell_remove_item_by_key(bump_world_get_cell(self, POINT(cx, cy)),
                                    item->key);
          }
        }
      }

      for (bump_coord cy = update_c2_rect.y; cy <= cb2; cy++) {
        cyOut = cy < update_c1_rect.y || cy > cb1;
        for (bump_coord cx = update_c2_rect.x; cx <= cr2; cx++) {
          if (cyOut || cx < update_c1_rect.x || cx > cr1) {
            cell_push_item(bump_world_get_cell(self, POINT(cx, cy)), item);
          }
        }
      }
    }

    memcpy(&item->rect, dest, sizeof(bump_rect));
  }
}

void bump_world_update_item(bump_world *self, bump_item *item,
                            double x, double y, double w, double h) {
  update_rect.x = x;
  update_rect.y = y;
  update_rect.h = h;
  update_rect.w = w;
  bump_world_update_item_rect(self, item, &update_rect);
}

/* Updates the position & size of the item with given key to the rectangle
   dest */
void bump_world_update_rect(bump_world *self, bump_key key, bump_rect *dest) {
  return bump_world_update_item_rect(self, bump_world_get_item(self, key),
                                     dest);
}

/* Updates the position & size of the item with given key to the
   co-ordinates & width/height given */
void bump_world_update(bump_world *self, bump_key key,
                       double x, double y, double w, double h) {
  return bump_world_update_item(self, bump_world_get_item(self, key),
                                x, y, w, h);
}

typedef struct {
  u64 visited;
  time_t visit_time;
  bump_filter_func *filter;
  void *closure;
} vf_closure;

bump_response visited_filter(bump_item* item, bump_item* other, void* closure) {
  vf_closure *vf = (vf_closure*)closure;

  if (other->visited_2 == vf->visited && other->visit_time_2 == vf->visit_time) {
    return r_nil;
  }

  return vf->filter(item, other, vf->closure);
}

void bump_world_check_item(bump_world *self, bump_item *item,
                           double goalX, double goalY,
                           double *actualX, double *actualY,
                           bump_filter_func *filter, void *closure,
                           bump_col_list *cols) {
  if (item == NULL) {
    return;
  }

  static u64 visited = VISITED_2;
  time_t visit_time = time(NULL);
  *actualX = goalX;
  *actualY = goalY;

  item->visited_2 = visited;
  item->visit_time_2 = visit_time;

  bump_col_list_clear(cols);
  bump_col_list_clear(self->util_projected);

  vf_closure vfc;
  vfc.visited = visited;
  vfc.visit_time = visit_time;
  vfc.filter = filter;
  vfc.closure = closure;

  memcpy(&update_rect, &item->rect, sizeof(bump_rect));

  bump_world_project(self, item, &update_rect, goalX, goalY,
                     visited_filter, &vfc, self->util_projected);

  while (self->util_projected->len > 0) {
    bump_col *col = &self->util_projected->data[0];

    bump_col_list_append(cols, col);
    col = &cols->data[0];

    col->other->visited_2 = visited;
    col->other->visit_time_2 = visit_time;

    bump_response_func *response = get_response_func(col->response);

    bump_col_list_clear(self->util_projected);
    response(self, col, &update_rect, *actualX, *actualY, visited_filter, &vfc,
             actualX, actualY, self->util_projected);
  }

  visited = UPDATE_VISITED(VISITED_2, visited);
}

void bump_world_check(bump_world *self, bump_key key, double goalX,
                      double goalY, double *actualX, double *actualY,
                      bump_filter_func *filter, void *closure,
                      bump_col_list *cols) {
  bump_world_check_item(self, bump_world_get_item(self, key), goalX, goalY,
                        actualX, actualY, filter, closure, cols);
}

void bump_world_move_item(bump_world *self, bump_item *item, double goalX,
                          double goalY, double *actualX, double *actualY,
                          bump_filter_func *filter, void *closure,
                          bump_col_list *cols) {
  if (item == NULL) {
    return;
  }

  bump_world_check_item(self, item, goalX, goalY, actualX, actualY, filter,
                        closure, cols);

  bump_world_update_item(self, item, *actualX, *actualY, item->rect.w,
                         item->rect.h);
}

void bump_world_move(bump_world *self, bump_key key, double goalX,
                     double goalY, double *actualX, double *actualY,
                     bump_filter_func *filter, void *closure,
                     bump_col_list *cols) {
  bump_world_move_item(self, bump_world_get_item(self, key), goalX, goalY,
                       actualX, actualY, filter, closure, cols);
}

/*

------------------------------------------
-- World
------------------------------------------

--- Query methods

function World:queryRect(x,y,w,h, filter)

  assertIsRect(x,y,w,h)

  local cl,ct,cw,ch = grid_toCellRect(self.cellSize, x,y,w,h)
  local dictItemsInCellRect = getDictItemsInCellRect(self, cl,ct,cw,ch)

  local items, len = {}, 0

  local rect
  for item,_ in pairs(dictItemsInCellRect) do
    rect = self.rects[item]
    if (not filter or filter(item))
    and rect_isIntersecting(x,y,w,h, rect.x, rect.y, rect.w, rect.h)
    then
      len = len + 1
      items[len] = item
    end
  end

  return items, len
end

function World:queryPoint(x,y, filter)
  local cx,cy = self:toCell(x,y)
  local dictItemsInCellRect = getDictItemsInCellRect(self, cx,cy,1,1)

  local items, len = {}, 0

  local rect
  for item,_ in pairs(dictItemsInCellRect) do
    rect = self.rects[item]
    if (not filter or filter(item))
    and rect_containsPoint(rect.x, rect.y, rect.w, rect.h, x, y)
    then
      len = len + 1
      items[len] = item
    end
  end

  return items, len
end

function World:querySegment(x1, y1, x2, y2, filter)
  local itemInfo, len = getInfoAboutItemsTouchedBySegment(self, x1, y1, x2, y2, filter)
  local items = {}
  for i=1, len do
    items[i] = itemInfo[i].item
  end
  return items, len
end

function World:querySegmentWithCoords(x1, y1, x2, y2, filter)
  local itemInfo, len = getInfoAboutItemsTouchedBySegment(self, x1, y1, x2, y2, filter)
  local dx, dy        = x2-x1, y2-y1
  local info, ti1, ti2
  for i=1, len do
    info  = itemInfo[i]
    ti1   = info.ti1
    ti2   = info.ti2

    info.weight  = nil
    info.x1      = x1 + dx * ti1
    info.y1      = y1 + dy * ti1
    info.x2      = x1 + dx * ti2
    info.y2      = y1 + dy * ti2
  end
  return itemInfo, len
end

*/
