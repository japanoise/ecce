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
#include "int.h"

/* You almost definitely want to change this */
typedef void *bump_item_data;
#define BUMP_ITEM_DATA_NONE NULL

/* You may want to change this. Or not. I'm not a cop. */
typedef i64 bump_key;

typedef enum {
  r_nil,
  r_touch,
  r_cross,
  r_slide,
  r_bounce,
} bump_response;

typedef struct {
  double x;
  double y;
  double w;
  double h;
} bump_rect;

struct bump_item_s;
typedef struct bump_item_s bump_item;

bump_key bump_item_get_key(bump_item*);
bump_item_data bump_item_get_data(bump_item*);
void bump_item_set_data(bump_item*, bump_item_data);
bump_rect *bump_item_get_rect(bump_item*);

typedef struct {
  bool overlaps;
  double ti;
  double move_x;
  double move_y;
  double normal_x;
  double normal_y;
  double touch_x;
  double touch_y;
  bump_response response;
  double response_x;
  double response_y;
  bump_item *item;
  bump_item *other;
} bump_col;

struct bump_world_s;
typedef struct bump_world_s bump_world;

bump_world *bump_world_create();

bump_world *bump_world_create_cellsize(double cellsize);

void bump_world_destroy(bump_world *w);

/*
 * In Lua, you can just merrily return lists from functions all day
 * long without care as to how they're allocated or who's going to
 * manage them or how they're going to be used.  In C we have to be a
 * little more careful.
 *
 * bump_col_list is a list of collisions.  Specifically, it's a
 * "stretchy" array with an explicit length.  The idea is that you'll
 * only allocate as many of these as you actually need (i.e. one,
 * unless you need to go through multiple collision lists at once),
 * and they'll grow only about as big as they think they need to be.
 */

typedef struct {
  bump_col *data;
  size_t len;

  /* Leave this alone, or you'll mess something up. */
  size_t cap;
} bump_col_list;

bump_col_list *bump_col_list_create(size_t initialCapacity);

/* Copies the data */
void bump_col_list_append(bump_col_list *l, bump_col *data);

/* "Clears" the list without deallocating anything.  This is much
   faster than destroying and recreating */
void bump_col_list_clear(bump_col_list *l);

/* Frees (deallocates) the list */
void bump_col_list_destroy(bump_col_list *l);

/* A very large number - used for some functions */
extern const double BUMP_HUGE;

/* Main methods */

/* Retrieves a previously added item */
bump_item *bump_world_get_item(bump_world *self, bump_key item_key);

/* Returns whether or not the world has the item. Functionally the
   same as bump_world_get_item(self, item_key) != NULL */
bool bump_world_has_item(bump_world *self, bump_key item_key);

/* Adds an item to a world, given a rectangle to put it in. Returns
   NULL if an item with that key has already been added, the new item
   otherwise. */
bump_item *bump_world_add_rect(bump_world *world, bump_key key,
                               bump_item_data data, bump_rect *r);

/* Adds an item to a world, given an x, y position & width, height
   size to put it in. Returns NULL if an item with that key has
   already been added, the new item otherwise. */
bump_item *bump_world_add(bump_world *world, bump_key key, bump_item_data data,
                          double x, double y, double w, double h);

/* Removes the item from the world. Returns the item's data - it does
   not free that data itself. */
bump_item_data bump_world_remove_item(bump_world *self, bump_item *item);

/* Removes the item with given key from the world. Returns the item's
   data - it does not free that data itself. */
bump_item_data bump_world_remove(bump_world *self, bump_key key);

/* Updates the position & size of the given item to the rectangle
   dest */
void bump_world_update_item_rect(bump_world *self, bump_item *item,
                                 bump_rect *dest);

/* Updates the position & size of the given item to the co-ordinates &
   width/height given */
void bump_world_update_item(bump_world *self, bump_item *item, double x,
                            double y, double w, double h);

/* Updates the position & size of the item with given key to the rectangle
   dest */
void bump_world_update_rect(bump_world *self, bump_key key, bump_rect *dest);

/* Updates the position & size of the item with given key to the
   co-ordinates & width/height given */
void bump_world_update(bump_world *self, bump_key key, double x, double y,
                       double w, double h);

/* Given two items (and a closure), return the type of collision
   response to use */
typedef bump_response(bump_filter_func)(bump_item*, bump_item*, void*);

/* default filter (returns slide every time) */
bump_response default_filter(bump_item *item, bump_item *other, void *closure);

/* Checks for collisions in item's path to (goalX,goalY); puts those
   collisions into cols, and the actual (x,y) position of item into
   (actualX, actualY). */
void bump_world_check_item(bump_world *self, bump_item *item, double goalX,
                           double goalY, double *actualX, double *actualY,
                           bump_filter_func *filter, void *closure,
                           bump_col_list *cols);

/* Checks for collisions in the item with given key's path to
   (goalX,goalY); puts those collisions into cols, and the actual
   (x,y) position of item into (actualX, actualY). */
void bump_world_check(bump_world *self, bump_key key, double goalX,
                      double goalY, double *actualX, double *actualY,
                      bump_filter_func *filter, void *closure,
                      bump_col_list *cols);

/* Checks for collisions in item's path to (goalX,goalY); puts those
   collisions into cols, and the actual (x,y) position of item into
   (actualX, actualY). Moves the item to (actualX, actualY) */
void bump_world_move_item(bump_world *self, bump_item *item, double goalX,
                           double goalY, double *actualX, double *actualY,
                           bump_filter_func *filter, void *closure,
                           bump_col_list *cols);

/* Checks for collisions in the item with given key's path to
   (goalX,goalY); puts those collisions into cols, and the actual
   (x,y) position of item into (actualX, actualY). Moves the item to
   (actualX, actualY) */
void bump_world_move(bump_world *self, bump_key key, double goalX,
                      double goalY, double *actualX, double *actualY,
                      bump_filter_func *filter, void *closure,
                      bump_col_list *cols);
