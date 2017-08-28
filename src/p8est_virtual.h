/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file p8est_virtual.h
 *
 * Layer of virtual quadrants in coarse quadrants at refinement boundaries for
 * 2:1 balanced trees in 3D.
 *
 * \ingroup p8est
 */

#ifndef P8EST_VIRTUAL_H
#define P8EST_VIRTUAL_H

#include <p8est.h>

SC_EXTERN_C_BEGIN;

/** Struct for storing virtual quadrants at refinement boundaries in coarse
 * quadrants for algorithms that rely on interaction to only take place between
 * same-sized quadrants.
 *
 * btype stores which quadrants we consider to decide if quadrant contains
 * virtual quadrants.
 *
 * If a quadrant contains virtual quadrants is stored in virtual_qflags for
 * local quadrants and virtual_gflags for ghost quadrants.
 *
 * For each level, an array of local quadrant numbers is stored in
 * virtual_qlevels if the respective quadrant contains virtual quadrants. Note,
 * that the respective index is inserted only once and added to the list of
 * quadrants on the level of the virtual quadrants.
 * For ghost quadrants the same information is stored in virtual_glevels.
 * Both arrays (virtual_qlevels and virtual_glevels) array are NULL by default
 * and can be enabled using \ref p4est_virtual_new_ext.
 *
 * For each local quadrant we store how many quadrants, real and virtual,
 * precede it on its level in Morton ordering in quad_qreal_offset. The same
 * information is stored for the first virtual quadrant if the quadrant contains
 * virtual quadrants in quad_qvirtual_offset. If the current quadrant does not
 * contain virtual quadrants, -1 is stored.
 */
typedef struct {
  p8est_connect_type_t btype;           /**< which neighbors are considered in
                                             virtual */

  int8_t             *virtual_qflags;   /**< stores if quad has virtual quads */
  int8_t             *virtual_gflags;   /**< stores if ghost has virtual
                                             quads */
  sc_array_t         *virtual_qlevels;  /**< stores lists of per-level
                                             virtual quads, NULL by default */
  sc_array_t         *virtual_glevels;  /**< stores lists of per-level
                                             virtual ghosts, NULL by default */
  p4est_locidx_t     *quad_qreal_offset;        /**< stores number of quadrants
                                                     preceding current quadrant
                                                     on its level, either real
                                                     or virtual. */
  p4est_locidx_t     *quad_qvirtual_offset;     /**< stores number of quadrants
                                                     preceding first virtual
                                                     sub-quadrant of current
                                                     quadrant on the level of
                                                     its virtual sub-quads if it
                                                     has any, either real or
                                                     virtual.  Else -1. */
  p4est_locidx_t     *quad_greal_offset;        /**< stores number of quadrants
                                                     preceeding current ghost
                                                     on its level, either real
                                                     or virtual. */
  p4est_locidx_t     *quad_gvirtual_offset;     /**< stores number of quadrants
                                                     preceding first virtual
                                                     sub-quadrant of current
                                                     ghost on the level of its
                                                     virtual sub-quads if it has
                                                     any, either real or
                                                     virtual.  Else -1. */
} p8est_virtual_t;

/* *INDENT-OFF* */
/** Store for each virtual subquad the index in Morton-order of the neighboring
 * virtual subquad in the direction of a face index. Note, that we are only
 * interested in external neighbors of the same size as the virtual quadrants,
 * i.e. the neighboring quad needs to be hanging. For choosing the appropriate
 * subquad, we have to perform a face corner transformation.
 * For internal neighbors we store the respective virtual quadrant index 0 .. 7,
 * external face neighbors are assigned values s = 8 .. 31 with
 *   s = 8 + 6 * k + f
 * where k is the index of subquad (0 <= k < 4) of the local face we have to
 * find the adjacent index in the opposite tree in and f the index of the local
 * face we have to query.
 */
extern const int    p8est_face_virtual_neighbors_inside[P8EST_CHILDREN]
                                                       [P8EST_FACES];

/** Store for each virtual subquad the index in Morton-order of the neighboring
 * virtual subquad in the direction of an edge index. Note, that we are only
 * interested in external neighbors of the same size as the virtual quadrants,
 * i.e. the neighboring quad needs to be hanging.
 * Additionally note, that it may be necessary to query face neighbors to find
 * the respective virtual quadrant, e.g. when querying the neighbor of subquad
 * v0 in direction of edge e1.
 * For choosing the appropriate subquad, we have to perform a face corner or
 * edge corner transformation, depending on which neighbor has been queried.
 * For internal neighbors we store the respective virtual quadrant index 0 .. 7,
 * external face neighbors are assigned values s = 8 .. 31 with
 *   s = 8 + 6 * k + f
 * where k is the index of the subquad (0 <= k < 4) of the local face we have to
 * find the adjacent index in the opposite tree in and f the index of the local
 * face we have to query.
 * External edge neighbors are assigned values s = 32 .. 55 with
 *   s = 32 + 12 * k + e
 * where k is the index of the subquad (0 <= k < 2) of the local edge we have to
 * find the adjacent index in the opposite tree in and e the index of the local
 * edge we have to query.
 */
extern const int    p8est_edge_virtual_neighbors_inside[P8EST_CHILDREN]
                                                       [P8EST_EDGES];

/** Store for each virtual subquad the index in Morton-order of the neighboring
 * virtual subquad in the direction of a corner index. Note, that we are only
 * interested in external neighbors of the same size as the virtual quadrants,
 * i.e. the neighboring quad needs to be hanging.
 * Additionally note, that it may be necessary to query face or edge neighbors
 * to find the respective virtual quadrant, e.g. when querying the neighbor of
 * subquad v0 in direction of corner c1.
 * For choosing the appropriate subquad, we may have to perform a face or edge
 * corner transformation, depending on which neighbor has been queried.
 * For internal neighbors we store the respective virtual quadrant index 0 .. 7,
 * external face neighbors are assigned values s = 8 .. 31 with
 *   s = 8 + 6 * k + f
 * where k is the index of subquad (0 <= k < 4) of the local face we have to
 * find the adjacent index in the opposite tree in and f the index of the local
 * face we have to query.
 * External edge neighbors are assigned values s = 32 .. 55 with
 *   s = 32 + 12 * k + e
 * where k is the index of the subquad (0 <= k < 2) of the local edge we have to
 * find the adjacent index in the opposite tree in and e the index of the local
 * edge we have to query.
 * External corner neighbors are assigned values s = 56 .. 63 with
 *   s = 56 + c
 * where c is the respective corner index.
 */
extern const int    p8est_corner_virtual_neighbors_inside[P8EST_CHILDREN]
                                                         [P8EST_CHILDREN];
/* *INDENT-ON* */

SC_EXTERN_C_END;

#endif /* P8EST_VIRTUAL_H */
