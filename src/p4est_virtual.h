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

/** \file p4est_virtual.h
 *
 * Layer of virtual quadrants in coarse quadrants at refinement boundaries for
 * 2:1 balanced trees in 2D.
 *
 * \ingroup p4est
 */

#ifndef P4EST_VIRTUAL_H
#define P4EST_VIRTUAL_H

#include <p4est.h>

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
typedef struct
{
  p4est_connect_type_t btype;           /**< which neighbors are considered in
                                             virtual */

  int8_t             *virtual_qflags;   /**< stores if quad has virtual quads */
  int8_t             *virtual_gflags;   /**< stores if ghost has virtual
                                             quads */
  sc_array_t         *virtual_qlevels;  /**< stores lists of per-level
                                             virtual quads, NULL by default,
                                             can be enabled in
                                             \ref p4est_virtual_new_ext */
  sc_array_t         *virtual_glevels;  /**< stores lists of per-level
                                             virtual ghosts, NULL by default,
                                             can be enabled in
                                             \ref p4est_virtual_new_ext */
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
} p4est_virtual_t;

/* -------------------------------------------------------------------------- */
/* |                             Ghost exchange                             | */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* |                            Neighbor search                             | */
/* -------------------------------------------------------------------------- */

/* *INDENT-OFF* */
/** Store for each virtual subquad the index in Morton-order of the neighboring
 * virtual subquad in the direction of a face index. Note, that we are only
 * interested in external neighbors of the same size as the virtual quadrants,
 * i.e. the neighboring quad needs to be hanging. For choosing the appropriate
 * subquad, we have to perform a face corner transformation.
 * For internal neighbors we store the respective virtual quadrant index 0 .. 3,
 * external face neighbors are assigned values s = 4 .. 11 with
 *   s = 4 + 4 * k + f
 * where k is the index of subquad (0 <= k < 2) of the local face we have to
 * find the adjacent index in the opposite tree in and f the index of the local
 * face we have to query.
 */
extern const int    p4est_face_virtual_neighbors_inside[P4EST_CHILDREN]
                                                       [P4EST_FACES];

/** Store for each virtual subquad the index in Morton-order of the neighboring
 * virtual subquad in the direction of a corner index. Note, that we are only
 * interested in external neighbors of the same size as the virtual quadrants,
 * i.e. the neighboring quad needs to be hanging.
 * Additionally note, that it may be necessary to query face neighbors to find
 * the respective virtual quadrant, e.g. when querying the neighbor of subquad
 * v0 in direction of corner c1.
 * For choosing the appropriate subquad, we have to perform a face corner
 * transformation, depending on which neighbor has been queried.
 * For internal neighbors we store the respective virtual quadrant index 0 .. 3,
 * external face neighbors are assigned values s = 4 .. 11 with
 *   s = 4 + 4 * k + f
 * where k is the index of the subquad (0 <= k < 2) of the local face we have to
 * find the adjacent index in the opposite tree in and f the index of the local
 * face we have to query.
 * External corner neighbors are assigned values s = 12 .. 15 with
 *   s = 12 + c
 * where c is the respective corner index.
 */
extern const int    p4est_corner_virtual_neighbors_inside[P4EST_CHILDREN]
                                                         [P4EST_CHILDREN];
/* *INDENT-ON* */

SC_EXTERN_C_END;

#endif /* P4EST_VIRTUAL_H */
