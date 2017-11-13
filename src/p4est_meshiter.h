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

#ifndef P4EST_MESHITERATE_H
#define P4EST_MESHITERATE_H

#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_virtual.h>

SC_EXTERN_C_BEGIN;

#define P4EST_MESHITER_DONE 12

/** Specify which quadrants to traverse w.r.t. local and ghost quads
 */
typedef enum
{
  P4EST_TRAVERSE_LOCAL,
  P4EST_TRAVERSE_LOCALGHOST,
  P4EST_TRAVERSE_GHOST
}
p4est_meshiter_localghost_t;

/** Specify which quadrants to traverse w.r.t. real and virtual quads
 */
typedef enum
{
  P4EST_TRAVERSE_REAL,
  P4EST_TRAVERSE_REALVIRTUAL,
  P4EST_TRAVERSE_VIRTUAL
}
p4est_meshiter_realvirt_t;

/** Specify which quadrants to traverse w.r.t. the parallel boundary
 */
typedef enum
{
  P4EST_TRAVERSE_PARBOUNDINNER,
  P4EST_TRAVERSE_PARALLEL_BOUNDARY,
  P4EST_TRAVERSE_INNER
}
p4est_meshiter_parallelboundary_t;

/** An iterator based on the information stored in p4est_mesh
 * The general idea how to use it is to create an iterator for a specific level
 * and a specific behaviour regarding which cells to visit.
 * Then traverse the grid using something like
 * do {
 *   ret_val = p4est_mesh_iterator_next(iterator);
 *   function_for_current_cell (iterator);
 * } (while ret_val != P4EST_MESHITER_DONE);
 * In the function for the current cell neighbor entities may be accessed up to
 * the specified connect type \a btype.
 *
 * p4est, ghost and mesh hold pointers to the respective p4est structures from
 * which the iterator is built. btype must match with the values used to create
 * ghost and mesh.
 *
 * The values traverse_ghost, traverse_virtuals, and traverse_parallel_boundary
 * must be within 0 .. 2.
 * In case of traverse_ghost a value of 0 indicates to only traverse local
 * quadrants, a value of 2 to only traverse ghost quadrants, and a value of 1 to
 * traverse both, local and ghost quadrants.
 * In case of traverse_virtuals a value of 0 indicates to only traverse real
 * quadrants, a value of 2 to only traverse virtual quadrants, and a value of 1
 * to traverse both, real and virtual quadrants.
 * In case of traverse_parallel_boundary a value of 0 indicates to traverse all
 * quadrants, a value of 1 to only traverse mirror quadrants, i.e. local
 * quadrants along the parallel process boundary, and a value of 2 to only
 * traverse inner quadrants, i.e. all local quadrants that are not directly
 * adjacent to the parallel process boundary.
 * Note that these flags may be combined in an arbitrary way.
 *
 * neighbor_qid, neighbor_vid, and neighbor_enc are sc_arrays used to store the
 * results of neighborsearch queries.
 *
 * The internal status of the mesh-based iterator consists of 3 parts:
 *
 * 1.) Global information
 * Global information is created once and must not be touched from outside
 * during the iteration over the selected set of quadrants.
 * The number of real local quadrants is stored in real_on_level_local, the
 * number of ghosts is stored in real_on_level_ghost. If real quadrants or
 * either ghost or local quadrants are not traversed the respective values are
 * set to 0.
 * The number of virtual local quadrants is stored in virtual_on_level_local,
 * the number of ghosts is stored in virtual_on_level_ghost. If virtual
 * quadrants or either ghost or local quadrants are not traversed the respective
 * values are set to 0.
 * The level of quadrants which the iterator visits is stored in current_level.
 * To keep track of the iterators status 2 counters are stored:
 * The number of traversed real cells is stored in processed_real, the number of
 * traversed cells containing virtual subquadrants is stored in
 * processed_virtual.
 * CAUTION: processed_virtual only increments after a complete block of virtual
 * subquadrants has been processed, i.e. the storage index is to be multiplied
 * by 4.
 * The maximum number of iterations is stored in maximum_index. As we have to
 * call p4est_meshiter_next() to change subquadrants, the number of subquadrants
 * is counted here instead of the p4est_quadrant_t that contain the respective 4
 * virtual subquadrants.
 *
 * 2.) Information about the currently processed quadrant
 * The number of quadrants the iterator already has visited is stored in
 * current_index. This is the value that is compared to maximum_index to decide
 * if the iterator is done.
 * current_qid holds the quadrant id of the currently processed quadrant. If we
 * are currently processing a set of virtual quadrants the quadrants virtual id
 * is stored in current_vid. If a real quadrant is processed, this value is set
 * to -1.
 * If the our neighbor is bigger, the current quadrant is part of a set of
 * hanging quadrants which are indexed in Morton-order. To find the appropriate
 * virtual subquadrant of the neighbor to interact with, current_subquad stores
 * the index of the current quadrant in the list of hanging quadrants.
 * If local and ghost quadrants should be traversed then ghost quadrants are
 * traversed after local quadrants. current_is_ghost stores 0 if there are local
 * quadrants left to process and 1 if all local quadrants are already traverses.
 *
 * 3.) Information about a neighboring quadrant in a specific direction.
 * The information about the neighboring quadrant of size current_level is
 * obtained from the mesh and stored in the same way as for the current
 * quadrant. Additionally we store the direction in which the neighbor is
 * situated w.r.t. the current quadrant:
 *   0 .. 3 neighbor is in direction of face f_i
 *   4 .. 7 neighbor is in direction of corner c_(i-4)
 * As the neighbor may be rotated w.r.t. the current quadrant the neighboring
 * quadrant's orientation is stored in neighbor_orientation and the entity in
 * which the current quadrant is positioned w.r.t. the neighboring quadrant is
 * stored in neighbor_entity_index, whose value must be interpreteted the same
 * way as the value of neighbor_direction.
 */
typedef struct
{
  /** forest information */
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  p4est_connect_type_t btype;

  /** internal neighbor access operators */
  sc_array_t          neighbor_qids;       /**< array of type int */
  sc_array_t          neighbor_vids;       /**< array of type int */
  sc_array_t          neighbor_encs;        /**< array of type int */

  /** flags to include or exclude certain cells */
  int8_t              traverse_ghosts;     /**< traverse ghosts:
                                                0: only local quads
                                                1: local and ghost quads,
                                                2: only ghost quads */
  int8_t              traverse_virtuals;   /**< traverse virtual quadrants:
                                                0: only real quads
                                                1: real and virtual quads
                                                2: only virtual quads */
  int8_t              traverse_parallel_boundary; /**< traverse parallel
                                                       boundary:
                                                       0: all quads,
                                                       1: quads on boundary,
                                                       2: inner quads */

  /** internal state */
  /** global */
  int                 real_on_level_local;
  int                 real_on_level_ghost;
  int                 virtual_on_level_local;
  int                 virtual_on_level_ghost;
  int                 current_level;
  p4est_locidx_t      processed_real;
  p4est_locidx_t      processed_virtual;
  int                 maximum_index;

  /** current */
  p4est_locidx_t      current_index;
  p4est_locidx_t      current_qid;
  int                 current_vid;
  int                 current_subquad;
  int                 current_is_ghost;

  /** current neighbor, only set when current not ghost */
  p4est_locidx_t      neighbor_qid;
  int                 neighbor_is_ghost;
  int                 neighbor_direction;
  int                 neighbor_entity_index;
  int                 neighbor_orientation;
  int                 neighbor_vid;
  int                 neighbor_subquad;
}
p4est_meshiter_t;

/** Create a p4est_mesh_iter structure for a given level.
 * \param [in] p4est   A forest that is 2:1 balanced up to the specified
 *                     codimension of neighbors.
 * \param [in] ghost   The ghost layer created from the provided p4est.
 * \param [in] mesh    The mesh created from the provided p4est and ghost.
 *                     CAUTION: mesh must contain quad_level and ghost_level
 *                              fields.
 * \param [in] virtual_quads  Virtual quadrants at refinement boundaries.
 *                            CAUTION: virtual_quadrants must contain
 *                                     virtual_qlevels and virtual_glevels
 *                                     fields.
 * \param [in] level   The level to be traversed
 * \param [in] btype   Defines the highest codimension of neighbors
 */
p4est_meshiter_t   *p4est_meshiter_new (p4est_t * p4est,
                                        p4est_ghost_t * ghost,
                                        p4est_mesh_t * mesh,
                                        p4est_virtual_t * virtual_quads,
                                        int level,
                                        p4est_connect_type_t btype);

/** Free the allocated mesh_iter struct
 * \param [in] mesh_iter     The meshiter struct to be freed
 */
void                p4est_meshiter_destroy (p4est_meshiter_t * mesh_iter);

/** Move iterator one quad further
 * \param[in][out] mesh_iter  The mesh iterator to be moved
 * \returns   0 if all went well or P4EST_MESHITER_DONE if the last quad was
 *            already processed
 */
int                 p4est_meshiter_next (p4est_meshiter_t * mesh_iter);

/** Populate the information about the neighboring quad of the current quad in
 * the specified direction.
 * \param [in][out] mesh_iter  The mesh iterator to set the information in.
 * \param [in]      direction  The direction in which a neighboring quad shall
 *                             be found:
 *                             0 .. 3 direction of face f_i
 *                             4 .. 8 direction of corner c_(i-4)
 */
int                 p4est_meshiter_set_neighbor_quad_info (p4est_meshiter_t *
                                                           mesh_iter,
                                                           int direction);

/** Get storage id of current quadrant's payload in a per-level data structure.
 * \param[in]       mesh_iter  The mesh iterator
 */
int                 p4est_meshiter_get_current_storage_id (p4est_meshiter_t *
                                                           mesh_iter);

/** Get storage id of neighboring quadrant's payload in a per-level data
 * structure.
 * \param[in]       mesh_iter  The mesh iterator
 */
int                 p4est_meshiter_get_neighbor_storage_id (p4est_meshiter_t *
                                                            mesh_iter);

SC_EXTERN_C_END;

#endif /* P4EST_MESHITERATE_H */
