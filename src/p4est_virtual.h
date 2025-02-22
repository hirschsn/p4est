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

#ifdef P4_TO_P8
#error "Including a p4est header with P4_TO_P8 defined"
#endif /* P4_TO_P8 */

#include <p4est_ghost.h>
#include <p4est_mesh.h>

SC_EXTERN_C_BEGIN;

/* -------------------------------------------------------------------------- */
/* |                           Virtual quadrants                            | */
/* -------------------------------------------------------------------------- */
/** Struct for storing virtual quadrants at refinement boundaries in coarse
 * quadrants for algorithms that rely on interaction to only take place between
 * same-sized quadrants.
 *
 * btype stores which quadrants we consider to decide if quadrant contains
 * virtual quadrants.
 *
 * If a quadrant contains virtual quadrants is stored in virtual_qflags for
 * local quadrants and virtual_gflags for ghost quadrants.  If the quadrant
 * contains virtual quadrants the number of coarse quadrants at a refinement
 * boundary in Morton order, i.e. the number of quadrants hosting virtual
 * quadrants up to the current quadrant. If the quadrant does not contain
 * virtual quadrants, we store -1.
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
 * These arrays are NULL by default and can be enabled using \ref
 * p4est_virtual_new_ext.
 */
typedef struct p4est_virtual
{
  p4est_locidx_t      local_num_quadrants;
  p4est_locidx_t      ghost_num_quadrants;
  p4est_locidx_t      local_num_virtuals;
  p4est_locidx_t      ghost_num_virtuals;

  p4est_connect_type_t btype;           /**< which neighbors are considered in
                                             virtual */

  p4est_locidx_t     *virtual_qflags;   /**< stores if quad has virtual quads
                                             and how many local quadrants
                                             containing virtual quadrants are
                                             preceding it in Morton-order */
  p4est_locidx_t     *virtual_gflags;   /**< stores if ghost has virtual quads
                                             and how many ghosts containing
                                             virtual quadrants are preceding
                                             it in Morton-order */
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
                                                     or virtual.  NULL by
                                                     default, can be enabled in
                                                     \ref p4est_virtual_new_ext */
  p4est_locidx_t     *quad_qvirtual_offset;     /**< stores number of quadrants
                                                     preceding first virtual
                                                     sub-quadrant of current
                                                     quadrant on the level of
                                                     its virtual sub-quads if it
                                                     has any, either real or
                                                     virtual.  Else -1.  NULL
                                                     by default, can be enabled
                                                     in \ref
                                                     p4est_virtual_new_ext */
  p4est_locidx_t     *quad_greal_offset;        /**< stores number of quadrants
                                                     preceding current ghost
                                                     on its level, either real
                                                     or virtual.  NULL by
                                                     default, can be enabled in
                                                     \ref p4est_virtual_new_ext */
  p4est_locidx_t     *quad_gvirtual_offset;     /**< stores number of quadrants
                                                     preceding first virtual
                                                     sub-quadrant of current
                                                     ghost on the level of its
                                                     virtual sub-quads if it has
                                                     any, either real or
                                                     virtual.  Else -1.  NULL
                                                     by default, can be enabled
                                                     in \ref
                                                     p4est_virtual_new_ext */
} p4est_virtual_t;

/** Create a p4est_virtual structure.
 * This function does not populate virtual_qlevels and virtual_glevels fields.
 * To populate them, use \ref p4est_virtual_new_ext.
 * \param[in] p4est    A forest that is 2:1 balanced up to \b btype.
 * \param[in] ghost    The ghost layer created from the provided p4est.
 * \param[in] mesh     The neighboring lookup tables created from provided p4est
 *                     and ghost.
 * \param[in] btype    The highest codimension of neighbors to consider for
 *                     embedding virtual quadrants.
 * \return             A fully allocated structure of virtual quadrants.
 */
p4est_virtual_t    *p4est_virtual_new (p4est_t * p4est, p4est_ghost_t * ghost,
                                       p4est_mesh_t * mesh,
                                       p4est_connect_type_t btype);

/** Destroy a p4est_virtual structure.
 * \param [in] virtual_quads   Virtual structure previously created by
 *                             p4est_virtual_new or p4est_virtual_new_ext.
 */
void                p4est_virtual_destroy (p4est_virtual_t * virtual_quads);

/** Calculate the memory usage of the virtual structure.
 * \param[in] virtual_quads   Virtual structure
 * \return                    Memory used in bytes.
 */
size_t              p4est_virtual_memory_used (p4est_virtual_t *
                                               virtual_quads);

/* -------------------------------------------------------------------------- */
/* |                             Ghost exchange                             | */
/* -------------------------------------------------------------------------- */

/** Store all required information for exchanging ghost data including that of
 * virtual quadrants.  The maximum codimension which we consider is stored in
 * btype.
 * mirror_proc_virtuals is structured and accessed in the same way as
 * mirror_proc_mirrors in p8est_ghost_t.
 * If either no virtual quadrants are embedded or no virtual quadrants are
 * expected 0 is stored, else 1.
 */
typedef struct p4est_virtual_ghost
{
  p4est_connect_type_t btype;
  int8_t             *mirror_proc_virtuals;
} p4est_virtual_ghost_t;

/** Transient storage for asynchronous ghost exchange. */
typedef struct p4est_virtual_ghost_exchange
{
  int                 is_levels;
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_virtual_t    *virtual_quads;
  p4est_virtual_ghost_t *virtual_ghost;
  int                 level;
  size_t              data_size;
  void              **ghost_data;
  sc_array_t          requests, sbuffers;
} p4est_virtual_ghost_exchange_t;

/** Extend p8est_ghost such that payload can be exchanged including virtual
 * quadrants.
 * \param[in]      p4est            The forest, must be 2:1 balanced.
 * \param[in]      ghost            The forest's ghost layer.
 * \param[in]      mesh             Lookup tables encoding neighbor information.
 *                 CAUTION: mirror_qid array needs to be populated.
 * \param[in]      virtual_quads    Which quadrants host virtual quadrants, both
 *                                  local and ghost quadrants.
 * \param[in]      btype            The maximum co-dimension of neighbors that
 *                                  will be considered.
 * \return                          Struct containing all necessary information
 *                                  for exchanging data including virtual quads.
 */
p4est_virtual_ghost_t *p4est_virtual_ghost_new (p4est_t * p4est,
                                                p4est_ghost_t * ghost,
                                                p4est_mesh_t * mesh,
                                                p4est_virtual_t *
                                                virtual_quads,
                                                p4est_connect_type_t btype);

/** Free a virtual_ghost struct previously allocated using
 * \ref p4est_virtual_ghost_new.
 * \param[in]     virtual_ghost    The struct to be freed.
 */
void                p4est_virtual_ghost_destroy (p4est_virtual_ghost_t *
                                                 virtual_ghost);

/** Return memory size in bytes that is occupied from p8est_ghost_virtual_t
 * struct.
 * CAUTION: Does not work yet.
 * \param[in] virtual_ghost   Virtual ghost structure.
 * \return                    Memory used in bytes.
 */
size_t              p4est_virtual_ghost_memory_used (p4est_virtual_ghost_t *
                                                     virtual_ghost);

/* *INDENT-OFF* */
/** Convenience function for exchanging contiguous data including that of
 * virtual quadrants.  Using this function prevents overlapping computation and
 * communication.
 * \param[in]      p4est           The forest.
 * \param[in]      ghost           Ghost layer.
 * \param[in]      mesh            Neighbor information.
 * \param[in]      virtual_quads   Virtual quadrants.
 * \param[in]      virtual_ghost   Extension of ghost layer by virtual quadrants
 * \param[in]      data_size       Size of one local data element
 * \param[in]      mirror_data     Contiguous data of local quadrants.
 * \param[in][out] ghost_data      Contiguous data storage that will be
 *                                 populated during ghost exchange.
 * \return                         Transient storage
 */
void                p4est_virtual_ghost_exchange_data
  (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
   p4est_virtual_t * virtual_quads, p4est_virtual_ghost_t * virtual_ghost,
   size_t data_size, void *mirror_data, void *ghost_data);

/** Initiate message transfer by posting asynchronous sends and receive
 * operations for contiguous data.
 * \param[in]      p4est           The forest.
 * \param[in]      ghost           Ghost layer.
 * \param[in]      mesh            Neighbor information.
 * \param[in]      virtual_quads   Virtual quadrants.
 * \param[in]      virtual_ghost   Extension of ghost layer by virtual quadrants
 * \param[in]      data_size       Size of one local data element
 * \param[in]      mirror_data     Contiguous data of local quadrants.
 * \param[in][out] ghost_data      Contiguous data storage that will be
 *                                 populated during ghost exchange.
 * \return                         Transient storage
 */
p4est_virtual_ghost_exchange_t
  * p4est_virtual_ghost_exchange_data_begin
  (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
   p4est_virtual_t * virtual_quads, p4est_virtual_ghost_t * virtual_ghost,
   size_t data_size, void *mirror_data, void *ghost_data);

/** Finish message transfer and clean up transient storage for contiguous data.
 * \param[in] exc      Transient storage of data transfer.  Will be freed during
 *                     execution.
 */
void               p4est_virtual_ghost_exchange_data_end
  (p4est_virtual_ghost_exchange_t * exc);

/** Convenience function for exchanging per-level data including that of
 * virtual quadrants.  Using this function prevents overlapping computation and
 * communication.
 * \param[in]      p4est           The forest.
 * \param[in]      ghost           Ghost layer.
 * \param[in]      mesh            Neighbor information.
 * \param[in]      virtual_quads   Virtual quadrants.
 * \param[in]      virtual_ghost   Extension of ghost layer by virtual quadrants
 * \param[in]      level           The level of quadrants whose data will be
 *                                 exchanged.
 * \param[in]      data_size       Size of one local data element
 * \param[in]      mirror_data     Per-level data of local quadrants.
 * \param[in][out] ghost_data      Per-level data storage that will be populated
 *                                 during ghost exchange.
 * \return                         Transient storage
 */
void                p4est_virtual_ghost_exchange_data_level
  (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
   p4est_virtual_t * virtual_quads, p4est_virtual_ghost_t * virtual_ghost,
   int level, size_t data_size, void **mirror_data, void **ghost_data);

/** Initiate message transfer by posting asynchronous sends and receive
 * operations for per-level data.
 * \param[in]      p4est           The forest.
 * \param[in]      ghost           Ghost layer.
 * \param[in]      mesh            Neighbor information.
 * \param[in]      virtual_quads   Virtual quadrants.
 * \param[in]      virtual_ghost   Extension of ghost layer by virtual quadrants
 * \param[in]      level           The level of quadrants whose data will be
 *                                 exchanged.
 * \param[in]      data_size       Size of one local data element
 * \param[in]      mirror_data     Per-level data of local quadrants.
 * \param[in][out] ghost_data      Per-level data storage that will be populated
 *                                 during ghost exchange.
 * \return                         Transient storage
 * \return                    Transient storage
 */
p4est_virtual_ghost_exchange_t
  * p4est_virtual_ghost_exchange_data_level_begin
  (p4est_t * p4est, p4est_ghost_t * ghost, p4est_mesh_t * mesh,
   p4est_virtual_t * virtual_quads, p4est_virtual_ghost_t * virtual_ghost,
   int level, size_t data_size, void **mirror_data, void **ghost_data);

/** Finish message transfer and clean up transient storage for per-level data.
 * \param[in] exc      Transient storage of data transfer.  Will be freed during
 *                     execution.
 */
void                p4est_virtual_ghost_exchange_data_level_end
  (p4est_virtual_ghost_exchange_t * exc);
/* *INDENT-ON* */

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
extern const int    p4est_virtual_face_neighbors_search_opts[P4EST_CHILDREN]
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
extern const int    p4est_virtual_corner_neighbors_search_opts[P4EST_CHILDREN]
                                                              [P4EST_CHILDREN];
/* *INDENT-ON* */

/** Neighbor lookup including virtual quadrants.  This function only returns
 * neighbors of the same size for virtual and real quadrants.
 * \param[in]      p4est          The forest, must be 2:1 balanced.
 * \param[in]      ghost          Ghost layer.
 * \param[in]      mesh           Neighboring information of real quadrants.
 * \param[in]      virtual_quads  Embedded virtual quadrants at refinement
 *                                boundaries.
 * \param[in]      qid            Local quadrant id of current quadrant.
 * \param[in]      vid            Virtual id of current quadrant.  -1 for real
 *                                quadrants.
 * \param[in]      dir            Direction in which to search for neighbors:
 *                                  0 .. 3 neighbor(-s) across face i,
 *                                  4 .. 7 neighbor(-s) across corner i-4.
 * \param    [out] n_encs         Array containing encodings for neighboring
 *                                quadrants as it is described in
 *                                \ref p4est_mesh_t.
 *                                Array must be empty and allocated for ints.
 * \param    [out] n_qids         Array containing neighboring quadrant ids.
 *                                Array must be empty and allocated for
 *                                p4est_locidx_t.
 * \param    [out] n_vids         Array containing neighboring quadrant's
 *                                virtual ids.  Real quadrants obtain vid -1.
 *                                Array must be empty and allocated for ints.
 */
int                 p4est_virtual_get_neighbor (p4est_t * p4est,
                                                p4est_ghost_t * ghost,
                                                p4est_mesh_t * mesh,
                                                p4est_virtual_t *
                                                virtual_quads,
                                                p4est_locidx_t qid, int vid,
                                                int dir,
                                                sc_array_t * n_encs,
                                                sc_array_t * n_qids,
                                                sc_array_t * n_vids);
SC_EXTERN_C_END;

#endif /* P4EST_VIRTUAL_H */
