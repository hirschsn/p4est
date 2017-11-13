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

#ifndef P4_TO_P8
#include <p4est_meshiter.h>
#include <p4est_extended.h>
#else /* P4_TO_P8 */
#include <p8est_meshiter.h>
#include <p8est_extended.h>
#endif /* P4_TO_P8 */

/** Function to decode entries of p4/8est_face_virtual_neighbors_inside,
 * p8est_edge_virtual_neigbors_inside, and
 * p4/8est_corner_virtual_neighbors_inside.
 * \param [in]        internal_neighbor_id  Entry to be decoded.
 * \param     [out]   direction             Give encoded direction for calling
 *                                          p4est_mesh_get_neighbors.
 * \param     [out]   subentity             Subindex to transform to adjacent
 *                                          tree for obtaining correct hanging
 *                                          quadrant.
 */
static int
decode_internal_neighbor_id (int internal_neighbor_id, int *direction,
                             int *subentity)
{
  P4EST_ASSERT (P4EST_CHILDREN <= internal_neighbor_id);
#ifdef P4_TO_P8
  /** faces */
  const int           face2edge = 32;
  const int           edge2corner = 56;
  const int           corner2invalid = 64;
  if (P4EST_CHILDREN <= internal_neighbor_id
      && internal_neighbor_id < face2edge) {
    *direction = (internal_neighbor_id - P4EST_CHILDREN) % P4EST_FACES;
    *subentity = (internal_neighbor_id - P4EST_CHILDREN) / P4EST_FACES;
  }
  else if (face2edge <= internal_neighbor_id
           && internal_neighbor_id < edge2corner) {
    *direction =
      P4EST_FACES + ((internal_neighbor_id - face2edge) % P8EST_EDGES);
    *subentity = (internal_neighbor_id - face2edge) / P8EST_EDGES;
  }
  else if (edge2corner <= internal_neighbor_id
           && internal_neighbor_id < corner2invalid) {
    *direction =
      P4EST_FACES + P8EST_EDGES + internal_neighbor_id - edge2corner;
    *subentity = -1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }

#else /* P4_TO_P8 */
  const int           face2corner = 12;
  const int           corner2invalid = 16;
  if (P4EST_CHILDREN <= internal_neighbor_id
      && internal_neighbor_id < face2corner) {
    *direction = (internal_neighbor_id - P4EST_CHILDREN) % P4EST_FACES;
    *subentity = (internal_neighbor_id - P4EST_CHILDREN) / P4EST_FACES;
  }
  else if (face2corner <= internal_neighbor_id
           && internal_neighbor_id < corner2invalid) {
    *direction = P4EST_FACES + internal_neighbor_id - face2corner;
    *subentity = -1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#endif /* P4_TO_P8 */

  return 0;
}

static int
set_limits (int direction, p4est_connect_type_t btype, int *l_same_size,
            int *u_same_size, int *l_double_size, int *u_double_size,
            int *l_half_size, int *u_half_size, int *n_neighbor_entities,
            int *n_hanging_quads)
{
#ifndef P4_TO_P8
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = *l_double_size = 8;
    *u_double_size = *l_half_size = 24;
    *u_half_size = 32;
    *n_neighbor_entities = P4EST_FACES;
    *n_hanging_quads = P4EST_HALF;
  }
  else if (P4EST_FACES <= direction
           && direction < (P4EST_FACES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_same_size = 0;
    *u_same_size = 4;
    *l_double_size = *u_double_size = *l_half_size = *u_half_size = 4;
    *n_neighbor_entities = P4EST_CHILDREN;
    *n_hanging_quads = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#else /* !P4_TO_P8 */
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = *l_double_size = 24;
    *u_double_size = *l_half_size = 120;
    *u_half_size = 144;
    *n_neighbor_entities = P4EST_FACES;
    *n_hanging_quads = P4EST_HALF;
  }
  else if (P4EST_FACES <= direction
           && direction < (P4EST_FACES + P8EST_EDGES)) {
    P4EST_ASSERT (btype >= P8EST_CONNECT_EDGE);
    *l_same_size = 0;
    *u_same_size = *l_double_size = 24;
    *u_double_size = *l_half_size = 72;
    *u_half_size = 96;
    *n_neighbor_entities = P8EST_EDGES;
    *n_hanging_quads = 2;
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= direction
           && direction < (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_same_size = 0;
    *u_same_size = 8;
    *l_double_size = *u_double_size = *l_half_size = *u_half_size = 8;
    *n_neighbor_entities = P4EST_CHILDREN;
    *n_hanging_quads = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#endif /* !P4_TO_P8 */

  return 0;
}

/** Get opposite face-, (edge-), or corner index
 * \param [in]   dir  Encoded direction.
 *
 * \returns Decoded index of opposite face (, edge,) or corner.
 */
static int
get_opposite (int dir)
{
  if (0 <= dir && dir < P4EST_FACES) {
    return dir ^ 1;
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= dir && dir < (P4EST_FACES + P8EST_EDGES)) {
    return (dir - P4EST_FACES) ^ 3;
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= dir
           && dir < (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    return (dir - (P4EST_FACES + P8EST_EDGES)) ^ 7;
  }
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= dir && dir < (P4EST_FACES + P4EST_CHILDREN)) {
    return (dir - P4EST_FACES) ^ 3;
  }
#endif /* P4_TO_P8 */
  else {
    SC_ABORT_NOT_REACHED ();
  }
  return -42;
}

/** Get neighbor virtual id in case current quadrant is hanging wrt. neighbor.
 * \param[in] mesh_iter   Instance of the current mesh based iterator with set
 *                        internal status for current and neighboring quadrant.
 */
static int
get_vid (p4est_meshiter_t * mesh_iter)
{
  P4EST_ASSERT (mesh_iter->neighbor_subquad != -2);
  int                 ret;

  if (0 <= mesh_iter->neighbor_direction
      && mesh_iter->neighbor_direction < P4EST_FACES) {
    /** faces */
    if (0 <= mesh_iter->neighbor_subquad
        && mesh_iter->neighbor_subquad < P4EST_HALF) {
      /** find neighboring vid from own subquad ID */
      ret =
        p4est_connectivity_face_neighbor_face_corner_orientation
        (mesh_iter->neighbor_subquad, mesh_iter->neighbor_direction,
         mesh_iter->neighbor_entity_index, mesh_iter->neighbor_orientation);
      ret = p4est_face_corners[mesh_iter->neighbor_entity_index][ret];
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction < (P4EST_FACES + P8EST_EDGES)) {
    /** edges */
    if (0 <= mesh_iter->neighbor_subquad && mesh_iter->neighbor_subquad < 2) {
      /** find neighboring vid from own subquad ID */
      ret =
        p8est_connectivity_edge_neighbor_edge_corner_orientation
        (mesh_iter->neighbor_subquad,
         mesh_iter->neighbor_direction - P4EST_FACES,
         mesh_iter->neighbor_entity_index, mesh_iter->neighbor_orientation);
      ret =
        p8est_edge_corners[mesh_iter->neighbor_entity_index -
                           P4EST_FACES][ret];
    }
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    /** corners */
    ret = mesh_iter->neighbor_entity_index;
  }
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P4EST_CHILDREN)) {
    ret = mesh_iter->neighbor_entity_index;
  }
#endif /* P4_TO_P8 */
  else {
    SC_ABORT_NOT_REACHED ();
  }

  return ret;
}

/* CAUTION: Only valid if called from a virtual quadrant */
static int
get_subentity (p4est_meshiter_t * mesh_iter)
{
  P4EST_ASSERT (mesh_iter->current_subquad == -2);
  P4EST_ASSERT (0 <= mesh_iter->current_vid
                && mesh_iter->current_vid < P4EST_CHILDREN);

  int                 ret;
  if (0 <= mesh_iter->neighbor_direction
      && mesh_iter->neighbor_direction < P4EST_FACES) {
    /** find neighboring vid from own vid across face */
    ret =
      p4est_connectivity_face_neighbor_corner_orientation
      (mesh_iter->current_vid, mesh_iter->neighbor_direction,
       mesh_iter->neighbor_entity_index, mesh_iter->neighbor_orientation);
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction < (P4EST_FACES + P8EST_EDGES)) {
    /** find neighboring vid from own vid across edge */
    ret =
      p8est_connectivity_edge_neighbor_corner_orientation
      (mesh_iter->current_vid, mesh_iter->neighbor_direction - P4EST_FACES,
       mesh_iter->neighbor_entity_index, mesh_iter->neighbor_orientation);
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    /** return corner directly */
    ret = mesh_iter->neighbor_entity_index;
  }
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P4EST_CHILDREN)) {
    /** return corner directly */
    ret = mesh_iter->neighbor_entity_index;
  }
#endif /* P4_TO_P8 */
  else {
    SC_ABORT_NOT_REACHED ();
  }

  return ret;
}

/** Decode encoding obtained in neighbor search
 * \param[in][out] mesh_iter     Instance of mesh-based iterator. Neighbor
 *                               information will be set according to decoding.
 * \param[in]      enc           The normalized encoding, i.e. 0 based and no
 *                               longer containing ghost status, i.e. 0 <= enc
 * \param[in]      n_entities    Number of faces, edges, or corners, depending
 *                               on the direction the neighbor has been looked
 *                               up.
 * \param[in]      l_same_size   Lower bound for encoding a neighbor of same
 *                               size.
 * \param[in]      u_same_size   Upper bound for encoding a neighbor of same
 *                               size.
 * \param[in]      l_double_size Lower bound for encoding a neighbor of double
 *                               size.
 * \param[in]      u_double_size Upper bound for encoding a neighbor of double
 *                               size.
 * \param[in]      l_half_size   Lower bound for encoding a neighbor of half
 *                               size.
 * \param[in]      u_half_size   Upper bound for encoding a neighbor of half
 *                               size.
 */
static int
p4est_meshiter_decode_encoding (p4est_meshiter_t * mesh_iter, int enc,
                                int n_entities, int l_same_size,
                                int u_same_size, int l_double_size,
                                int u_double_size, int l_half_size,
                                int u_half_size)
{
  P4EST_ASSERT (u_same_size == l_double_size && u_double_size == l_half_size);
  P4EST_ASSERT (l_same_size == 0);
  P4EST_ASSERT (l_same_size < l_double_size);
  P4EST_ASSERT (l_double_size <= l_half_size);
  P4EST_ASSERT (0 <= enc);

  if (l_same_size <= enc && enc < u_same_size) {
    mesh_iter->neighbor_orientation = enc / n_entities;
    mesh_iter->neighbor_entity_index = enc % n_entities;
  }
  else if (l_double_size <= enc && enc < u_double_size) {
    int                 e = enc;
    e -= l_double_size;
    mesh_iter->neighbor_subquad = e / l_double_size;
    e -= (l_double_size * mesh_iter->neighbor_subquad);
    mesh_iter->neighbor_orientation = e / n_entities;
    mesh_iter->neighbor_entity_index = e % n_entities;
  }
  else if (l_half_size <= enc && enc < u_half_size) {
    int                 e = enc;
    e -= l_half_size;
    mesh_iter->neighbor_orientation = e / n_entities;
    mesh_iter->neighbor_entity_index = e % n_entities;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }

#ifdef P4EST_ENABLE_DEBUG
  int8_t              upper_bnd;
  if (n_entities == P4EST_FACES) {
    upper_bnd = P4EST_HALF;
  }
#ifdef P4_TO_P8
  else if (n_entities == P8EST_EDGES) {
    upper_bnd = 2;
  }
#endif /* P4_TO_P8 */
  else if (n_entities == P4EST_CHILDREN) {
    upper_bnd = 1;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }

  P4EST_ASSERT (0 <= mesh_iter->neighbor_orientation &&
                mesh_iter->neighbor_orientation < upper_bnd);
#endif /* P4EST_ENABLE_DEBUG */
  P4EST_ASSERT (0 <= mesh_iter->neighbor_entity_index &&
                mesh_iter->neighbor_entity_index < n_entities);

  int                 entity_offset;
  if (0 <= mesh_iter->neighbor_direction
      && mesh_iter->neighbor_direction < P4EST_FACES) {
    entity_offset = 0;
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction < (P4EST_FACES + P8EST_EDGES)) {
    entity_offset = P4EST_FACES;
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    entity_offset = P4EST_FACES + P8EST_EDGES;
  }
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= mesh_iter->neighbor_direction
           && mesh_iter->neighbor_direction <
           (P4EST_FACES + P4EST_CHILDREN)) {
    entity_offset = P4EST_FACES;
  }
#endif /* P4_TO_P8 */
  else {
    SC_ABORT_NOT_REACHED ();
  }
  mesh_iter->neighbor_entity_index += entity_offset;

  return 0;
}

/** Discard all neighbor information
 *
 * \param [in][out]  mesh_iter   mesh_iter whose neighbor info is unset
 */
static int
discard_neighbors (p4est_meshiter_t * mesh_iter)
{
  sc_array_truncate (&mesh_iter->neighbor_ids);
  sc_array_truncate (&mesh_iter->neighbor_enc);

  mesh_iter->neighbor_orientation = -1;
  mesh_iter->neighbor_qid = -1;
  mesh_iter->neighbor_subquad = -2;
  mesh_iter->neighbor_vid = -1;
  mesh_iter->neighbor_is_ghost = -1;
  mesh_iter->neighbor_direction = -1;
  mesh_iter->neighbor_entity_index = -1;

  return 0;
}

p4est_meshiter_t   *
p4est_meshiter_new (p4est_t * p4est, p4est_ghost_t * ghost,
                    p4est_mesh_t * mesh, int level,
                    p4est_connect_type_t btype)
{
  return p4est_meshiter_new_ext (p4est, ghost, mesh, level, btype,
                                 P4EST_TRAVERSE_LOCAL, P4EST_TRAVERSE_REAL,
                                 P4EST_TRAVERSE_PARBOUNDINNER);
}

p4est_meshiter_t   *
p4est_meshiter_new_ext (p4est_t * p4est, p4est_ghost_t * ghost,
                        p4est_mesh_t * mesh, int level,
                        p4est_connect_type_t btype,
                        p4est_meshiter_localghost_t traverse_ghosts,
                        p4est_meshiter_realvirt_t traverse_virtuals,
                        p4est_meshiter_parallelboundary_t
                        traverse_parallel_boundary)
{
  p4est_meshiter_t   *mesh_iter;
  mesh_iter = P4EST_ALLOC_ZERO (p4est_meshiter_t, 1);

  /* pass parameters */
  mesh_iter->p4est = p4est;
  mesh_iter->ghost = ghost;
  mesh_iter->mesh = mesh;
  mesh_iter->btype = btype;
  mesh_iter->current_level = level;
  mesh_iter->traverse_ghosts = traverse_ghosts;
  mesh_iter->traverse_virtuals = traverse_virtuals;
  mesh_iter->traverse_parallel_boundary = traverse_parallel_boundary;

  /* sanity checks that passed parameters are OK */
  P4EST_ASSERT (btype <= mesh_iter->ghost->btype);
  P4EST_ASSERT (mesh_iter->mesh->quad_level != NULL);

#ifdef P4_TO_P8
  if (P8EST_CONNECT_EDGE <= btype) {
    P4EST_ASSERT (mesh_iter->mesh->quad_to_edge != NULL);
  }
#endif /* P4_TO_P8 */
  if (P4EST_CONNECT_CORNER == btype) {
    P4EST_ASSERT (mesh_iter->mesh->quad_to_corner != NULL);
  }

  /* set internal state */
  mesh_iter->real_on_level_local =
    ((traverse_ghosts == 2) || (traverse_virtuals == 2) ? 0 :
     (mesh_iter->mesh->quad_level + mesh_iter->current_level)->elem_count);

  mesh_iter->real_on_level_ghost =
    ((traverse_ghosts == 0) || (traverse_virtuals == 2) ? 0 :
     (mesh_iter->mesh->ghost_level + mesh_iter->current_level)->elem_count);

  mesh_iter->virtual_on_level_local =
    ((traverse_ghosts == 2) || (traverse_virtuals == 0) ? 0 :
     (mesh_iter->mesh->virtual_qlevels +
      mesh_iter->current_level)->elem_count);

  mesh_iter->virtual_on_level_ghost =
    ((traverse_ghosts == 0) || (traverse_virtuals == 0) ? 0 :
     (mesh_iter->mesh->virtual_glevels +
      mesh_iter->current_level)->elem_count);

  mesh_iter->maximum_index =
    mesh_iter->real_on_level_local + mesh_iter->real_on_level_ghost +
    P4EST_CHILDREN * mesh_iter->virtual_on_level_local +
    P4EST_CHILDREN * mesh_iter->virtual_on_level_ghost;

  sc_array_init (&mesh_iter->neighbor_ids, sizeof (int));
  sc_array_init (&mesh_iter->neighbor_enc, sizeof (int));

  mesh_iter->processed_real = 0;
  mesh_iter->processed_virtual = 0;

  mesh_iter->current_index = -1;
  mesh_iter->current_vid = -1;
  mesh_iter->current_qid = -1;
  mesh_iter->current_subquad = -1;
  mesh_iter->current_is_ghost = 0;

  discard_neighbors (mesh_iter);

  return mesh_iter;
}

void
p4est_meshiter_destroy (p4est_meshiter_t * mesh_iter)
{
  sc_array_reset (&mesh_iter->neighbor_ids);
  sc_array_reset (&mesh_iter->neighbor_enc);
  P4EST_FREE (mesh_iter);
}

int
p4est_meshiter_next (p4est_meshiter_t * mesh_iter)
{
  int                 found_valid_quad = 0;
  int                 virt_qid, real_qid;
  int                 max_virt, max_real;

  /** unset neighbor information */
  discard_neighbors (mesh_iter);

  /** sanity checks */
  P4EST_ASSERT (mesh_iter->current_index <= mesh_iter->maximum_index);
  P4EST_ASSERT (-1 <= mesh_iter->current_vid
                && mesh_iter->current_vid < P4EST_CHILDREN);

  /** handle virtual quadrants */
  if (0 <= mesh_iter->current_vid && mesh_iter->current_vid < P4EST_CHILDREN) {
    ++mesh_iter->current_vid;
    ++mesh_iter->current_index;
  }

  /** handle if we are done processing a set of virtual quadrants */
  if (mesh_iter->current_vid == P4EST_CHILDREN) {
    mesh_iter->current_vid = -1;
    ++mesh_iter->processed_virtual;
    --mesh_iter->current_index;
  }

  /** if we are currently processing a set of virtual quadrants, we have
   * already found the next quadrant */
  if (0 < mesh_iter->current_vid)
    return 0;

  /** check if we have processed all local quadrants and switch on ghost
      if so */
  if ((mesh_iter->real_on_level_local + mesh_iter->virtual_on_level_local) <=
      (mesh_iter->processed_real + mesh_iter->processed_virtual)) {
    mesh_iter->current_is_ghost = 1;
  }

  /** Prepare to find a new quadrant id
   * There are 21 different cases from which set to choose the next quadrant:
   *  0.) all real local quadrants
   *  1.) real local quadrants at parallel boundary
   *  2.) real local inner quadrants
   *  3.) real or virtual local quadrants
   *  4.) real or virtual local quadrants at parallel boundary
   *  5.) real or virtual local inner quadrants
   *  6.) virtual local quadrants
   *  7.) virtual local quadrants at parallel boundary
   *  8.) virtual local inner quadrants
   *  9.) all real local or ghost quadrants
   * 10.) real local or ghost quadrants at parallel boundary
   * 11.) real local or ghost inner quadrants
   * 12.) real or virtual local or ghost quadrants
   * 13.) real or virtual local or ghost quadrants at parallel boundary
   * 14.) real or virtual local or ghost inner quadrants
   * 15.) virtual local or ghost quadrants
   * 16.) virtual local or ghost quadrants at parallel boundary
   * 17.) virtual local or ghost inner quadrants
   * 18.) real ghost quadrants
   * 19.) real or virtual ghost quadrants
   * 20.) virtual ghost quadrants
   *
   * As we traverse ghost quadrants after local quadrants if they are to be
   * traversed, there are at most 2 sets to choose from:
   * a) real and virtual local quadrants
   * b) real and virtual ghost quadrants
   * The constraints imposed by the parallel boundary will be imposed by
   * skipping the quadrants if they are not part of the respective set of
   * quadrants.
   *
   * To distinguish we create an encoding e, with
   *   e := 9 * trav_ghosts + 3 * trav_virtuals + trav_parallel_boundary
   * and map the cases 0 .. 17 sketched above to the respective encoding. The
   * cases 18, 19, and 20 are mapped to encoding ranges 18..20, 21..23, and
   * 24..26. The parallel_boundary flag is not needed here, because all ghost
   * quadrants are part of the parallel boundary.
   */
  int                 encode_case =
    9 * mesh_iter->traverse_ghosts + 3 * mesh_iter->traverse_virtuals +
    mesh_iter->traverse_parallel_boundary;

  do {
    ++mesh_iter->current_index;

    /** check if we are done */
    if (mesh_iter->maximum_index <= mesh_iter->current_index) {
      return P4EST_MESHITER_DONE;
    }

    /** Find new quadid */
    switch (encode_case) {
    case 0:                    /* local, real, inner + pb */
    case 1:                    /* local, real, pb */
    case 2:                    /* local, real, inner */
      /** next quadrant is next local real quadrant */
      mesh_iter->current_qid =
        *(int *) sc_array_index (mesh_iter->mesh->quad_level +
                                 mesh_iter->current_level,
                                 mesh_iter->processed_real);
      ++mesh_iter->processed_real;
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      break;
    case 3:
    case 4:
    case 5:
      /** next quadrant to choose from set a) */
      max_real =
        (mesh_iter->mesh->quad_level + mesh_iter->current_level)->elem_count;
      max_virt =
        (mesh_iter->mesh->virtual_qlevels +
         mesh_iter->current_level)->elem_count;
      if (mesh_iter->processed_real < max_real) {
        real_qid = *(int *) sc_array_index (mesh_iter->mesh->quad_level +
                                            mesh_iter->current_level,
                                            mesh_iter->processed_real);
      }
      else {
        real_qid = INT_MAX;
      }
      if (mesh_iter->processed_virtual < max_virt) {
        virt_qid = *(int *) sc_array_index (mesh_iter->mesh->virtual_qlevels +
                                            mesh_iter->current_level,
                                            mesh_iter->processed_virtual);
      }
      else {
        virt_qid = INT_MAX;
      }
      P4EST_ASSERT (virt_qid != real_qid);

      if (real_qid < virt_qid) {
        mesh_iter->current_qid = real_qid;
        ++mesh_iter->processed_real;
      }
      else {
        mesh_iter->current_qid = virt_qid;
        mesh_iter->current_vid = 0;
      }
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      else {
        mesh_iter->current_vid = -1;
      }
      break;
    case 6:
    case 7:
    case 8:
      /** next quadrant is next virtual quadrant */
      mesh_iter->current_qid =
        *(int *) sc_array_index (mesh_iter->mesh->virtual_qlevels +
                                 mesh_iter->current_level,
                                 mesh_iter->processed_virtual);
      mesh_iter->current_vid = 0;
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      else {
        mesh_iter->current_vid = -1;
      }
      break;
    case 9:
    case 10:
    case 11:
      /** while available: next quadrant is next local real quadrant, then next
          ghost real quadrant */
      if (!mesh_iter->current_is_ghost) {
        mesh_iter->current_qid =
          *(int *) sc_array_index (mesh_iter->mesh->quad_level +
                                   mesh_iter->current_level,
                                   mesh_iter->processed_real);
        ++mesh_iter->processed_real;
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER)
            ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY)
             && (1 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
      }
      else {
        mesh_iter->current_qid =
          *(int *) sc_array_index (mesh_iter->mesh->ghost_level +
                                   mesh_iter->current_level,
                                   mesh_iter->processed_real -
                                   mesh_iter->real_on_level_local);
        ++mesh_iter->processed_real;
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER)
            ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY)
             && (1 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
      }
      break;
    case 12:
    case 13:
    case 14:
      /** while available: next quadrant to choose from set a), then from set b) */
      if (!mesh_iter->current_is_ghost) {
        max_real =
          (mesh_iter->mesh->quad_level +
           mesh_iter->current_level)->elem_count;
        max_virt =
          (mesh_iter->mesh->virtual_qlevels +
           mesh_iter->current_level)->elem_count;
        if (mesh_iter->processed_real < max_real) {
          real_qid = *(int *) sc_array_index (mesh_iter->mesh->quad_level +
                                              mesh_iter->current_level,
                                              mesh_iter->processed_real);
        }
        else {
          real_qid = INT_MAX;
        }
        if (mesh_iter->processed_virtual < max_virt) {
          virt_qid =
            *(int *) sc_array_index (mesh_iter->mesh->virtual_qlevels +
                                     mesh_iter->current_level,
                                     mesh_iter->processed_virtual);
        }
        else {
          virt_qid = INT_MAX;
        }
        P4EST_ASSERT (real_qid != virt_qid);
        if (real_qid < virt_qid) {
          mesh_iter->current_qid = real_qid;
          ++mesh_iter->processed_real;
        }
        else {
          mesh_iter->current_qid = virt_qid;
          mesh_iter->current_vid = 0;
        }
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER)
            ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY)
             && (1 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
        else {
          mesh_iter->current_vid = -1;
        }
      }
      else {
        max_real =
          (mesh_iter->mesh->ghost_level +
           mesh_iter->current_level)->elem_count;
        max_virt =
          (mesh_iter->mesh->virtual_glevels +
           mesh_iter->current_level)->elem_count;
        if ((mesh_iter->processed_real - mesh_iter->real_on_level_local) <
            max_real) {
          real_qid =
            *(int *) sc_array_index (mesh_iter->mesh->ghost_level +
                                     mesh_iter->current_level,
                                     mesh_iter->processed_real -
                                     mesh_iter->real_on_level_local);
        }
        else {
          real_qid = INT_MAX;
        }
        if ((mesh_iter->processed_virtual -
             mesh_iter->virtual_on_level_local) < max_virt) {
          virt_qid =
            *(int *) sc_array_index (mesh_iter->mesh->virtual_glevels +
                                     mesh_iter->current_level,
                                     mesh_iter->processed_virtual -
                                     mesh_iter->virtual_on_level_local);
        }
        else {
          virt_qid = INT_MAX;
        }
        P4EST_ASSERT (real_qid != virt_qid);
        if (real_qid < virt_qid) {
          mesh_iter->current_qid = real_qid;
          ++mesh_iter->processed_real;
        }
        else {
          mesh_iter->current_qid = virt_qid;
          mesh_iter->current_vid = 0;
        }
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER)
            ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY)
             && (1 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
        else {
          mesh_iter->current_vid = -1;
        }
      }
      break;
    case 15:
    case 16:
    case 17:
      /** while available: next quadrant is next local virtual quadrant, then
       * next ghost virtual quadrant */
      if (!mesh_iter->current_is_ghost) {
        mesh_iter->current_qid =
          *(int *) sc_array_index (mesh_iter->mesh->virtual_qlevels +
                                   mesh_iter->current_level,
                                   mesh_iter->processed_virtual);
        mesh_iter->current_vid = 0;
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER) ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY) &&
             (1 ==
              mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
        else {
          mesh_iter->current_vid = -1;
        }
      }
      else {
        mesh_iter->current_qid =
          *(int *) sc_array_index (mesh_iter->mesh->virtual_glevels +
                                   mesh_iter->current_level,
                                   mesh_iter->processed_virtual -
                                   mesh_iter->virtual_on_level_local);
        mesh_iter->current_vid = 0;
        /** sanity checks for parallel boundary */
        if ((mesh_iter->traverse_parallel_boundary ==
             P4EST_TRAVERSE_PARBOUNDINNER)
            ||
            ((mesh_iter->traverse_parallel_boundary ==
              P4EST_TRAVERSE_PARALLEL_BOUNDARY)
             && (1 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
            ||
            ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
             && (0 ==
                 mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                    current_qid]))) {
          found_valid_quad = 1;
        }
        else {
          mesh_iter->current_vid = -1;
        }
      }
      break;
    case 18:
    case 19:
    case 20:
      /** next quadrant is next ghost real quadrant */
      mesh_iter->current_qid =
        *(int *) sc_array_index (mesh_iter->mesh->ghost_level +
                                 mesh_iter->current_level,
                                 mesh_iter->processed_real -
                                 mesh_iter->real_on_level_local);
      ++mesh_iter->processed_real;
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      break;
    case 21:
    case 22:
    case 23:
      /** next quadrant to choose from set b) */
      max_real =
        (mesh_iter->mesh->ghost_level + mesh_iter->current_level)->elem_count;
      max_virt =
        (mesh_iter->mesh->virtual_glevels +
         mesh_iter->current_level)->elem_count;
      if (mesh_iter->processed_real - mesh_iter->real_on_level_local <
          max_real) {

        real_qid = *(int *) sc_array_index (mesh_iter->mesh->ghost_level +
                                            mesh_iter->current_level,
                                            mesh_iter->processed_real -
                                            mesh_iter->real_on_level_local);
      }
      else {
        real_qid = INT_MAX;
      }
      if (mesh_iter->processed_virtual - mesh_iter->virtual_on_level_local <
          max_virt) {
        virt_qid =
          *(int *) sc_array_index (mesh_iter->mesh->virtual_glevels +
                                   mesh_iter->current_level,
                                   mesh_iter->processed_virtual -
                                   mesh_iter->virtual_on_level_local);
      }
      else {
        virt_qid = INT_MAX;
      }
      P4EST_ASSERT (real_qid != virt_qid);
      if (real_qid < virt_qid) {
        mesh_iter->current_qid = real_qid;
        ++mesh_iter->processed_real;
      }
      else {
        mesh_iter->current_qid = virt_qid;
        mesh_iter->current_vid = 0;
      }
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      else {
        mesh_iter->current_vid = -1;
      }
      break;
    case 24:
    case 25:
    case 26:
      /** next quadrant is next ghost virtual quadrant */
      mesh_iter->current_qid =
        *(int *) sc_array_index (mesh_iter->mesh->virtual_glevels +
                                 mesh_iter->current_level,
                                 mesh_iter->processed_virtual -
                                 mesh_iter->virtual_on_level_local);
      mesh_iter->current_vid = 0;
      /** sanity checks for parallel boundary */
      if ((mesh_iter->traverse_parallel_boundary ==
           P4EST_TRAVERSE_PARBOUNDINNER)
          ||
          ((mesh_iter->traverse_parallel_boundary ==
            P4EST_TRAVERSE_PARALLEL_BOUNDARY)
           && (1 ==
               mesh_iter->mesh->parallel_boundary[mesh_iter->current_qid]))
          || ((mesh_iter->traverse_parallel_boundary == P4EST_TRAVERSE_INNER)
              && (0 ==
                  mesh_iter->mesh->parallel_boundary[mesh_iter->
                                                     current_qid]))) {
        found_valid_quad = 1;
      }
      else {
        mesh_iter->current_vid = -1;
      }
      break;
    default:
      SC_ABORT_NOT_REACHED ();
    }
  } while (found_valid_quad == 0);

  return 0;
}

int
p4est_meshiter_set_neighbor_quad_info (p4est_meshiter_t * mesh_iter,
                                       int direction)
{
  int                 i;
  int                 n_neighbor_entities, n_hanging_quads;
  int                 l_same_size, u_same_size;
  int                 l_double_size, u_double_size;
  int                 l_half_size, u_half_size;

  /* reset previously set neighbor information */
  discard_neighbors (mesh_iter);

  /* set constants */
  set_limits (direction, mesh_iter->btype, &l_same_size, &u_same_size,
              &l_double_size, &u_double_size, &l_half_size, &u_half_size,
              &n_neighbor_entities, &n_hanging_quads);

  /** set direction */
  mesh_iter->neighbor_direction = direction;

  /** current quadrant is a real quadrant */
  int                 curr_neighbor_enc;

  p4est_mesh_get_neighbors (mesh_iter->p4est, mesh_iter->ghost,
                            mesh_iter->mesh, mesh_iter->current_qid,
                            mesh_iter->current_vid, direction, 1, NULL,
                            &mesh_iter->neighbor_enc,
                            &mesh_iter->neighbor_ids);

  P4EST_ASSERT ((0 == mesh_iter->neighbor_ids.elem_count) ||
                (2 * mesh_iter->neighbor_enc.elem_count ==
                 mesh_iter->neighbor_ids.elem_count));
  if (0 == mesh_iter->neighbor_enc.elem_count) {
    discard_neighbors (mesh_iter);
  }
  else {
    P4EST_ASSERT (1 == mesh_iter->neighbor_enc.elem_count);

    curr_neighbor_enc = *(int *) sc_array_index (&mesh_iter->neighbor_enc, 0);
    mesh_iter->neighbor_is_ghost = (curr_neighbor_enc < 0);

    /* normalize encoding */
    curr_neighbor_enc = abs (curr_neighbor_enc) - 1;
    p4est_meshiter_decode_encoding (mesh_iter, curr_neighbor_enc,
                                    n_neighbor_entities, l_same_size,
                                    u_same_size, l_double_size, u_double_size,
                                    l_half_size, u_half_size);

    mesh_iter->neighbor_qid =
      *(int *) sc_array_index (&mesh_iter->neighbor_ids, 0);
    mesh_iter->neighbor_vid =
      *(int *) sc_array_index (&mesh_iter->neighbor_ids, 1);
  }

  return 0;
}

int
p4est_meshiter_get_current_storage_id (p4est_meshiter_t * mesh_iter)
{
  int                 sid, max_id;
  max_id = mesh_iter->current_is_ghost ?
    (mesh_iter->mesh->ghost_level + mesh_iter->current_level)->elem_count +
    P4EST_CHILDREN * (mesh_iter->mesh->virtual_glevels +
                      mesh_iter->current_level)->elem_count :
    (mesh_iter->mesh->quad_level + mesh_iter->current_level)->elem_count +
    P4EST_CHILDREN * (mesh_iter->mesh->virtual_qlevels +
                      mesh_iter->current_level)->elem_count;

  if (-1 < mesh_iter->current_vid) {
    P4EST_ASSERT (mesh_iter->traverse_virtuals != P4EST_TRAVERSE_REAL);
    P4EST_ASSERT (0 <= mesh_iter->current_vid
                  && mesh_iter->current_vid < P4EST_CHILDREN);

    /** virtual */
    if (!mesh_iter->current_is_ghost) {
      /* local */
      P4EST_ASSERT (mesh_iter->traverse_ghosts != P4EST_TRAVERSE_GHOST);
      sid = mesh_iter->mesh->quad_qvirtual_offset[mesh_iter->current_qid] +
        mesh_iter->current_vid;
    }
    else {
      /** ghost */
      P4EST_ASSERT (mesh_iter->traverse_ghosts != P4EST_TRAVERSE_LOCAL);
      sid = mesh_iter->mesh->quad_gvirtual_offset[mesh_iter->current_qid] +
        mesh_iter->current_vid;
    }
  }
  else {
    P4EST_ASSERT (mesh_iter->traverse_virtuals != P4EST_TRAVERSE_VIRTUAL);
    P4EST_ASSERT (mesh_iter->current_vid == -1);
    /** real */
    if (!mesh_iter->current_is_ghost) {
      /** local */
      P4EST_ASSERT (mesh_iter->traverse_ghosts != P4EST_TRAVERSE_GHOST);
      sid = mesh_iter->mesh->quad_qreal_offset[mesh_iter->current_qid];
    }
    else {
      /** ghost */
      P4EST_ASSERT (mesh_iter->traverse_ghosts != P4EST_TRAVERSE_LOCAL);
      sid = mesh_iter->mesh->quad_greal_offset[mesh_iter->current_qid];
    }
  }
  P4EST_ASSERT (0 <= sid && sid < max_id);

  return sid;
}

int
p4est_meshiter_get_neighbor_storage_id (p4est_meshiter_t * mesh_iter)
{
  int                 sid, max_id;
  max_id = mesh_iter->neighbor_is_ghost ?
    (mesh_iter->mesh->ghost_level + mesh_iter->current_level)->elem_count +
    P4EST_CHILDREN * (mesh_iter->mesh->virtual_glevels +
                      mesh_iter->current_level)->elem_count :
    (mesh_iter->mesh->quad_level + mesh_iter->current_level)->elem_count +
    P4EST_CHILDREN * (mesh_iter->mesh->virtual_qlevels +
                      mesh_iter->current_level)->elem_count;

  if (-1 < mesh_iter->neighbor_vid) {
    P4EST_ASSERT (mesh_iter->traverse_virtuals != P4EST_TRAVERSE_REAL);
    P4EST_ASSERT (0 <= mesh_iter->neighbor_vid
                  && mesh_iter->neighbor_vid < P4EST_CHILDREN);

    /** virtual */
    if (!mesh_iter->neighbor_is_ghost) {
      /** local */
      sid = mesh_iter->mesh->quad_qvirtual_offset[mesh_iter->neighbor_qid] +
        mesh_iter->neighbor_vid;
    }
    else {
      /** ghost */
      sid = mesh_iter->mesh->quad_gvirtual_offset[mesh_iter->neighbor_qid] +
        mesh_iter->neighbor_vid;
    }
  }
  else {
    P4EST_ASSERT (mesh_iter->traverse_virtuals != P4EST_TRAVERSE_VIRTUAL);
    P4EST_ASSERT (mesh_iter->neighbor_vid == -1);
    /** real */
    if (!mesh_iter->neighbor_is_ghost) {
      /** local */
      sid = mesh_iter->mesh->quad_qreal_offset[mesh_iter->neighbor_qid];
    }
    else {
      /** ghost */
      sid = mesh_iter->mesh->quad_greal_offset[mesh_iter->neighbor_qid];
    }
  }
  P4EST_ASSERT (0 <= sid && sid < max_id);

  return sid;
}
