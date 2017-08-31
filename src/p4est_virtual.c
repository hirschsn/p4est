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
#include <p4est_virtual.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#else /* !P4_TO_P8 */
#include <p8est_virtual.h>
#include <p8est_extended.h>
#endif /* !P4_TO_P8 */

#ifndef P4_TO_P8

/* *INDENT-OFF* */
const int           p4est_face_virtual_neighbors_inside[P4EST_CHILDREN]
                                                       [P4EST_FACES] =
{{  4,  1,  6,  2 },
 {  0,  5, 10,  3 },
 {  8,  3,  0,  7 },
 {  2,  9,  1, 11 }};

const int           p4est_corner_virtual_neighbors_inside[P4EST_CHILDREN]
                                                         [P4EST_CHILDREN] =
{{ 12, 10,  8,  3 },
 {  6, 13,  2,  9 },
 {  4,  1, 14, 11 },
 {  0,  5,  7, 15 }};
/* *INDENT-ON* */

#endif /* P4_TO_P8 */

/** Determine if qid needs to contain virtual quadrants for inner quadrants,
 * i.e. quadrants that are not mirrors.
 * This function can potentially exit the loop earlier than \ref
 * has_virtuals_parallel_boundary, because it needs not decide if ghost
 * quadrants need virtual quadrants as well.
 * For using this optimization mesh needs a populated parallel_boundary array.
 * \param    [out] virtual   Virtual structure to populate.
 * \param[in]      ghost     Current ghost-layer.
 * \param[in]      mesh      Mesh structure containing neighbor information.
 * \param[in]      qid       Local quadrant index for which to check.
 * \param[in][out] lq_per_level_real   Number of real quadrants per level
 *                                     processed up to current quadrant.
 * \param[in][out] lq_per_level_virt   Number of virtual quadrants per level
                                       created up to current quadrant.
 * \param[in]      quads     Empty container for more efficient neighbor search,
 *                           contains *p4est_quadrant_t.
 */
static int
has_virtuals_inner (p4est_virtual_t * virtual, p4est_t * p4est,
                    p4est_ghost_t * ghost, p4est_mesh_t * mesh, int qid,
                    p4est_locidx_t * lq_per_level_real,
                    p4est_locidx_t * lq_per_level_virt, int *last_virtual,
                    sc_array_t * quads)
{
  int                 i, imax, j;
  int                 has_virtuals = 0;
  p4est_quadrant_t   *curr_quad = p4est_mesh_get_quadrant (p4est, mesh, qid);
  p4est_quadrant_t   *neighbor;
  int                 level = curr_quad->level;
  p4est_locidx_t     *insert_locidx_t;

  switch (virtual->btype) {
  case P4EST_CONNECT_FACE:
    imax = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    imax = P4EST_FACES + P8EST_EDGES;
    break;
#endif /* P4_TO_P8 */
  case P4EST_CONNECT_FULL:
    /* *INDENT-OFF* */
    imax = P4EST_FACES +
#ifdef P4_TO_P8
           P8EST_EDGES +
#endif /* P4_TO_P8 */
           P4EST_CHILDREN;
    /* *INDENT-ON* */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  for (has_virtuals = 0, i = 0; !has_virtuals && i < imax; ++i) {
    sc_array_truncate (quads);

    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid, i, quads, NULL, NULL);
    for (j = 0; j < quads->elem_count; ++j) {
      neighbor = *(p4est_quadrant_t **) sc_array_index (quads, j);
      if (level < neighbor->level) {
        has_virtuals = 1;
        break;
      }
    }
  }
  virtual->quad_qreal_offset[qid] =
    lq_per_level_real[level] + P4EST_CHILDREN * lq_per_level_virt[level];
  ++lq_per_level_real[level];
  if (has_virtuals) {
    virtual->virtual_qflags[qid] = ++(*last_virtual);
    virtual->quad_qvirtual_offset[qid] =
      lq_per_level_real[level + 1] +
      P4EST_CHILDREN * lq_per_level_virt[level + 1];
    insert_locidx_t =
      (p4est_locidx_t *) sc_array_push (virtual->virtual_qlevels +
                                        (level + 1));
    *insert_locidx_t = qid;
  }

  return 0;
}

/** Determine if qid needs to contain virtual quadrants for quadrants that are
 * either mirrors or for a mesh without parallel_boundary array.
 * This function always checks all neighbors, because it has to decide if ghost
 * quadrants need virtual quadrants.
 * \param    [out] virtual   Virtual structure to populate.
 * \param[in]      ghost     Current ghost-layer.
 * \param[in]      mesh      Mesh structure containing neighbor information.
 * \param[in]      qid       Local quadrant index for which to check.
 * \param[in][out] lq_per_level_real   Number of real quadrants per level
 *                                     processed up to current quadrant.
 * \param[in][out] lq_per_level_virt   Number of virtual quadrants per level
                                       created up to current quadrant.
 * \param[in]      quads     Empty container for more efficient neighbor search,
 *                           contains *p4est_quadrant_t.
 * \param[in]      qids      Empty container for more efficient neighbor search,
 *                           contains p4est_locidx_t.
 */
static int
has_virtuals_parallel_boundary (p4est_virtual_t * virtual, p4est_t * p4est,
                                p4est_ghost_t * ghost, p4est_mesh_t * mesh,
                                p4est_locidx_t qid,
                                p4est_locidx_t * lq_per_level_real,
                                p4est_locidx_t * lq_per_level_virt,
                                int *last_virtual, sc_array_t * quads,
                                sc_array_t * qids)
{
  int                 i, imax, j;
  p4est_locidx_t      lq, gq;
  int                 level;
  int                 has_virtuals = 0;
  p4est_quadrant_t   *curr_quad = p4est_mesh_get_quadrant (p4est, mesh, qid);
  p4est_quadrant_t   *neighbor;
  p4est_locidx_t      neighbor_qid;
  p4est_locidx_t     *insert_locidx_t;

  lq = virtual->local_num_quadrants;
  gq = virtual->ghost_num_quadrants;

  level = curr_quad->level;

  switch (virtual->btype) {
  case P4EST_CONNECT_FACE:
    imax = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    imax = P4EST_FACES + P8EST_EDGES;
#endif /* P4_TO_P8 */
    break;
  case P4EST_CONNECT_FULL:
    /* *INDENT-OFF* */
    imax = P4EST_FACES +
#ifdef P4_TO_P8
           P8EST_EDGES +
#endif /* P4_TO_P8 */
           P4EST_CHILDREN;
    /* *INDENT-ON* */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* check if virtual quadrants need to be created */
  for (i = 0; i < imax; ++i) {
    sc_array_truncate (quads);
    sc_array_truncate (qids);

    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid, i, quads, NULL, qids);
    for (j = 0; j < quads->elem_count; ++j) {
      neighbor = *(p4est_quadrant_t **) sc_array_index (quads, j);
      neighbor_qid = *(int *) sc_array_index (qids, j);
      if (level < neighbor->level) {
        has_virtuals = 1;
      }
      else if ((lq <= neighbor_qid) && (neighbor_qid <= (lq + gq))
               && (neighbor->level < level)) {
        virtual->virtual_gflags[neighbor_qid - lq] = 1;
      }
    }
  }
  virtual->quad_qreal_offset[qid] =
    lq_per_level_real[level] + P4EST_CHILDREN * lq_per_level_virt[level];
  ++lq_per_level_real[level];
  if (has_virtuals) {
    virtual->virtual_qflags[qid] = ++(*last_virtual);
    virtual->quad_qvirtual_offset[qid] =
      lq_per_level_real[level + 1] +
      P4EST_CHILDREN * lq_per_level_virt[level + 1];
    insert_locidx_t =
      (p4est_locidx_t *) sc_array_push (virtual->virtual_qlevels +
                                        (level + 1));
    *insert_locidx_t = qid;
  }

  return 0;
}

p4est_virtual_t    *
p4est_virtual_new (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_mesh_t * mesh, p4est_connect_type_t btype)
{
  return p4est_virtual_new_ext (p4est, ghost, mesh, btype, 0);
}

p4est_virtual_t    *
p4est_virtual_new_ext (p4est_t * p4est, p4est_ghost_t * ghost,
                       p4est_mesh_t * mesh, p4est_connect_type_t btype,
                       int compute_level_lists)
{
  p4est_virtual_t    *virtual;
  int                 last_virtual_index = -1;
  p4est_locidx_t      level, quad;
  sc_array_t         *quads, *qids;
  p4est_locidx_t      lq, gq;
  p4est_locidx_t     *lq_per_level_real, *lq_per_level_virt;
  p4est_locidx_t     *gq_per_level_real, *gq_per_level_virt;
  p4est_locidx_t     *insert_locidx_t;
  p4est_quadrant_t   *ghost_quad;

  quads = sc_array_new (sizeof (p4est_quadrant_t *));
  qids = sc_array_new (sizeof (int));
  lq_per_level_real = P4EST_ALLOC_ZERO (p4est_locidx_t, P4EST_QMAXLEVEL + 1);
  lq_per_level_virt = P4EST_ALLOC_ZERO (p4est_locidx_t, P4EST_QMAXLEVEL + 1);
  gq_per_level_real = P4EST_ALLOC_ZERO (p4est_locidx_t, P4EST_QMAXLEVEL + 1);
  gq_per_level_virt = P4EST_ALLOC_ZERO (p4est_locidx_t, P4EST_QMAXLEVEL + 1);

  /* check if input conditions are met */
  P4EST_ASSERT (p4est_is_balanced (p4est, btype));
  P4EST_ASSERT (btype <= mesh->btype);

  virtual = P4EST_ALLOC_ZERO (p4est_virtual_t, 1);

  virtual->btype = btype;
  virtual->local_num_quadrants = lq = mesh->local_num_quadrants;
  virtual->ghost_num_quadrants = gq = mesh->ghost_num_quadrants;
  virtual->virtual_qflags = P4EST_ALLOC (p4est_locidx_t, lq);
  virtual->virtual_gflags = P4EST_ALLOC (p4est_locidx_t, gq);
  memset (virtual->virtual_qflags, (char) -1, lq * sizeof (p4est_locidx_t));
  memset (virtual->virtual_gflags, (char) -1, gq * sizeof (p4est_locidx_t));

  if (compute_level_lists) {
    virtual->quad_qreal_offset = P4EST_ALLOC (p4est_locidx_t, lq);
    virtual->quad_qvirtual_offset = P4EST_ALLOC (p4est_locidx_t, lq);
    virtual->quad_greal_offset = P4EST_ALLOC (p4est_locidx_t, gq);
    virtual->quad_gvirtual_offset = P4EST_ALLOC (p4est_locidx_t, gq);
    memset (virtual->quad_qreal_offset, (char) -1,
            lq * sizeof (p4est_locidx_t));
    memset (virtual->quad_qvirtual_offset, (char) -1,
            lq * sizeof (p4est_locidx_t));
    memset (virtual->quad_greal_offset, (char) -1,
            gq * sizeof (p4est_locidx_t));
    memset (virtual->quad_gvirtual_offset, (char) -1,
            gq * sizeof (p4est_locidx_t));

    virtual->virtual_qlevels = P4EST_ALLOC (sc_array_t, P4EST_QMAXLEVEL + 1);
    virtual->virtual_glevels = P4EST_ALLOC (sc_array_t, P4EST_QMAXLEVEL + 1);
    for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
      sc_array_init (virtual->virtual_qlevels + level,
                     sizeof (p4est_locidx_t));
      sc_array_init (virtual->virtual_glevels + level,
                     sizeof (p4est_locidx_t));
    }
  }

  for (quad = 0; quad < lq; ++quad) {
    sc_array_truncate (quads);
    sc_array_truncate (qids);
    if (mesh->parallel_boundary && -1 == mesh->parallel_boundary[quad]) {
      has_virtuals_inner (virtual, p4est, ghost, mesh, quad,
                          lq_per_level_real, lq_per_level_virt,
                          &last_virtual_index, quads);
    }
    else {
      has_virtuals_parallel_boundary (virtual, p4est, ghost, mesh, quad,
                                      lq_per_level_real, lq_per_level_virt,
                                      &last_virtual_index, quads, qids);
    }
  }

  last_virtual_index = 0;
  /* set gflags and create level and offset arrays if necessary */
  for (quad = 0; quad < gq; ++quad) {
    if (virtual->quad_qreal_offset) {
      ghost_quad = p4est_quadrant_array_index (&ghost->ghosts, quad);
      level = ghost_quad->level;

      virtual->quad_greal_offset[quad] =
        gq_per_level_real[level] + P4EST_CHILDREN * gq_per_level_virt[level];
      ++gq_per_level_real[level];
    }
    if (virtual->virtual_gflags[quad] != -1) {
      virtual->virtual_gflags[quad] = last_virtual_index;
      ++last_virtual_index;
      if (virtual->quad_qreal_offset) {
        virtual->quad_gvirtual_offset[quad] =
          gq_per_level_real[level + 1] +
          P4EST_CHILDREN * gq_per_level_virt[level + 1];
        ++gq_per_level_virt[level + 1];
        insert_locidx_t =
          (p4est_locidx_t *) sc_array_push (virtual->virtual_glevels +
                                            (level + 1));
        *insert_locidx_t = quad;
      }
    }
  }

  P4EST_FREE (lq_per_level_real);
  P4EST_FREE (lq_per_level_virt);
  P4EST_FREE (gq_per_level_real);
  P4EST_FREE (gq_per_level_virt);

  sc_array_destroy (quads);
  sc_array_destroy (qids);

  return virtual;
}

void
p4est_virtual_destroy (p4est_virtual_t * virtual)
{
  int                 i;

  P4EST_FREE (virtual->virtual_qflags);
  P4EST_FREE (virtual->virtual_gflags);
  if (virtual->quad_qreal_offset != NULL) {
    P4EST_FREE (virtual->quad_qreal_offset);
    P4EST_FREE (virtual->quad_qvirtual_offset);
    P4EST_FREE (virtual->quad_greal_offset);
    P4EST_FREE (virtual->quad_gvirtual_offset);

    for (i = 0; i < P4EST_QMAXLEVEL + 1; ++i) {
      sc_array_reset (virtual->virtual_qlevels + i);
      sc_array_reset (virtual->virtual_glevels + i);
    }
    P4EST_FREE (virtual->virtual_qlevels);
    P4EST_FREE (virtual->virtual_glevels);
  }
  P4EST_FREE (virtual);
}

size_t
p4est_virtual_memory_used (p4est_virtual_t * virtual)
{
  size_t              lqz, ngz;
  int                 level;
  size_t              mem_flags = 0;
  size_t              mem_offset = 0;
  size_t              mem_levels = 0;
  size_t              all_mem;

  lqz = (size_t) virtual->local_num_quadrants;
  ngz = (size_t) virtual->ghost_num_quadrants;

  mem_flags = (lqz + ngz) * sizeof (p4est_locidx_t);
  if (virtual->quad_qreal_offset) {
    mem_offset = 2 * (lqz + ngz) * sizeof (p4est_locidx_t);
    mem_levels = 2 * sizeof (sc_array_t) * (P4EST_QMAXLEVEL + 1);
    for (level = 0; level <= P4EST_QMAXLEVEL; ++level) {
      mem_levels +=
        sc_array_memory_used (virtual->virtual_qlevels + level, 0);
      mem_levels +=
        sc_array_memory_used (virtual->virtual_glevels + level, 0);
    }
  }

  all_mem = mem_flags + mem_offset + mem_levels + sizeof (p4est_virtual_t);

  return all_mem;
}
