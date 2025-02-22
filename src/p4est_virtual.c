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
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#else /* !P4_TO_P8 */
#include <p8est_virtual.h>
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#endif /* !P4_TO_P8 */

/* -------------------------------------------------------------------------- */
/* |                           Virtual quadrants                            | */
/* -------------------------------------------------------------------------- */
/** Determine if qid needs to contain virtual quadrants.
 * \param    [out] virtual_quads       Virtual structure to populate.
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
has_virtuals (p4est_virtual_t * virtual_quads, p4est_t * p4est,
              p4est_ghost_t * ghost, p4est_mesh_t * mesh, int qid,
              p4est_locidx_t * lq_per_level_real,
              p4est_locidx_t * lq_per_level_virt, int *last_virtual,
              sc_array_t * quads, sc_array_t * qids)
{
  int                 i, imax, j;
  int                 has_virtuals = 0;
  p4est_quadrant_t   *curr_quad = p4est_mesh_get_quadrant (p4est, mesh, qid);
  p4est_quadrant_t   *neighbor;
  int                 level = curr_quad->level;
  p4est_locidx_t     *insert_locidx_t, neighbor_qid;

  switch (virtual_quads->btype) {
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

  for (has_virtuals = virtual_quads->virtual_qflags[qid], i = 0;
       has_virtuals == -1 && i < imax; ++i) {
    sc_array_truncate (quads);
    sc_array_truncate (qids);
    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid, i, quads, NULL, qids);

    for (j = 0; j < (int) quads->elem_count; ++j) {
      neighbor = *(p4est_quadrant_t **) sc_array_index (quads, j);
      if (level < neighbor->level) {
        has_virtuals = 1;
        break;
      }
      else if (neighbor->level < level) {
        neighbor_qid = *(int *) sc_array_index (qids, j);
        if (neighbor_qid < virtual_quads->local_num_quadrants &&
            qid < neighbor_qid) {
          virtual_quads->virtual_qflags[neighbor_qid] = 1;
        }
      }
    }
  }
  if (virtual_quads->quad_qreal_offset) {
    virtual_quads->quad_qreal_offset[qid] =
      lq_per_level_real[level] + P4EST_CHILDREN * lq_per_level_virt[level];
    ++lq_per_level_real[level];
  }
  if (has_virtuals != -1) {
    virtual_quads->virtual_qflags[qid] = ++(*last_virtual);
    virtual_quads->local_num_virtuals += P4EST_CHILDREN;
    if (virtual_quads->quad_qreal_offset) {
      virtual_quads->quad_qvirtual_offset[qid] =
        lq_per_level_real[level + 1] +
        P4EST_CHILDREN * lq_per_level_virt[level + 1];
      ++lq_per_level_virt[level + 1];
      insert_locidx_t =
        (p4est_locidx_t *) sc_array_push (virtual_quads->virtual_qlevels +
                                          (level + 1));
      *insert_locidx_t = qid;
    }
  }

  return 0;
}

static void
populate_ghost_flags(p4est_virtual_t *virtual_quads, p4est_t *p4est,
                     p4est_ghost_t *ghost, p4est_mesh_t *mesh,
                     sc_array_t *quads, sc_array_t *qids)
{
  int i, imax, j;
  int level, mirror;
  p4est_locidx_t mirror_qid, neighbor_qid;
  p4est_quadrant_t *curr_quad, *neighbor;
  p4est_locidx_t lq = virtual_quads->local_num_quadrants;
  p4est_locidx_t gq = virtual_quads->ghost_num_quadrants;

  switch (virtual_quads->btype) {
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

  for (mirror = 0; mirror < ghost->mirrors.elem_count; ++mirror) {
    mirror_qid = mesh->mirror_qid[mirror];
    curr_quad = p4est_mesh_get_quadrant(p4est, mesh, mirror_qid);
    level = curr_quad->level;
    for (i = 0; i < imax; ++i) {
      sc_array_truncate (quads);
      sc_array_truncate (qids);
      p4est_mesh_get_neighbors(p4est, ghost, mesh, mirror_qid, i, quads, NULL,
                               qids);

      for (j = 0; j < (int) quads->elem_count; ++j) {
        neighbor_qid = *(int*) sc_array_index(qids, j);
        neighbor_qid -= lq;
        if (0 <= neighbor_qid && neighbor_qid < gq) {
          neighbor = *(p4est_quadrant_t **) sc_array_index (quads, j);
          if (neighbor->level < level) {
            virtual_quads->virtual_gflags[neighbor_qid] = 1;
          }
        }
      }
    }
  }
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
  p4est_virtual_t    *virtual_quads;
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

  virtual_quads = P4EST_ALLOC_ZERO (p4est_virtual_t, 1);

  virtual_quads->btype = btype;
  virtual_quads->local_num_quadrants = lq = mesh->local_num_quadrants;
  virtual_quads->ghost_num_quadrants = gq = mesh->ghost_num_quadrants;
  virtual_quads->virtual_qflags = P4EST_ALLOC (p4est_locidx_t, lq);
  virtual_quads->virtual_gflags = P4EST_ALLOC (p4est_locidx_t, gq);
  memset (virtual_quads->virtual_qflags, (char) -1,
          lq * sizeof (p4est_locidx_t));
  memset (virtual_quads->virtual_gflags, (char) -1,
          gq * sizeof (p4est_locidx_t));

  if (compute_level_lists) {
    virtual_quads->quad_qreal_offset = P4EST_ALLOC (p4est_locidx_t, lq);
    virtual_quads->quad_qvirtual_offset = P4EST_ALLOC (p4est_locidx_t, lq);
    virtual_quads->quad_greal_offset = P4EST_ALLOC (p4est_locidx_t, gq);
    virtual_quads->quad_gvirtual_offset = P4EST_ALLOC (p4est_locidx_t, gq);
    memset (virtual_quads->quad_qreal_offset, (char) -1,
            lq * sizeof (p4est_locidx_t));
    memset (virtual_quads->quad_qvirtual_offset, (char) -1,
            lq * sizeof (p4est_locidx_t));
    memset (virtual_quads->quad_greal_offset, (char) -1,
            gq * sizeof (p4est_locidx_t));
    memset (virtual_quads->quad_gvirtual_offset, (char) -1,
            gq * sizeof (p4est_locidx_t));

    virtual_quads->virtual_qlevels =
      P4EST_ALLOC (sc_array_t, P4EST_QMAXLEVEL + 1);
    virtual_quads->virtual_glevels =
      P4EST_ALLOC (sc_array_t, P4EST_QMAXLEVEL + 1);
    for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
      sc_array_init (virtual_quads->virtual_qlevels + level,
                     sizeof (p4est_locidx_t));
      sc_array_init (virtual_quads->virtual_glevels + level,
                     sizeof (p4est_locidx_t));
    }
  }

  for (quad = 0; quad < lq; ++quad) {
    has_virtuals (virtual_quads, p4est, ghost, mesh, quad, lq_per_level_real,
                  lq_per_level_virt, &last_virtual_index, quads, qids);
  }

  populate_ghost_flags(virtual_quads, p4est, ghost, mesh, quads, qids);

  last_virtual_index = 0;
  /* set gflags and create level and offset arrays if necessary */
  for (quad = 0; quad < gq; ++quad) {
    ghost_quad = p4est_quadrant_array_index (&ghost->ghosts, quad);
    level = ghost_quad->level;
    if (virtual_quads->quad_qreal_offset) {
      virtual_quads->quad_greal_offset[quad] =
        gq_per_level_real[level] + P4EST_CHILDREN * gq_per_level_virt[level];
      ++gq_per_level_real[level];
    }
    if (virtual_quads->virtual_gflags[quad] != -1) {
      virtual_quads->virtual_gflags[quad] = last_virtual_index;
      ++last_virtual_index;
      virtual_quads->ghost_num_virtuals += P4EST_CHILDREN;
      if (virtual_quads->quad_qreal_offset) {
        virtual_quads->quad_gvirtual_offset[quad] =
          gq_per_level_real[level + 1] +
          P4EST_CHILDREN * gq_per_level_virt[level + 1];
        ++gq_per_level_virt[level + 1];
        insert_locidx_t =
          (p4est_locidx_t *) sc_array_push (virtual_quads->virtual_glevels +
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

  return virtual_quads;
}

void
p4est_virtual_destroy (p4est_virtual_t * virtual_quads)
{
  int                 i;

  P4EST_FREE (virtual_quads->virtual_qflags);
  P4EST_FREE (virtual_quads->virtual_gflags);
  if (virtual_quads->quad_qreal_offset != NULL) {
    P4EST_FREE (virtual_quads->quad_qreal_offset);
    P4EST_FREE (virtual_quads->quad_qvirtual_offset);
    P4EST_FREE (virtual_quads->quad_greal_offset);
    P4EST_FREE (virtual_quads->quad_gvirtual_offset);

    for (i = 0; i < P4EST_QMAXLEVEL + 1; ++i) {
      sc_array_reset (virtual_quads->virtual_qlevels + i);
      sc_array_reset (virtual_quads->virtual_glevels + i);
    }
    P4EST_FREE (virtual_quads->virtual_qlevels);
    P4EST_FREE (virtual_quads->virtual_glevels);
  }
  P4EST_FREE (virtual_quads);
}

size_t
p4est_virtual_memory_used (p4est_virtual_t * virtual_quads)
{
  size_t              lqz, ngz;
  int                 level;
  size_t              mem_flags = 0;
  size_t              mem_offset = 0;
  size_t              mem_levels = 0;
  size_t              all_mem;

  lqz = (size_t) virtual_quads->local_num_quadrants;
  ngz = (size_t) virtual_quads->ghost_num_quadrants;

  mem_flags = (lqz + ngz) * sizeof (p4est_locidx_t);
  if (virtual_quads->quad_qreal_offset) {
    mem_offset = 2 * (lqz + ngz) * sizeof (p4est_locidx_t);
    mem_levels = 2 * sizeof (sc_array_t) * (P4EST_QMAXLEVEL + 1);
    for (level = 0; level <= P4EST_QMAXLEVEL; ++level) {
      mem_levels +=
        sc_array_memory_used (virtual_quads->virtual_qlevels + level, 0);
      mem_levels +=
        sc_array_memory_used (virtual_quads->virtual_glevels + level, 0);
    }
  }

  all_mem = mem_flags + mem_offset + mem_levels + sizeof (p4est_virtual_t);

  return all_mem;
}

/* -------------------------------------------------------------------------- */
/* |                             Ghost exchange                             | */
/* -------------------------------------------------------------------------- */
p4est_virtual_ghost_t *
p4est_virtual_ghost_new (p4est_t * p4est, p4est_ghost_t * ghost,
                         p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                         p4est_connect_type_t btype)
{
  int                 proc;
  p4est_locidx_t      lq, gq;
  p4est_locidx_t      offset_begin, offset_end;
  p4est_locidx_t      mirror_idx, mirror_qid;
  p4est_locidx_t      neighbor_qid;
  p4est_quadrant_t   *curr_quad, *neighbor_quad;
  int                 n, neighbor_idx, max_neighbor_idx;
  sc_array_t         *nqid, *nquad;
  p4est_virtual_ghost_t *virtual_ghost;

  lq = mesh->local_num_quadrants;
  gq = mesh->ghost_num_quadrants;
  nqid = sc_array_new (sizeof (p4est_locidx_t));
  nquad = sc_array_new (sizeof (p4est_quadrant_t *));

  virtual_ghost = P4EST_ALLOC_ZERO (p4est_virtual_ghost_t, 1);
  P4EST_ASSERT (btype <= virtual_quads->btype);
  virtual_ghost->btype = btype;
  virtual_ghost->mirror_proc_virtuals =
    P4EST_ALLOC_ZERO (int8_t, ghost->mirror_proc_offsets[p4est->mpisize]);

  switch (btype) {
  case P4EST_CONNECT_FACE:
    max_neighbor_idx = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    max_neighbor_idx = P4EST_FACES + P8EST_EDGES;
    break;
#endif /* P4_TO_P8 */
  case P4EST_CONNECT_FULL:
      /* *INDENT-OFF* */
      max_neighbor_idx = P4EST_FACES +
#ifdef P4_TO_P8
                         P8EST_EDGES +
#endif /* P4_TO_P8 */
                         P4EST_CHILDREN;
      /* *INDENT-ON* */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* populate mirror_proc_virtuals:
   * Iterate ghost->mirror_proc_mirrors for each process.  Consider for each
   * mirror index hosting virtual quadrants if its ghost neighbors on the
   * respective neighbor rank are half-sized w.r.t. that mirror.  In this case
   * the neighboring process places virtual quads which we have to send.
   */
  for (proc = 0; proc < p4est->mpisize; ++proc) {
    offset_begin = ghost->mirror_proc_offsets[proc];
    offset_end = ghost->mirror_proc_offsets[proc + 1];
    for (mirror_idx = offset_begin; mirror_idx < offset_end; ++mirror_idx) {
      mirror_qid = mesh->mirror_qid[ghost->mirror_proc_mirrors[mirror_idx]];
      curr_quad = p4est_mesh_get_quadrant (p4est, mesh, mirror_qid);
      if (-1 < virtual_quads->virtual_qflags[mirror_qid]) {
        for (neighbor_idx = 0;
             virtual_ghost->mirror_proc_virtuals[mirror_idx] != 1
             && neighbor_idx < max_neighbor_idx; ++neighbor_idx) {
          sc_array_truncate (nqid);
          sc_array_truncate (nquad);
          p4est_mesh_get_neighbors (p4est, ghost, mesh, mirror_qid,
                                    neighbor_idx, nquad, NULL, nqid);
          for (n = 0; n < nqid->elem_count; ++n) {
            neighbor_qid = *(p4est_locidx_t *) sc_array_index (nqid, n);
            if (lq <= neighbor_qid && neighbor_qid < (lq + gq)) {
              neighbor_qid -= lq;
              if (mesh->ghost_to_proc[neighbor_qid] == proc) {
                neighbor_quad =
                  *(p4est_quadrant_t **) sc_array_index (nquad, n);
                if (neighbor_quad->level > curr_quad->level) {
                  virtual_ghost->mirror_proc_virtuals[mirror_idx] = 1;
                }
              }
            }
          }
        }
      }
    }
  }

  sc_array_destroy(nqid);
  sc_array_destroy(nquad);

  return virtual_ghost;
}

void
p4est_virtual_ghost_destroy (p4est_virtual_ghost_t * virtual_ghost)
{
  P4EST_FREE (virtual_ghost->mirror_proc_virtuals);
  P4EST_FREE (virtual_ghost);
}

size_t
p4est_virtual_ghost_memory_used (p4est_virtual_ghost_t * virtual_ghost)
{
  return 0;
}

void
p4est_virtual_ghost_exchange_data (p4est_t * p4est, p4est_ghost_t * ghost,
                                   p4est_mesh_t * mesh,
                                   p4est_virtual_t * virtual_quads,
                                   p4est_virtual_ghost_t * virtual_ghost,
                                   size_t data_size, void *mirror_data,
                                   void *ghost_data)
{
  p4est_virtual_ghost_exchange_data_end
    (p4est_virtual_ghost_exchange_data_begin
     (p4est, ghost, mesh, virtual_quads, virtual_ghost, data_size,
      mirror_data, ghost_data));
}

p4est_virtual_ghost_exchange_t *
p4est_virtual_ghost_exchange_data_begin (p4est_t * p4est,
                                         p4est_ghost_t * ghost,
                                         p4est_mesh_t * mesh,
                                         p4est_virtual_t * virtual_quads,
                                         p4est_virtual_ghost_t *
                                         virtual_ghost, size_t data_size,
                                         void *mirror_data, void *ghost_data)
{
  const int           num_procs = p4est->mpisize;
  int                 mpiret;
  int                 g, q;
  char               *mem, **sbuf;
  p4est_locidx_t      ng_excl, ng_incl, ng, virt_offset;
  p4est_locidx_t      mirr;
  p4est_locidx_t     *mirr_offset;
  sc_MPI_Request     *r;
  p4est_virtual_ghost_exchange_t *exc;

  /* initialize transient storage */
  exc = P4EST_ALLOC_ZERO (p4est_virtual_ghost_exchange_t, 1);
  exc->p4est = p4est;
  exc->ghost = ghost;
  exc->virtual_quads = virtual_quads;
  exc->virtual_ghost = virtual_ghost;
  exc->level = -1;
  exc->data_size = data_size;
  exc->ghost_data = ghost_data;
  sc_array_init (&exc->requests, sizeof (sc_MPI_Request));
  sc_array_init (&exc->sbuffers, sizeof (char *));

  /* return early if there is nothing to do */
  if (data_size == 0) {
    return exc;
  }

  /* receive data from other processors */
  ng_excl = 0;
  virt_offset = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = ghost->proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      /* check if we expect virtual quadrants and increase size of receive-
         buffer accordingly. */
      for (g = ng_excl; g < ng_incl; ++g) {
        if (virtual_quads->virtual_gflags[g] != -1) {
          ng += P4EST_CHILDREN;
        }
      }
      r = (sc_MPI_Request *) sc_array_push (&exc->requests);
      mpiret =
        sc_MPI_Irecv ((char *) ghost_data +
                      (virt_offset + ng_excl) * data_size, ng * data_size,
                      sc_MPI_BYTE, q, P4EST_COMM_GHOST_EXCHANGE,
                      p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      /* add number of virtual quadrants to offset */
      virt_offset += ng - ng_incl + ng_excl;
      ng_excl = ng_incl;
    }
  }
  P4EST_ASSERT (ng_excl == (p4est_locidx_t) ghost->ghosts.elem_count);
  P4EST_ASSERT (virtual_quads->ghost_num_virtuals == virt_offset);

  /* send data to other processors */
  /* first determine where to find data of each mirror */
  mirr_offset = P4EST_ALLOC (p4est_locidx_t, ghost->mirrors.elem_count);
  memset (mirr_offset, (char) -1,
          ghost->mirrors.elem_count * sizeof (p4est_locidx_t));
  virt_offset = 0;
  mirr = 0;
  for (q = 0; q < p4est->local_num_quadrants; ++q) {
    if (q == mesh->mirror_qid[mirr]) {
      mirr_offset[mirr] = q + virt_offset;
      ++mirr;
    }
    if (-1 != virtual_quads->virtual_qflags[q]) {
      virt_offset += P4EST_CHILDREN;
    }
  }
  ng_excl = 0;
  for (q = 0; q < num_procs; ++q) {
    ng_incl = ghost->mirror_proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    if (ng > 0) {
      for (g = ng_excl; g < ng_incl; ++g) {
        if (virtual_ghost->mirror_proc_virtuals[g]) {
          ng += P4EST_CHILDREN;
        }
      }
      /* every peer populates its own send buffer */
      sbuf = (char **) sc_array_push (&exc->sbuffers);
      mem = *sbuf = P4EST_ALLOC (char, ng * data_size);
      for (g = ng_excl; g < ng_incl; ++g) {
        if (virtual_ghost->mirror_proc_virtuals[g]) {
          mirr = ghost->mirror_proc_mirrors[g];
          memcpy (mem, (char *) mirror_data + mirr_offset[mirr] * data_size,
                  (1 + P4EST_CHILDREN) * data_size);
          mem += ((1 + P4EST_CHILDREN) * data_size);
        }
        else {
          mirr = ghost->mirror_proc_mirrors[g];
          memcpy (mem, (char *) mirror_data + mirr_offset[mirr] * data_size,
                  data_size);
          mem += data_size;
        }
      }
      r = (sc_MPI_Request *) sc_array_push (&exc->requests);
      mpiret = sc_MPI_Isend (*sbuf, ng * data_size, sc_MPI_BYTE, q,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);

      virt_offset = ng - ng_incl + ng_excl;
      ng_excl = ng_incl;
    }
  }
  P4EST_FREE (mirr_offset);

  /* we are done posting the messages */
  return exc;
}

void
p4est_virtual_ghost_exchange_data_end (p4est_virtual_ghost_exchange_t * exc)
{
  int                 mpiret;
  size_t              zz;
  char              **sbuf;

  /* don't confuse it with p4est_virtual_ghost_exchange_levels_end */
  P4EST_ASSERT (!exc->is_levels);

  /* wait for messages to complete and clean up */
  mpiret = sc_MPI_Waitall (exc->requests.elem_count, (sc_MPI_Request *)
                           exc->requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_reset (&exc->requests);
  for (zz = 0; zz < exc->sbuffers.elem_count; ++zz) {
    sbuf = (char **) sc_array_index (&exc->sbuffers, zz);
    P4EST_FREE (*sbuf);
  }
  sc_array_reset (&exc->sbuffers);
  P4EST_FREE (exc);

  return;
}

void
p4est_virtual_ghost_exchange_data_level (p4est_t * p4est,
                                         p4est_ghost_t * ghost,
                                         p4est_mesh_t * mesh,
                                         p4est_virtual_t * virtual_quads,
                                         p4est_virtual_ghost_t *
                                         virtual_ghost, int level,
                                         size_t data_size, void **mirror_data,
                                         void **ghost_data)
{
  p4est_virtual_ghost_exchange_data_level_end
    (p4est_virtual_ghost_exchange_data_level_begin
     (p4est, ghost, mesh, virtual_quads, virtual_ghost, level, data_size,
      mirror_data, ghost_data));
}

p4est_virtual_ghost_exchange_t *
p4est_virtual_ghost_exchange_data_level_begin (p4est_t * p4est,
                                               p4est_ghost_t * ghost,
                                               p4est_mesh_t * mesh,
                                               p4est_virtual_t *
                                               virtual_quads,
                                               p4est_virtual_ghost_t *
                                               virtual_ghost, int level,
                                               size_t data_size,
                                               void **mirror_data,
                                               void **ghost_data)
{
  const int           num_procs = p4est->mpisize;
  int                 mpiret;
  int                 q;
  char               *mem, **sbuf;
  p4est_locidx_t      ng_excl, ng_incl, ng, theg;
  p4est_locidx_t      offset, first_real, first_virt;
  p4est_locidx_t      lmatches;
  p4est_locidx_t      ghost_idx, mirror_idx, mirror_qid;
  p4est_quadrant_t   *m;
  sc_MPI_Request     *r;
  p4est_virtual_ghost_exchange_t *exc;

  P4EST_ASSERT (mesh->ghost_level != NULL);
  P4EST_ASSERT (mesh->quad_level != NULL);
  P4EST_ASSERT (virtual_quads->quad_qreal_offset != NULL);
  P4EST_ASSERT (virtual_quads->quad_qvirtual_offset != NULL);
  P4EST_ASSERT (virtual_quads->quad_greal_offset != NULL);
  P4EST_ASSERT (virtual_quads->quad_gvirtual_offset != NULL);

  /* initialize transient storage */
  exc = P4EST_ALLOC_ZERO (p4est_virtual_ghost_exchange_t, 1);
  exc->p4est = p4est;
  exc->ghost = ghost;
  exc->virtual_quads = virtual_quads;
  exc->virtual_ghost = virtual_ghost;
  exc->level = level;
  exc->is_levels = 1;
  exc->data_size = data_size;
  exc->ghost_data = ghost_data;
  sc_array_init (&exc->requests, sizeof (sc_MPI_Request));
  sc_array_init (&exc->sbuffers, sizeof (char *));

  P4EST_ASSERT (0 <= level && level <= P4EST_QMAXLEVEL);

  /* return early if there is nothing to do */
  if (data_size == 0) {
    return exc;
  }

  /* receive data from other processors */
  for (ng_excl = 0, q = 0; q < num_procs; ++q) {
    first_real = -1;
    first_virt = -1;
    ng_incl = ghost->proc_offsets[q + 1];
    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);

    /* if there are any ghosts to be received check if they match the level
     * criterion (both real and virtual) */
    if (ng > 0) {
      P4EST_ASSERT (q != p4est->mpirank);
      for (lmatches = 0, theg = 0;
           theg < (mesh->ghost_level + level)->elem_count; ++theg) {
        ghost_idx =
          *(p4est_locidx_t *) sc_array_index (mesh->ghost_level + level,
                                              theg);
        if (ng_excl <= ghost_idx && ghost_idx < ng_incl) {
          if (first_real == -1) {
            first_real = ghost_idx;
          }
          ++lmatches;
        }
      }
      for (theg = 0;
           theg < (virtual_quads->virtual_glevels + level)->elem_count;
           ++theg) {
        ghost_idx =
          *(p4est_locidx_t *) sc_array_index (virtual_quads->virtual_glevels +
                                              level, theg);
        if (ng_excl <= ghost_idx && ghost_idx < ng_incl) {
          if (first_virt == -1) {
            first_virt = ghost_idx;
          }
          lmatches += P4EST_CHILDREN;
        }
      }

      if (lmatches > 0) {
        P4EST_ASSERT (first_virt != -1 || first_real != -1);
        P4EST_ASSERT (first_real != first_virt);

        r = (sc_MPI_Request *) sc_array_push (&exc->requests);
        /* We can use the ghost data memory as is */
        /* 4 cases to distinguish:
         * a) no virtual quads
         * b) no real quads
         * c) real and virtual quads, first quad is real
         * d) real and virtual quads, first quad is virtual
         * leads to 2 different cases: offset is determined from real quadrant
         * (a, c) or from virtual quadrant (b, d)
         */
        if ((first_virt == -1)
            || ((first_real < first_virt) && (first_real != -1))) {
          offset = virtual_quads->quad_greal_offset[first_real];
        }
        else if ((first_real == -1)
                 || ((first_virt < first_real) && (first_virt != -1))) {
          offset = virtual_quads->quad_gvirtual_offset[first_virt];
        }
        else {
          SC_ABORT_NOT_REACHED ();
        }
        mpiret =
          sc_MPI_Irecv ((char *) ghost_data[level] + offset * data_size,
                        lmatches * data_size, sc_MPI_BYTE, q, (200 + level),
                        p4est->mpicomm, r);
        SC_CHECK_MPI (mpiret);
      }
    }
    ng_excl = ng_incl;
  }

  /* send data to other processors */
  for (ng_excl = 0, q = 0; q < num_procs; ++q) {
    ng_incl = ghost->mirror_proc_offsets[q + 1];

    ng = ng_incl - ng_excl;
    P4EST_ASSERT (ng >= 0);
    /* if there is anything to send to the neighboring process for all levels */
    if (ng > 0) {
      /* count how many quadrants match the level criterion:
       * either they have a matching level or are one level coarser holding
       * virtual quadrants (which the neighboring process has to expect)
       */
      for (lmatches = 0, theg = ng_excl; theg < ng_incl; ++theg) {
        mirror_idx = ghost->mirror_proc_mirrors[theg];
        mirror_qid = mesh->mirror_qid[mirror_idx];
        m = p4est_quadrant_array_index (&ghost->mirrors, mirror_idx);
        if (m->level == level) {
          ++lmatches;
        }
        else if (m->level == (level - 1)
                 && virtual_ghost->mirror_proc_virtuals[theg]) {
          P4EST_ASSERT (-1 != virtual_quads->virtual_qflags[mirror_qid]);
          lmatches += P4EST_CHILDREN;
        }
      }
      if (0 < lmatches) {
        /* every peer populates its own send buffer */
        sbuf = (char **) sc_array_push (&exc->sbuffers);
        mem = *sbuf = P4EST_ALLOC (char, lmatches * data_size);

        for (theg = ng_excl; theg < ng_incl; ++theg) {
          mirror_idx = ghost->mirror_proc_mirrors[theg];
          m = p4est_quadrant_array_index (&ghost->mirrors, mirror_idx);
          mirror_qid = mesh->mirror_qid[mirror_idx];
          if (m->level == level) {
            offset = virtual_quads->quad_qreal_offset[mirror_qid];
            P4EST_ASSERT (0 <= offset && offset <
                          ((mesh->quad_level + level)->elem_count +
                           P4EST_CHILDREN * (virtual_quads->virtual_qlevels +
                                             level)->elem_count));
            memcpy (mem, (char *) mirror_data[level] + data_size * offset,
                    data_size);
            mem += data_size;
          }
          else if (m->level == (level - 1)
                   && virtual_ghost->mirror_proc_virtuals[theg]) {
            offset = virtual_quads->quad_qvirtual_offset[mirror_qid];
            P4EST_ASSERT (0 <= offset && offset <
                          ((mesh->quad_level + level)->elem_count +
                           P4EST_CHILDREN * (virtual_quads->virtual_qlevels +
                                             level)->elem_count));
            memcpy (mem, (char *) mirror_data[level] + data_size * offset,
                    P4EST_CHILDREN * data_size);
            mem += P4EST_CHILDREN * data_size;
          }
        }
        r = (sc_MPI_Request *) sc_array_push (&exc->requests);
        mpiret = sc_MPI_Isend (*sbuf, lmatches * data_size, sc_MPI_BYTE, q,
                               (200 + level), p4est->mpicomm, r);
        SC_CHECK_MPI (mpiret);
      }
      ng_excl = ng_incl;
    }
  }

  /* we are done posting messages */
  return exc;
}

void
p4est_virtual_ghost_exchange_data_level_end (p4est_virtual_ghost_exchange_t *
                                             exc)
{
  int                 mpiret;
  size_t              zz;
  char              **sbuf;

  /* don't confuse it with p4est_virtual_ghost_exchange_end */
  P4EST_ASSERT (exc->is_levels);

  /* wait for messages to complete and clean up */
  mpiret = sc_MPI_Waitall (exc->requests.elem_count, (sc_MPI_Request *)
                           exc->requests.array, sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_reset (&exc->requests);
  for (zz = 0; zz < exc->sbuffers.elem_count; ++zz) {
    sbuf = (char **) sc_array_index (&exc->sbuffers, zz);
    P4EST_FREE (*sbuf);
  }
  sc_array_reset (&exc->sbuffers);

  P4EST_FREE (exc);

  return;
}

/* -------------------------------------------------------------------------- */
/* |                            Neighbor search                             | */
/* -------------------------------------------------------------------------- */

#ifndef P4_TO_P8
/* *INDENT-OFF* */
const int           p4est_virtual_face_neighbors_search_opts[P4EST_CHILDREN]
                                                            [P4EST_FACES] =
{{  4,  1,  6,  2 },
 {  0,  5, 10,  3 },
 {  8,  3,  0,  7 },
 {  2,  9,  1, 11 }};

const int           p4est_virtual_corner_neighbors_search_opts[P4EST_CHILDREN]
                                                              [P4EST_CHILDREN] =
{{ 12, 10,  8,  3 },
 {  6, 13,  2,  9 },
 {  4,  1, 14, 11 },
 {  0,  5,  7, 15 }};
/* *INDENT-ON* */
#endif /* P4_TO_P8 */

/* -------------------------------------------------------------------------- */
/* |                       Neighbor search: Bit Magic                       | */
/* -------------------------------------------------------------------------- */
typedef int         (*is_max_t) (const int);

/** Determine if corner index is situated at x = 0 or x = 1 in local
 * coordinates.
 */
static int
is_max_x (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return c & 1;
}

/** Determine if corner index is situated at y = 0 or y = 1 in local
 * coordinates.
 */
static int
is_max_y (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return (is_max_x (c >> 1));
}

#ifdef P4_TO_P8
/** Determine if corner index is situated at z = 0 or z = 1 in local
 * coordinates.
 */
static int
is_max_z (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return (is_max_x (c >> 2));
}
#endif /* P4_TO_P8 */

static is_max_t     is_max[P4EST_DIM] = {
  is_max_x,
  is_max_y
#ifdef P4_TO_P8
    , is_max_z
#endif /* P4_TO_P8 */
};

#ifdef P4_TO_P8
typedef int         (*mask_dir_t) (const int);

/** Obtain index 0 .. 3 indicating the position of corner wrt. of its y and z
 * components in Morton-order where y is incremented first.
 */
static int
mask_x (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return (c >> 1);
}

/** Obtain index 0 .. 3 indicating the position of corner wrt. of its x and z
 * components in Morton-order where x is incremented first.
 */
static int
mask_y (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return (((c >> 1) & 2) + (c & 1));
}

/** Obtain index 0 .. 3 indicating the position of corner wrt. of its x and y
 * components in Morton-order where x is incremented first.
 */
static int
mask_z (const int c)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  return c & 3;
}

static mask_dir_t   masks[P4EST_DIM] = { mask_x, mask_y, mask_z };
#endif /* P4_TO_P8 */

/** Get face adjacent to the corner in a specific spatial direction.
 * \param     [in]  c        Corner index 0 .. 2^(dim - 1).
 * \param     [in]  d        Spatial direction:
 *                              d = 0 .. x-axis
 *                              d = 1 .. y-axis
 *                              d = 2 .. z-axis (3D only)
 * \return    face index     0 .. (2 * dim)
 */
static int
get_adjacent_face (const int c, const int d)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= d && d < P4EST_DIM);
  return 2 * d + is_max[d] (c);
}

#ifdef P4_TO_P8
/** Get face opposite to the corner in a specific spatial direction.
 * CAUTION: Works for 2D as well but is unused
 * \param     [in]  c        Corner index 0 .. 2^(dim - 1).
 * \param     [in]  d        Spatial direction:
 *                              d = 0 .. x-axis
 *                              d = 1 .. y-axis
 *                              d = 2 .. z-axis (3D only)
 * \return    face index     0 .. (2 * dim)
 */
static int
get_opposite_face (const int c, const int d)
{
  int                 oc = c ^ (P4EST_CHILDREN - 1);
  return get_adjacent_face (oc, d);
}

/** Get edge index adjacent to corner parallel to a specific spatial direction.
 * \param     [in]  c        Corner index 0 .. 2^(dim - 1).
 * \param     [in]  d        Spatial direction:
 *                              d = 0 .. x-axis
 *                              d = 1 .. y-axis
 *                              d = 2 .. z-axis
 * \return    edge index     0 .. 12
 */
static int
get_adjacent_edge (const int c, const int d)
{
  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= d && d < P4EST_DIM);
  return 4 * d + masks[d] (c);
}

#endif /* P4_TO_P8 */

/** insert one element into result arrays
 * \param     [out] n_encs   Result array of int containing the neighbor's
 *                           encoding.  If NULL n_enc is ignored.
 * \param     [out] n_qids   Result array of int containing the neighbor's
 *                           quadrant index.  If NULL n_qid is ignored.
 * \param     [out] n_vids   Result array of int containing the neighbor's
 *                           virtual index.  If NULL n_vid is ignored.
 * \param [in]      n_enc    Element that will be appended to n_encs.
 * \param [in]      n_qid    Element that will be appended to qids.
 * \param [in]      n_vid    Element that will be appended to vids.
 */
static int
insert_neighbor_elements (sc_array_t * n_encs, sc_array_t * n_qids,
                          sc_array_t * n_vids, int8_t n_enc,
                          p4est_locidx_t n_qid, int n_vid)
{
  int                *insert_int;

  if (n_encs) {
    insert_int = (int *) sc_array_push (n_encs);
    *insert_int = n_enc;
  }
  if (n_qids) {
    insert_int = (int *) sc_array_push (n_qids);
    *insert_int = n_qid;
  }
  if (n_vids) {
    insert_int = (int *) sc_array_push (n_vids);
    *insert_int = n_vid;
  }

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
  return -1;
}

#ifdef P4_TO_P8
/** Set constants to get opposite edges and corners across a face, originating
 * from an edge.  Note, that we are not interested in diagonally opposite
 * elements.
 * \param [in]      edge_dir
 * \param [in]      face_dir
 * \param     [out] edge_offset
 * \param     [out] corner_offset
 */
static int
set_xor_constants_edge (int edge_dir, int face_dir, int *edge_offset,
                        int *corner_offset)
{
  /** edges parallel to x axis */
  if (0 <= edge_dir && edge_dir < 4) {
    /** faces normal to y axis */
    if (2 <= face_dir && face_dir < 4) {
      *edge_offset = 2;
      *corner_offset = 4;
    }
    /** faces normal to z axis */
    else if (4 <= face_dir && face_dir < 6) {
      *edge_offset = 1;
      *corner_offset = 2;
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }
  /** edges parallel to y axis */
  else if (4 <= edge_dir && edge_dir < 8) {
    /** faces normal to x axis */
    if (0 <= face_dir && face_dir < 2) {
      *edge_offset = 2;
      *corner_offset = 4;
    }
    /** faces normal to z axis */
    else if (4 <= face_dir && face_dir < 6) {
      *edge_offset = 1;
      *corner_offset = 1;
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }
  /** edges parallel to z axis */
  else if (8 <= edge_dir && edge_dir < 12) {
    /** faces normal to x axis */
    if (0 <= face_dir && face_dir < 2) {
      *edge_offset = 2;
      *corner_offset = 2;
    }
    /** faces normal to y axis */
    else if (2 <= face_dir && face_dir < 4) {
      *edge_offset = 1;
      *corner_offset = 1;
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  return 0;
}
#endif /* P4_TO_P8 */

/** Go in the direction of a diagonally opposite corner within a face.
 * \param [in]      face_dir    Face index 0 .. P4EST_FACES
 * \param     [out] offset      XOR operand for obtaining respective corner
 *                              index in diagonally opposite direction within
 *                              the current face within the current quadrant.
 */
static int
set_xor_constants_corner_inside (int face_dir, int *offset)
{
  int                 i;
  *offset = 0;

  for (i = 0; i < P4EST_DIM; ++i) {
    if (face_dir / 2 != i) {
      *offset += (1 << i);
    }
  }

  return 0;
}

#ifdef P4_TO_P8
/** Go ortogonal to the face.  Works for 2D as well but is unused.
 * \param [in]      face_dir    Face index 0 .. P4EST_FACES
 * \param     [out] offset      XOR operand for obtaining respective corner
 *                              index in direction normal to the face within
 *                              quadrant.
 */
static int
set_xor_constants_corner_orth (int face_dir, int *offset)
{
  int                 i;
  *offset = 0;

  for (i = 0; i < P4EST_DIM; ++i) {
    if (face_dir / 2 == i) {
      *offset += (1 << i);
    }
  }

  return 0;
}

/** Transform a hanging quadrant subindex to a corresponding index inside of a
 * virtual quadrant.
 * This is a very similar transformation compared to
 * p4est_connectivity_face_neighbor_face_corner_orientation, however the
 * permutation sets for a permutation ref of 1 or 2 is swapped wrt. this
 * function.
 *
 * \param[in] fc     Subindex 0 .. 3 of the virtual quadrant
 * \param[in] f      Face index 0 .. 5 in which we search for face neighbors
 * \param[in] nf     Face index 0 .. 5 which is adjacent to f in the
 *                   neighboring quadrant
 * \param[in] o      Orientation 0 .. 3 between both faces
 * \returns          Equivalent virtual index 0 .. 7 to current subindex
 */
int
edge_transform_subindex (int fc, int f, int nf, int o)
{
  int                 pref, pset, nfc;

  P4EST_ASSERT (0 <= fc && fc < P4EST_HALF);
  P4EST_ASSERT (0 <= f && f < P4EST_FACES);
  P4EST_ASSERT (0 <= nf && nf < P4EST_FACES);
  P4EST_ASSERT (0 <= o && o < P4EST_HALF);

  /* *INDENT-OFF* */
  const int hanging_corner_permutations[8][4] =
    {{ 0, 1, 2, 3 },
     { 0, 2, 1, 3 },
     { 1, 0, 3, 2 },
     { 1, 3, 0, 2 },
     { 2, 0, 3, 1 },
     { 2, 3, 0, 1 },
     { 3, 1, 2, 0 },
     { 3, 2, 1, 0 }};
  const int hanging_corner_permutation_refs[3][4] =
    {{ 1, 2, 5, 6 },
     { 0, 4, 3, 7 },
     { 0, 3, 4, 7 }};
  /* *INDENT-ON* */

  pref = p8est_face_permutation_refs[f][nf];
  pset = hanging_corner_permutation_refs[pref][o];
  nfc = hanging_corner_permutations[pset][fc];

  P4EST_ASSERT (0 <= nfc && nfc < P4EST_HALF);

  return p8est_face_corners[f][nfc];
}
#endif /* P4_TO_P8 */

/** Search the correct virtual quadrant within a double-sized neighbor for real
 * quadrants.
 * \param[in]      direction
 * \param[in][out] encoding
 * \param    [out] vid
 */
static int
get_real_neighbor_vid (int dir, int *encoding, int *vid)
{
  int                 n_orientation, n_subquad, n_entity;
  int                 l_same_size, u_same_size;
  int                 l_double_size, u_double_size;
  int                 l_half_size, u_half_size;

  if (0 <= dir && dir < P4EST_FACES) {
    /* get orientation from encoding */
#ifdef P4_TO_P8
    l_same_size = 0;
    u_same_size = 4 * P4EST_FACES;
    l_double_size = u_same_size;
    u_double_size = 120;
    l_half_size = -24;
    u_half_size = l_same_size;
#else /* P4_TO_P8 */
    l_same_size = 0;
    u_same_size = 2 * P4EST_FACES;
    l_double_size = u_same_size;
    u_double_size = 24;
    l_half_size = -8;
    u_half_size = l_same_size;
#endif /* P4_TO_P8 */
    p4est_mesh_decode_encoding (*encoding, P4EST_FACES, l_same_size,
                                u_same_size, l_double_size, u_double_size,
                                l_half_size, u_half_size, &n_subquad,
                                &n_orientation, &n_entity);

    *encoding = P4EST_FACES * n_orientation + n_entity;
    *vid = p4est_face_corners[n_entity][n_subquad];
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= dir && dir < P4EST_FACES + P8EST_EDGES) {
    l_same_size = 0;
    u_same_size = 2 * P8EST_EDGES;
    l_double_size = u_same_size;
    u_double_size = 72;
    l_half_size = -24;
    u_half_size = l_same_size;

    p4est_mesh_decode_encoding (*encoding, P8EST_EDGES, l_same_size,
                                u_same_size, l_double_size, u_double_size,
                                l_half_size, u_half_size, &n_subquad,
                                &n_orientation, &n_entity);
    *encoding = P8EST_EDGES * n_orientation + n_entity;
    *vid =
      p8est_connectivity_edge_neighbor_edge_corner (n_subquad, n_orientation);
    *vid = p8est_edge_corners[n_entity][*vid];
  }
  else if (P4EST_FACES + P8EST_EDGES <= dir
           && dir < P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= dir && dir < P4EST_FACES + P4EST_CHILDREN)
#endif /* P4_TO_P8 */
  {
    /* for corners there is nothing to be done. vid = encoding */
    *vid = *encoding;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  return 0;
}

/** Get the corner neighbor of a real quadrant which is hanging across a face.
 * \param [in]      p4est          The forest.
 * \param [in]      ghost          Ghost layer.
 * \param [in]      mesh           Neighbor information.
 * \param [in]      virtual_quads  Virtual quadrants.
 * \param [in]      qid            Current quadrant id.
 * \param [in]      dir            Face index across which to search neighbors.
 * \param     [out] n_encs         Array of encodings to be populated.
 * \param     [out] n_qids         Array of qids to be populated.
 * \param     [out] n_vids         Array of vids to be populated.
 */
static int
get_corner_hanging_face (p4est_t * p4est, p4est_ghost_t * ghost,
                         p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                         p4est_locidx_t qid, int dir, sc_array_t * n_encs,
                         sc_array_t * n_qids, sc_array_t * n_vids)
{
#ifdef P4EST_ENABLE_DEBUG
  int                 nqid;
  p4est_quadrant_t   *q, *n;
#endif /* P4EST_ENABLE_DEBUG */

  int                 enc;
  int                 offset;
  int                *int_ins;
  int                 n_orientation, n_subquad, n_entity;
  int                 l_same_size, u_same_size;
  int                 l_double_size, u_double_size;
  int                 l_half_size, u_half_size;

  P4EST_ASSERT (n_vids->elem_count == 0);
  P4EST_ASSERT (n_encs->elem_count == 0);
  P4EST_ASSERT (n_qids->elem_count == 0);

#ifdef P4_TO_P8
  l_same_size = 0;
  u_same_size = 4 * P4EST_FACES;
  l_double_size = u_same_size;
  u_double_size = 120;
  l_half_size = -24;
  u_half_size = l_same_size;
#else /* P4_TO_P8 */
  l_same_size = 0;
  u_same_size = 2 * P4EST_FACES;
  l_double_size = u_same_size;
  u_double_size = 24;
  l_half_size = -8;
  u_half_size = l_same_size;
#endif /* P4_TO_P8 */

  p4est_mesh_get_neighbors (p4est, ghost, mesh, qid, dir, NULL, n_encs,
                            n_qids);

  P4EST_ASSERT (n_qids->elem_count == 1);
  P4EST_ASSERT (n_encs->elem_count == 1);
  enc = *(int *) sc_array_index (n_encs, 0);

#ifdef P4EST_ENABLE_DEBUG
  nqid = *(int *) sc_array_index (n_qids, 0);
  P4EST_ASSERT (0 <= nqid
                && nqid <
                (mesh->local_num_quadrants + mesh->ghost_num_quadrants));
  q = p4est_mesh_get_quadrant (p4est, mesh, qid);
  n =
    nqid < mesh->local_num_quadrants ? p4est_mesh_get_quadrant (p4est, mesh,
                                                                nqid) :
    p4est_quadrant_array_index (&ghost->ghosts,
                                nqid - mesh->local_num_quadrants);
  P4EST_ASSERT (n->level + 1 == q->level);
#endif /* P4EST_ENABLE_DEBUG */
  p4est_mesh_decode_encoding (enc, P4EST_FACES, l_same_size, u_same_size,
                              l_double_size, u_double_size, l_half_size,
                              u_half_size, &n_subquad, &n_orientation,
                              &n_entity);

  set_xor_constants_corner_inside (n_entity, &offset);
  enc = p4est_face_corners[n_entity][n_subquad];
  int_ins = sc_array_push (n_vids);
  *int_ins = enc ^ offset;
  int_ins = sc_array_index (n_encs, 0);
  *int_ins = enc;

  return 0;
}

#ifdef P4_TO_P8
/** Get the corner neighbor of a real quadrant which is hanging across a edge.
 * \param [in]      p4est          The forest.
 * \param [in]      ghost          Ghost layer.
 * \param [in]      mesh           Neighbor information.
 * \param [in]      virtual_quads  Virtual quadrants.
 * \param [in]      qid            Current quadrant id.
 * \param [in]      c_dir          Original direction of corner-neighbor query.
 * \param [in]      e_dir          Edge index across which to search neighbors.
 * \param [in]      sid            Current quadrant's sibling id 0 .. 7 in its
 *                                 ancestor quadrant.
 * \param     [out] n_encs         Array of encodings to be populated.
 * \param     [out] n_qids         Array of qids to be populated.
 * \param     [out] n_vids         Array of vids to be populated.
 */
static int
get_corner_hanging_edge (p4est_t * p4est, p4est_ghost_t * ghost,
                         p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                         p4est_locidx_t qid, int c_dir, int e_dir, int sid,
                         sc_array_t * n_encs, sc_array_t * n_qids,
                         sc_array_t * n_vids)
{
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  int                 search_dir;
  int                 tmp_subquad, tmp_ori, tmp_entity;
  p4est_quadrant_t   *n, c_n;
  p4est_topidx_t      n_tree_id, n_rank;
  p4est_tree_t       *n_tree;
  int                 n_sibling, n_enc, n_qid, n_vid;
  int                 l_same_size_edge, u_same_size_edge, l_double_size_edge;
  int                 u_double_size_edge, l_half_size_edge, u_half_size_edge;
  p4est_qcoord_t      inc;

  P4EST_ASSERT (0 <= e_dir && e_dir < P8EST_EDGES);
  l_half_size_edge = -24;
  u_half_size_edge = l_same_size_edge = 0;
  u_same_size_edge = l_double_size_edge = 24;
  u_double_size_edge = 72;

  p4est_mesh_get_neighbors (p4est, ghost, mesh, qid, e_dir + P4EST_FACES,
                            NULL, n_encs, n_qids);
  /* domain boundary: no edge neighbor present */
  if (0 == n_encs->elem_count) {
    return 0;
  }
  /* FIXME: This assertion is only valid in brick-like scenarios.  However, as
   * we are already part of a hanging corner, there can only be same or double
   * size neighbors.  In case of half-size neighbors the forest would not be
   * 2:1-balanced.*/
  P4EST_ASSERT (1 == n_qids->elem_count);

  /* check result of neighbor query */
  n_qid = *(int *) sc_array_index (n_qids, 0);
  P4EST_ASSERT (0 <= n_qid && n_qid < (lq + gq));

  n_enc = *(int *) sc_array_index (n_encs, 0);
  P4EST_ASSERT (l_same_size_edge <= n_enc && n_enc < u_double_size_edge);
  p4est_mesh_decode_encoding (n_enc, P8EST_EDGES, l_same_size_edge,
                              u_same_size_edge, l_double_size_edge,
                              u_double_size_edge, l_half_size_edge,
                              u_half_size_edge, &tmp_subquad, &tmp_ori,
                              &tmp_entity);
  sc_array_truncate (n_encs);

  /* same size neighbor */
  if (l_same_size_edge <= n_enc && n_enc < u_same_size_edge) {
    sc_array_truncate (n_qids);
    /* Find correct quadrant first.  The direction in which to search is found
     * by searching for the opposite face in the direction of the edge.
     */
    n = n_qid < lq ? p4est_mesh_get_quadrant (p4est, mesh, n_qid) :
      p4est_quadrant_array_index (&ghost->ghosts, (n_qid - lq));
    n_sibling = p4est_quadrant_child_id (n);
    search_dir = get_opposite_face (n_sibling, tmp_entity / 4);

    /* If neighboring quadrant is local, search the respective neighbor via
     * p4est_mesh in the respective direction.
     * Else, we have to construct the respective quadrant, determine its owner
     * rank as well as its tree id.  If the constructed quadrant is local we
     * perform a binary search among the local quadrants of that tree.  If the
     * owner rank is a different rank we search among the ghosts owned by that
     * rank. */
    if (0 <= n_qid && n_qid < lq) {
      p4est_mesh_get_neighbors (p4est, ghost, mesh, n_qid, search_dir, NULL,
                                NULL, n_qids);
      P4EST_ASSERT (n_qids->elem_count == 1);
    }
    else {
      /* construct quadrant */
      memcpy (&c_n, n, sizeof (p4est_quadrant_t));
      inc = P4EST_QUADRANT_LEN (n->level);
      /* now we only want to distinguish x, y, and z direction */
      search_dir *= 0.5;
      switch (search_dir) {
      case 0:
        if (is_max[search_dir] (n_sibling)) {
          c_n.x = c_n.x & (~inc);
        }
        else {
          c_n.x = c_n.x | inc;
        }
        break;
      case 1:
        if (is_max[search_dir] (n_sibling)) {
          c_n.y = c_n.y & (~inc);
        }
        else {
          c_n.y = c_n.y | inc;
        }
        break;
      case 2:
        if (is_max[search_dir] (n_sibling)) {
          c_n.z = c_n.z & (~inc);
        }
        else {
          c_n.z = c_n.z | inc;
        }
        break;
      default:
        SC_ABORT_NOT_REACHED ();
      }
      /* determine its tree id */
      n_tree_id = mesh->ghost_to_tree[n_qid - lq];
      n_rank = p4est_quadrant_find_owner (p4est, n_tree_id, -1, &c_n);
      P4EST_ASSERT (0 <= n_rank && n_rank < p4est->mpisize);
      if (n_rank == p4est->mpirank) {
        n_tree = p4est_tree_array_index (p4est->trees, n_tree_id);
        n_qid =
          sc_array_bsearch (&n_tree->quadrants, &c_n,
                            p4est_quadrant_compare) +
          n_tree->quadrants_offset;
        P4EST_ASSERT (0 <= n_qid && n_qid < lq);
      }
      else {
        n_qid = p4est_ghost_bsearch (ghost, n_rank, n_tree_id, &c_n);
        P4EST_ASSERT (0 <= n_qid && n_qid < gq);
        n_qid += lq;
      }
    }
    n_vid = -1;
    n_enc = p8est_connectivity_edge_neighbor_corner (sid, e_dir, tmp_entity,
                                                     tmp_ori);
    P4EST_ASSERT (0 <= n_enc && n_enc < P4EST_CHILDREN);
    if (n_qids->elem_count) {
      insert_neighbor_elements (n_encs, NULL, n_vids, n_enc, -1, n_vid);
    }
    else {
      insert_neighbor_elements (n_encs, n_qids, n_vids, n_enc, n_qid, n_vid);
    }
  }
  /* double size neighbor:
   * We have obtained the correct quadrant index.  This leaves calculating the
   * respective same-size encoding and determine the quadrant's vid.
   */
  else if (l_double_size_edge <= n_enc && n_enc < u_double_size_edge) {
    /* get neighbor's encoding and vid:
     * vid is obtained by transforming original search direction across edge
     * into neighboring quadrant.  The querying quadrant is seen from the
     * opposite corner on that edge which can be obtained by XOR with 2 to the
     * power of the index of the spatial direction the edge is parallel to. */
    P4EST_ASSERT (((0 <= n_qid && n_qid < lq)
                   && -1 != virtual_quads->virtual_qflags[n_qid])
                  || (lq <= n_qid && n_qid < (lq + gq)
                      && -1 != virtual_quads->virtual_gflags[n_qid - lq]));
    n_vid = p8est_connectivity_edge_neighbor_corner (c_dir, e_dir, tmp_entity,
                                                     tmp_ori);
    P4EST_ASSERT (0 <= n_vid && n_vid < P4EST_CHILDREN);
    n_enc = n_vid ^ (1 << (tmp_entity / 4));
    P4EST_ASSERT (0 <= n_enc && n_enc < P4EST_CHILDREN);
    insert_neighbor_elements (n_encs, NULL, n_vids, n_enc, -1, n_vid);
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }

  return 0;
}
#endif /* P4_TO_P8 */

/** Neighbor lookup for real quadrants that will only return same sized
 * quadrants as neighbors, real or virtual. That means compared to
 * \ref p4est_mesh_get_neighbors we obtain the same result for same-sized
 * quadrants.  For double-sized neighbors we find the respective virtual
 * quadrant while half-sized neighbors are discarded.
 * \param[in]      p4est          The forest.
 * \param[in]      ghost          Ghost layer
 * \param[in]      mesh           Neighbor information of real quadrants.
 * \param[in]      virtual_quads  Information which quadrants host virtual
 *                                quadrants.
 * \param[in]      qid            Local quadrant-id of the quadrant whose
 *                                neighbors we are searching.
 * \param[in]      dir            Direction we search for neighbors.
 *                                2D:
 *                                  0 ..  3 neighbor(-s) across face i,
 *                                  4 ..  7 neighbor(-s) across corner i-4.
 *                                3D:
 *                                  0 ..  5 neighbor(-s) across face i,
 *                                  6 .. 17 neighbor(-s) across edge i-6
 *                                 18 .. 25 neighbor(-s) across corner i-18.
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
static int
get_neighbor_real (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                   p4est_locidx_t qid, int dir, sc_array_t * n_encs,
                   sc_array_t * n_qids, sc_array_t * n_vids)
{
  int                 i;
  int                 offset = P4EST_FACES;
  int                 n_vid, n_enc;
  int                *int_ins;
  p4est_quadrant_t   *neighbor, *current;
  sc_array_t         *n_quads = sc_array_new (sizeof (p4est_quadrant_t *));
  int                 tmp_dir;
  int                 sid;
  int                 corner_type, corner_type_norm;
#ifdef P4_TO_P8
  int                 query_dir_face;
  int                 tmp_qid, tmp_enc, tmp_subquad, tmp_entity, tmp_ori;
  int                 tmp_subindex;
  int                 same = -1;
  int                 larger = -1;
  int                 edge_offset, corner_offset;
  int                 r_edge, r_corner, c_corner;
#endif /* P4_TO_P8 */
#ifdef P4_TO_P8
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  int                 l_same_size_face, u_same_size_face;
  int                 l_double_size_face, u_double_size_face;
  int                 l_half_size_face, u_half_size_face;
  l_half_size_face = -24;
  u_half_size_face = l_same_size_face = 0;
  u_same_size_face = l_double_size_face = 24;
  u_double_size_face = 120;
#endif /* P4_TO_P8 */
  current = p4est_mesh_get_quadrant (p4est, mesh, qid);
  /* perform regular neighbor search */
  p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                            dir, n_quads, n_encs, n_qids);
  /* inspect results and populate n_vids array */
  for (i = 0; i < n_encs->elem_count; ++i) {
    /* FIXME:
     * for now only brick-like scenarios are supported. To extend this to an
     * arbitrary connectivity, respective elements need to be invalidated or
     * deleted in the half size case instead of discarding all neighbors. */
    neighbor = *(p4est_quadrant_t **) sc_array_index (n_quads, i);
    n_enc = *(int *) sc_array_index (n_encs, i);
    /* same size */
    if (neighbor->level == current->level) {
      n_vid = -1;
      int_ins = (int *) sc_array_push (n_vids);
      *int_ins = n_vid;
    }
    /* double size:
     * substitute encoding and set correct virtual quadrant */
    else if (neighbor->level == current->level - 1) {
      get_real_neighbor_vid (dir, &n_enc, &n_vid);
      int_ins = (int *) sc_array_index (n_encs, i);
      *int_ins = n_enc;
      int_ins = (int *) sc_array_push (n_vids);
      *int_ins = n_vid;
    }
    /* half size brick and inner */
    else if ((neighbor->level == current->level + 1)
             && (((0 <= dir && dir < P4EST_FACES)
                  && (n_encs->elem_count == P4EST_HALF)) ||
#ifdef P4_TO_P8
                 ((P4EST_FACES <= dir && dir < P4EST_FACES + P8EST_EDGES)
                  && (n_encs->elem_count == 2))
                 ||
                 ((P4EST_FACES + P8EST_EDGES <= dir
                   && dir < P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)
                  && (n_encs->elem_count == 1))
#else /* P4_TO_P8 */
                 ((P4EST_FACES <= dir && dir < P4EST_FACES + P4EST_CHILDREN)
                  && (n_encs->elem_count == 1))
#endif /* P4_TO_P8 */
             )) {
      sc_array_truncate (n_encs);
      sc_array_truncate (n_qids);
      sc_array_destroy (n_quads);
      return 0;
    }
    /* half size non-brick.
     * CAUTION: Not working yet */
    else if (neighbor->level == current->level + 1) {
      /* FIXME: delete P4EST_HALF many quadrants from the respective arrays */
      SC_ABORTF
        ("Current quadrant %i is hanging at a non-brick tree-"
         "boundary. This is not yet implemented", qid);
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }

#ifdef P4_TO_P8
  /* Post process lookup of virtual edge neighbors of hanging edges across a
   * face boundary:
   * edge is touching the face of a large quadrant, ie. there is no matching
   * edge on the other side. By construction the neighboring quadrant must
   * contain virtual quadrants.
   *  Profile:                          Projection A          Projection B
   *   c .. current quad
   *   X .. edge across which
   *        neighbor is searched
   *      -------- ----                 -----------          ---------------
   *    /        /    /|               |     |     |        |          |    |
   *   /        /----/ |               |     |     |        |          |    |
   *  /        /|   /| /               |XXXXX|-----|        |          X----
   *  -------- -+--  |/|               |     |     |        |          |    |
   * |        | X--+-/ |               |  c  |     |        |          |  c |
   * |        |X   |/| /                -----------          ---------------
   * |        |----| |/
   * |        |    | /
   * |        |  c |/
   *  -------- ----
   * I.e. the current quadrant is one of 4 hanging face quadrants and the
   * neighbor currently queried is across one of the 4 edges highlighted by 'X'
   */
  if ((0 == n_encs->elem_count)
      && (P4EST_FACES <= dir && dir < (P4EST_FACES + P8EST_EDGES))
      && (mesh->quad_to_edge[P8EST_EDGES * qid + (dir - P4EST_FACES)] == -1)) {
    /** look in both directions and find we are really in a situation as
     * above. Find the bigger face and obtain the subindex of the current
     * quadrant wrt. this face.
     */

    for (i = 0; i < 2; ++i) {
      tmp_dir = P4EST_FACES * qid + p8est_edge_faces[dir - P4EST_FACES][i];
      if (l_same_size_face <= mesh->quad_to_face[tmp_dir]
          && mesh->quad_to_face[tmp_dir] < u_same_size_face) {
        same = i;
      }
      else if (l_double_size_face <= mesh->quad_to_face[tmp_dir]
               && mesh->quad_to_face[tmp_dir] < u_double_size_face) {
        larger = i;
        p4est_mesh_decode_encoding (mesh->quad_to_face[tmp_dir], P4EST_FACES,
                                    l_same_size_face, u_same_size_face,
                                    l_double_size_face, u_double_size_face,
                                    l_half_size_face, u_half_size_face,
                                    &tmp_subquad, &tmp_ori, &tmp_entity);
      }
      else if (l_half_size_face <= mesh->quad_to_face[tmp_dir]
               && mesh->quad_to_face[tmp_dir] < u_half_size_face) {
        same = larger = -1;
        break;
      }
      else {
        SC_ABORT_NOT_REACHED ();
      }
    }

    if (((0 == same) && (1 == larger))
        || ((1 == same) && (0 == larger))) {
      /* obtain qid of double size neighbor */
      tmp_dir = p8est_edge_faces[dir - P4EST_FACES][same];
      query_dir_face = p8est_edge_faces[dir - P4EST_FACES][larger];
      p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                                query_dir_face, NULL, n_encs, n_qids);

      /* by design we need to find a single neighbor */
      P4EST_ASSERT (1 == n_encs->elem_count);
      tmp_qid = *(int *) sc_array_index (n_qids, 0);
      P4EST_ASSERT (0 <= tmp_qid && tmp_qid < (lq + gq));
      tmp_enc = *(int *) sc_array_index (n_encs, 0);

      sc_array_truncate (n_encs);
      sc_array_truncate (n_qids);

      /* calculate vid:
       * obtain local face corner through transformation,
       * translate to corner index,
       * go in direction of sibling quadrant in local quadrant's face.
       */
      p4est_mesh_decode_encoding (tmp_enc, P4EST_FACES, l_same_size_face,
                                  u_same_size_face, l_double_size_face,
                                  u_double_size_face, l_half_size_face,
                                  u_half_size_face, &tmp_subquad, &tmp_ori,
                                  &tmp_entity);
      tmp_subindex =
        p4est_connectivity_face_neighbor_face_corner (tmp_subquad, tmp_entity,
                                                      query_dir_face,
                                                      tmp_ori);
      P4EST_ASSERT (0 <= tmp_subindex && tmp_subindex < P4EST_HALF);
      tmp_subindex = p4est_face_corners[query_dir_face][tmp_subindex];
      P4EST_ASSERT (0 <= tmp_subindex && tmp_subindex < P4EST_CHILDREN);

      set_xor_constants_corner_orth (tmp_dir, &corner_offset);
      tmp_subindex = tmp_subindex ^ corner_offset;
      P4EST_ASSERT (0 <= tmp_subindex && tmp_subindex < P4EST_CHILDREN);
      n_vid =
        p4est_connectivity_face_neighbor_corner (tmp_subindex,
                                                 query_dir_face, tmp_entity,
                                                 tmp_ori);
      P4EST_ASSERT (0 <= n_vid && n_vid < P4EST_CHILDREN);

      /* Obtain encoding:
       * Two informations are missing:
       * 1) across which edge the current quadrant is seen from the neighbor
       * 2) are the edges flipped or do they have the same orientation
       * The first part can be obtained by transforming our edge into the
       * neighboring quadrant.  The index we are looking for is the opposite
       * edge across the double size neighbor's face.
       * The orientation can be determined by transforming the vid once across
       * the face with known orientation and once across the edge with assumed
       * orientation 0.  We guessed correctly if both indices match.
       */
      r_edge =
        p8est_connectivity_face_neighbor_edge (dir - P4EST_FACES,
                                               query_dir_face, tmp_entity,
                                               tmp_ori);
      set_xor_constants_edge (r_edge, tmp_entity, &edge_offset,
                              &corner_offset);
      r_corner =
        p8est_connectivity_edge_neighbor_corner (n_vid, r_edge,
                                                 dir - P4EST_FACES, 0);
      r_edge = r_edge ^ edge_offset;

      c_corner =
        p4est_connectivity_face_neighbor_corner (n_vid, tmp_entity,
                                                 query_dir_face, tmp_ori);
      tmp_enc = (r_corner != c_corner) * P8EST_EDGES + r_edge;
      insert_neighbor_elements (n_encs, n_qids, n_vids, tmp_enc, tmp_qid,
                                n_vid);
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
    sc_array_destroy (n_quads);
    return 0;
  }
  offset += P8EST_EDGES;
#endif /* P4_TO_P8 */
  if ((0 == n_encs->elem_count)
      && (offset <= dir && dir < (offset + P4EST_CHILDREN))
      && (mesh->quad_to_corner[P4EST_CHILDREN * qid + (dir - offset)] == -1)) {
    /* understand if we are hanging across an edge or a face (3D only, in 2D we
     * must be hanging across a face).  The following table shows the meaning
     * of the value in the sense that the corner is hanging on a face normal to
     * or on an edge parallel to the given direction.
     *
     * |-------+-----------+-----------+-----------|
     * | value | face (2D) | face (3D) | edge (3D) |
     * |-------+-----------+-----------+-----------|
     * |     1 | y-dir     |           | x-dir     |
     * |     2 | x-dir     |           | y-dir     |
     * |     3 |           | z-dir     |           |
     * |     4 |           |           | z-dir     |
     * |     5 |           | y-dir     |           |
     * |     6 |           | x-dir     |           |
     * |-------+-----------+-----------+-----------|
     */
    sid = p4est_quadrant_child_id (current);
    corner_type = sid ^ (dir - offset);
    /* normalize values for 2D and 3D such that
     *  1 -> face in x-dir
     *  2 -> face in y-dir
     *  4 -> face in z-dir
     */
    corner_type_norm = (P4EST_CHILDREN - 1) ^ corner_type;

    /* face hanging */
    if ((corner_type_norm == 1) || (corner_type_norm == 2)
#ifdef P4_TO_P8
        || (corner_type_norm == 4)
#endif /* P4_TO_P8 */
      ) {
      i = 0;
      while (corner_type_norm != 1) {
        corner_type_norm = corner_type_norm >> 1;
        ++i;
      }
      tmp_dir = get_adjacent_face (sid, i);
      get_corner_hanging_face (p4est, ghost, mesh, virtual_quads, qid,
                               tmp_dir, n_encs, n_qids, n_vids);
    }
#ifdef P4_TO_P8
    /* edge hanging
     * from normalization:
     *  3 -> z-dir (= xy-plane)
     *  5 -> y-dir (= xz-plane)
     *  6 -> x-dir (= yz-plane)
     */
    else if ((corner_type_norm == 3) || (corner_type_norm == 5)
             || (corner_type_norm == 6)) {
      i = 0;
      while (corner_type != 1) {
        corner_type = corner_type >> 1;
        ++i;
      }
      tmp_dir = get_adjacent_edge (sid, i);
      get_corner_hanging_edge (p4est, ghost, mesh, virtual_quads, qid,
                               dir - offset, tmp_dir, sid, n_encs, n_qids,
                               n_vids);
    }
#endif /* P4_TO_P8 */
    else {
      SC_ABORT_NOT_REACHED ();
    }
  }

  /* check output */
  P4EST_ASSERT (n_qids->elem_count == n_encs->elem_count);
  P4EST_ASSERT (n_vids->elem_count == n_encs->elem_count);
  sc_array_destroy (n_quads);

  return 0;
}

/** Find neighboring quadrants across faces querying from a virtual quadrant
 * All parameters are defined exactly as in \ref p4est_virtual_get_neighbor.
 */
static int
get_virtual_face_neighbors (p4est_t * p4est, p4est_ghost_t * ghost,
                            p4est_mesh_t * mesh,
                            p4est_virtual_t * virtual_quads,
                            p4est_locidx_t qid, p4est_locidx_t vid, int dir,
                            sc_array_t * n_encs, sc_array_t * n_qids,
                            sc_array_t * n_vids)
{
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  int                 l_same_size, u_same_size, l_double_size, u_double_size;
  int                 l_half_size, u_half_size;
#ifdef P4EST_ENABLE_DEBUG
  p4est_quadrant_t   *curr_quad = p4est_mesh_get_quadrant (p4est, mesh, qid);
  p4est_quadrant_t   *quad;
#endif /* P4EST_ENABLE_DEBUG */
  p4est_locidx_t      neighbor_idx, neighbor_enc, neighbor_vid,
    neighbor_subid;
  int                 neighbor_subquad, neighbor_orientation;
  int                 neighbor_entity_index;
  int                 tmp_dir;
  p4est_locidx_t      quad_idx;
  p4est_locidx_t     *quad_ptr;
  P4EST_ASSERT (0 <= dir && dir < P4EST_FACES);
  tmp_dir = p4est_virtual_face_neighbors_search_opts[vid][dir];
  /** internal neighbor */
  if (0 <= tmp_dir && tmp_dir < P4EST_CHILDREN) {
    neighbor_idx = qid;
    neighbor_enc = get_opposite (dir);
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, neighbor_idx, tmp_dir);
    return 0;
  }

#ifdef P4_TO_P8
  l_same_size = 0;
  u_same_size = 4 * P4EST_FACES;
  l_double_size = u_same_size;
  u_double_size = 120;
  l_half_size = -24;
  u_half_size = l_same_size;
#else /* P4_TO_P8 */
  l_same_size = 0;
  u_same_size = 2 * P4EST_FACES;
  l_double_size = u_same_size;
  u_double_size = 24;
  l_half_size = -8;
  u_half_size = l_same_size;
#endif /* P4_TO_P8 */
  /** perform external neighbor search */
  neighbor_idx = mesh->quad_to_quad[P4EST_FACES * qid + dir];
  neighbor_enc = mesh->quad_to_face[P4EST_FACES * qid + dir];
  /** no neighbor present */
  if ((neighbor_idx == qid) && neighbor_enc == dir) {
    return 0;
  }

  /** same size neighbor:
   * Virtual quadrants look for same sized neighbors if the neighboring
   * quadrant contains virtual quadrants. In this case we have to select the
   * neighbor's proper virtual subquadrant.
   */
  if (l_same_size <= neighbor_enc && neighbor_enc < u_same_size) {
    if (neighbor_idx < lq) {
      if (-1 == virtual_quads->virtual_qflags[neighbor_idx]) {
        return 0;
      }
    }
    else if (lq <= neighbor_idx && neighbor_idx < (lq + gq)) {
      if (-1 == virtual_quads->virtual_gflags[neighbor_idx - lq]) {
        return 0;
      }
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }

    /** decode encoding to obtain orientation and neighboring face index */
    p4est_mesh_decode_encoding (neighbor_enc, P4EST_FACES,
                                l_same_size, u_same_size,
                                l_double_size, u_double_size,
                                l_half_size, u_half_size,
                                &neighbor_subquad,
                                &neighbor_orientation,
                                &neighbor_entity_index);
    /** find neighboring subquad from own vid across face */
    neighbor_vid =
      p4est_connectivity_face_neighbor_corner (vid, dir,
                                               neighbor_entity_index,
                                               neighbor_orientation);
    P4EST_ASSERT (0 <= neighbor_vid && neighbor_vid < P4EST_CHILDREN);
    if (neighbor_idx < lq) {
#ifdef P4EST_ENABLE_DEBUG
      quad = p4est_mesh_get_quadrant (p4est, mesh, neighbor_idx);
#endif /* P4EST_ENABLE_DEBUG */
    }
    else {
#ifdef P4EST_ENABLE_DEBUG
      quad = p4est_quadrant_array_index (&ghost->ghosts, neighbor_idx - lq);
#endif /* P4EST_ENABLE_DEBUG */
    }
#ifdef P4EST_ENABLE_DEBUG
    /* sanity check level */
    P4EST_ASSERT (quad->level == curr_quad->level);
#endif /* P4EST_ENABLE_DEBUG */
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, neighbor_idx, neighbor_vid);
  }
  /** double size neighbor:
   * Only real quadrants are looking for neighbors here. Do not take any action
   */
  else if (l_double_size <= neighbor_enc && neighbor_enc < u_double_size) {
    NULL;
  }
  /* half size:
   * search among the P4EST_HALF neighboring quadrants for the one touching the
   * current virtual quadrant.
   */
  else if (l_half_size <= neighbor_enc && neighbor_enc < u_half_size) {
    quad_ptr =
      (p4est_locidx_t *) sc_array_index (mesh->quad_to_half, neighbor_idx);
    p4est_mesh_decode_encoding (neighbor_enc, P4EST_FACES, l_same_size,
                                u_same_size, l_double_size, u_double_size,
                                l_half_size, u_half_size, &neighbor_subquad,
                                &neighbor_orientation,
                                &neighbor_entity_index);
    /** find neighboring subquad from own vid across face:
     * Lookup face corner index and transform it across face */
    neighbor_subid = p4est_corner_face_corners[vid][dir];
    P4EST_ASSERT (0 <= neighbor_subid && neighbor_subid < P4EST_HALF);
    quad_idx = quad_ptr[neighbor_subid];
    neighbor_vid = -1;
    neighbor_enc = P4EST_FACES * neighbor_orientation + neighbor_entity_index;
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, quad_idx, neighbor_vid);
  }
  return 0;
}

#ifdef P4_TO_P8
/** Find neighboring quadrants across edges querying from a virtual quadrant
 * All parameters are defined exactly as in \ref p4est_virtual_get_neighbor
 */
static int
get_virtual_edge_neighbors (p4est_t * p4est, p4est_ghost_t * ghost,
                            p4est_mesh_t * mesh,
                            p4est_virtual_t * virtual_quads,
                            p4est_locidx_t qid, p4est_locidx_t vid,
                            int dir, sc_array_t * n_encs,
                            sc_array_t * n_qids, sc_array_t * n_vids)
{
  int                 i;
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  int                 l_same_size, u_same_size, l_double_size, u_double_size;
  int                 l_half_size, u_half_size;
  int                 u_double_size_face;
  int                 l_half_size_face, u_half_size_face;
  p4est_locidx_t      n_adj_quads, offset;
  p4est_locidx_t      quad_idx;
  p4est_locidx_t      neighbor_idx, neighbor_enc;
  int                 neighbor_vid;
  int                 tmp_dir, tmp_subindex;
  int                 tmp_qid, tmp_enc;
  int                 tmp_subquad, tmp_ori, tmp_entity;
  int                 n_vid;
  int                 r_edge, r_corner, c_corner;
  int                 opp_edge_offset, opp_corner_offset;
  int                 tmp_corner;
  int                 edge_offset, corner_offset;
  P4EST_ASSERT (0 <= dir && dir < P8EST_EDGES);
  tmp_dir = p8est_virtual_edge_neighbors_search_opts[vid][dir];
  /** find internal neighbor of a virtual quadrant */
  if (tmp_dir < P4EST_CHILDREN) {
    neighbor_enc = get_opposite (dir + P4EST_FACES);
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, qid, tmp_dir);
    return 0;
  }
  l_same_size = 0;
  u_same_size = 2 * P8EST_EDGES;
  l_double_size = u_same_size;
  u_double_size = 72;
  l_half_size = -24;
  u_half_size = l_same_size;
  u_double_size_face = 120;
  l_half_size_face = -24;
  u_half_size_face = l_same_size;
  /** adapt search direction to grab a specific face neighbor */
  if (P4EST_CHILDREN <= tmp_dir && tmp_dir < 32) {
    tmp_dir -= P8EST_CHILDREN;
    tmp_subindex = tmp_dir / P4EST_FACES;
    tmp_dir = tmp_dir % P4EST_FACES;
    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                              tmp_dir, NULL, n_encs, n_qids);
    /** select the correct quadrant among hanging quadrants */
    if (n_encs->elem_count == P4EST_HALF) {
      tmp_qid = *(int *) sc_array_index (n_qids, tmp_subindex);
      tmp_enc = *(int *) sc_array_index (n_encs, tmp_subindex);
      /** decode encoding */
      p4est_mesh_decode_encoding (tmp_enc, P4EST_FACES,
                                  l_same_size, u_same_size,
                                  l_double_size, u_double_size_face,
                                  l_half_size_face, u_half_size_face,
                                  &tmp_subquad, &tmp_ori, &tmp_entity);
      sc_array_truncate (n_qids);
      sc_array_truncate (n_encs);
      /** set xor operants for obtaining opposite edges and corners across the
       * respective face we were querying */
      set_xor_constants_edge (dir, tmp_dir, &edge_offset, &corner_offset);
      /** calculate orientation for edge:
       * 1. calculate index of edge in the neighboring quadrant
       * 2. calculate corner index of hanging quadrant in the neighboring
       *    quadrant
       * 3. transfer matching corner index across edge assuming o = 0;
       *  a) corner indices are equal => o = 0
       *  b) else o = 1
       */
      P4EST_ASSERT (0 <= tmp_subindex && tmp_subindex < P4EST_HALF);
      r_edge =
        p8est_connectivity_face_neighbor_edge (dir ^
                                               edge_offset,
                                               tmp_dir, tmp_entity, tmp_ori);
      r_corner =
        p4est_connectivity_face_neighbor_corner (vid ^
                                                 corner_offset,
                                                 tmp_dir,
                                                 tmp_entity, tmp_ori);
      set_xor_constants_edge (r_edge, tmp_entity,
                              &opp_edge_offset, &opp_corner_offset);
      tmp_corner =
        p8est_connectivity_edge_neighbor_corner (vid,
                                                 dir ^
                                                 edge_offset,
                                                 r_edge ^ opp_edge_offset, 0);
      tmp_ori = (tmp_corner != r_corner);
      P4EST_ASSERT ((0 == tmp_ori)
                    || (r_corner ==
                        p8est_connectivity_edge_neighbor_corner
                        (vid, dir ^ edge_offset,
                         r_edge ^ opp_edge_offset, 1)));
      tmp_enc = tmp_ori * P8EST_EDGES + r_edge;
      insert_neighbor_elements (n_encs, n_qids, n_vids, tmp_enc, tmp_qid, -1);
    }
    /* among same sized neighbor: select correct virtual quadrant if it hosts
     * any */
    else if (1 == n_encs->elem_count) {
      tmp_qid = *(int *) sc_array_index (n_qids, 0);
      P4EST_ASSERT (0 <= tmp_qid && tmp_qid < (lq + gq));
      tmp_enc = *(int *) sc_array_index (n_encs, 0);
      sc_array_truncate (n_qids);
      sc_array_truncate (n_encs);

      /** Take a look at the encoding:
       * If the neighboring quadrant has the same size and contains virtual
       * quadrants the current quadrant has a neighbor in the neighboring
       * quadrant.
       */
      if (l_same_size <= tmp_enc && tmp_enc < u_same_size
          && (((0 <= tmp_qid) && (tmp_qid < lq)
               && -1 != virtual_quads->virtual_qflags[tmp_qid])
              || ((lq <= tmp_qid) && (tmp_qid < (lq + gq))
                  && -1 != virtual_quads->virtual_gflags[tmp_qid - lq]))) {
        /* obtain vid from face corner transformation:
         * transform face corner index of same size sibling neighbor to large
         * face and lookup corner index */
        p4est_mesh_decode_encoding (tmp_enc, P4EST_FACES, l_same_size,
                                    u_same_size, l_double_size,
                                    u_double_size_face, l_half_size_face,
                                    u_half_size_face, &tmp_subquad, &tmp_ori,
                                    &tmp_entity);
        n_vid =
          p4est_connectivity_face_neighbor_face_corner (tmp_subindex, tmp_dir,
                                                        tmp_entity, tmp_ori);
        P4EST_ASSERT (0 <= n_vid && n_vid < P4EST_HALF);
        n_vid = p4est_face_corners[tmp_entity][n_vid];
        P4EST_ASSERT (0 <= n_vid && n_vid < P4EST_CHILDREN);

        /* Obtain encoding:
         * Two informations are missing:
         * 1) across which edge the current quadrant is seen from the neighbor
         * 2) are the edges flipped or do they have the same orientation
         * The first part can be obtained by transforming our edge into the
         * neighboring quadrant.  The index we are looking for is the opposite
         * edge across the double size neighbor's face.
         * The orientation can be determined by transforming the vid once across
         * the face with known orientation and once across the edge with assumed
         * orientation 0.  We guessed correctly if both indices match.
         */
        r_edge =
          p8est_connectivity_face_neighbor_edge (dir, tmp_dir, tmp_entity,
                                                 tmp_ori);
        r_corner =
          p8est_connectivity_edge_neighbor_corner (n_vid, r_edge, dir, 0);

        set_xor_constants_edge (r_edge, tmp_entity, &edge_offset,
                                &corner_offset);
        r_edge = r_edge ^ edge_offset;

        c_corner =
          p4est_connectivity_face_neighbor_corner (n_vid, tmp_entity, tmp_dir,
                                                   tmp_ori);
        tmp_enc = (r_corner != c_corner) * P8EST_EDGES + r_edge;
        insert_neighbor_elements (n_encs, n_qids, n_vids, tmp_enc, tmp_qid,
                                  n_vid);
      }
    }
    /* this is actually not a hack, because trees cannot hang.  Therefore, we do
     * not have to consider any special cases and invalidate a subset of the
     * result from neighbor lookup. */
    else {
      P4EST_ASSERT (0 == n_encs->elem_count);
    }
    return 0;
  }

  /* external querying across edge, i.e. retain original search direction */
  P4EST_ASSERT (32 <= tmp_dir && tmp_dir < 56);
  neighbor_idx = mesh->quad_to_edge[P8EST_EDGES * qid + dir];
  /** Domain boundary: No neighbor present. */
  if (neighbor_idx == -3) {
    return 0;
  }

  /** neighbor has same size as parent cell (and orientation 0 across tree
   * boundaries): Collect virtual quadrant */
  if (((0 <= neighbor_idx) && (neighbor_idx < lq)
       && (virtual_quads->virtual_qflags[neighbor_idx] != -1))
      || ((lq <= neighbor_idx)
          && (neighbor_idx < (lq + gq))
          && (virtual_quads->virtual_gflags[neighbor_idx - lq] != -1))) {

    neighbor_enc = dir ^ 3;
    neighbor_vid = p8est_edge_edge_corners[dir][vid];
    P4EST_ASSERT (neighbor_vid == 0 || neighbor_vid == 1);
    neighbor_vid = p8est_edge_corners[dir ^ 3][neighbor_vid];
    P4EST_ASSERT (0 <= neighbor_vid && neighbor_vid < P4EST_CHILDREN);
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, neighbor_idx, neighbor_vid);
    return 0;
  }
  else if ((lq + gq) <= neighbor_idx) {
    /* anything else: different sized neighbors and same size neighbors across
     * tree boundaries with flipped edges.
     * There are two relevant cases:
     *  1) the neighbor quadrant hosts virtual quadrants and has the same size
     *     as current virtual quadrant's hosts
     *  2) there are two half size neighboring quadrants adjacent to the current
     *     host quadrant
     */

    /* normalize neighbor index */
    neighbor_idx -= (lq + gq);
    P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < mesh->local_num_edges);
    /* get offset and number of adjacent quads as loop boundaries */
    offset = *((p4est_locidx_t *)
               sc_array_index_int (mesh->edge_offset, neighbor_idx));
    n_adj_quads = *((p4est_locidx_t *)
                    sc_array_index_int (mesh->edge_offset, neighbor_idx + 1));
    for (i = offset; i < n_adj_quads; ++i) {
      quad_idx = *((p4est_locidx_t *)
                   sc_array_index_int (mesh->edge_quad, i));
      P4EST_ASSERT (0 <= quad_idx && quad_idx < (lq + gq));
      neighbor_enc = *((int8_t *)
                       sc_array_index_int (mesh->edge_edge, i));
      /* inspect encoding of obtained neighbor */
      if ((l_same_size <= neighbor_enc && neighbor_enc < u_same_size)
          && (((0 <= quad_idx && quad_idx < lq)
               && virtual_quads->virtual_qflags[quad_idx] != -1)
              || ((lq <= quad_idx && quad_idx < (lq + gq))
                  && virtual_quads->virtual_gflags[quad_idx - lq]))) {
        /* decode it and get virtual id in neighboring quadrant */
        p4est_mesh_decode_encoding (neighbor_enc, P8EST_EDGES,
                                    l_same_size, u_same_size,
                                    l_double_size, u_double_size,
                                    l_half_size, u_half_size,
                                    &tmp_subquad, &tmp_ori, &tmp_entity);
        /** tmp_dir must be between 32 and 55 */
        tmp_dir -= 32;
        tmp_subindex = tmp_dir / P8EST_EDGES;
        /** find neighboring vid from own subquad ID */
        neighbor_vid =
          p8est_connectivity_edge_neighbor_edge_corner
          (tmp_subindex, tmp_ori);
        neighbor_vid = p8est_edge_corners[tmp_entity][neighbor_vid];
        P4EST_ASSERT (0 <= neighbor_vid && neighbor_vid < P4EST_CHILDREN);
        /** convert encoding */
        neighbor_enc = P8EST_EDGES * tmp_ori + tmp_entity;
        insert_neighbor_elements (n_encs, n_qids, n_vids,
                                  neighbor_enc, quad_idx, neighbor_vid);
      }
      else if (l_half_size <= neighbor_enc && neighbor_enc < u_half_size) {
        p4est_mesh_decode_encoding (neighbor_enc, P8EST_EDGES,
                                    l_same_size, u_same_size,
                                    l_double_size, u_double_size,
                                    l_half_size, u_half_size,
                                    &tmp_subquad, &tmp_ori, &tmp_entity);
        /* find neighboring subquad from own vid across edge */
        tmp_subquad = p8est_corner_edge_corners[vid][dir];
        P4EST_ASSERT (0 <= tmp_subquad && tmp_subquad < P4EST_HALF);
        tmp_subquad =
          p8est_connectivity_edge_neighbor_edge_corner (tmp_subquad, tmp_ori);
        P4EST_ASSERT (0 <= tmp_subquad && tmp_subquad < 2);
        neighbor_enc = P8EST_EDGES * tmp_ori + tmp_entity;
        neighbor_vid = -1;
        if (i - offset == tmp_subquad) {
          insert_neighbor_elements (n_encs, n_qids, n_vids,
                                    neighbor_enc, quad_idx, -1);
        }
      }
    }
  }

  return 0;
}
#endif /* P4_TO_P8 */

/** Find neighboring quadrants across corners querying from a virtual quadrant
 * All parameters are defined exactly as in \ref p4est_virtual_get_neighbor
 */
static int
get_virtual_corner_neighbors (p4est_t * p4est,
                              p4est_ghost_t * ghost,
                              p4est_mesh_t * mesh,
                              p4est_virtual_t *
                              virtual_quads,
                              p4est_locidx_t qid,
                              p4est_locidx_t vid,
                              int dir,
                              sc_array_t * n_encs,
                              sc_array_t * n_qids, sc_array_t * n_vids)
{
  int                 i, *int_ins;
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  p4est_quadrant_t   *curr_quad = p4est_mesh_get_quadrant (p4est, mesh, qid);
  p4est_quadrant_t   *quad;
  p4est_locidx_t      neighbor_idx, neighbor_vid, neighbor_enc;
  int                 tmp_dir, tmp_subindex;
  int                 tmp_subquad, tmp_ori, tmp_entity;
  int                 l_same_size_face, u_same_size_face;
  int                 l_double_size_face, u_double_size_face;
  int                 l_half_size_face, u_half_size_face;
#ifdef P4_TO_P8
  int                 l_same_size_edge, u_same_size_edge;
  int                 l_double_size_edge, u_double_size_edge;
  int                 l_half_size_edge, u_half_size_edge;
#endif /* P4_TO_P8 */
  int                 c_corner, offset;
  P4EST_ASSERT (0 <= dir && dir < P4EST_CHILDREN);

  /* TODO:
   * Change this according to ideas sketched above in get_neighbor_real.
   * Find the respective corner type and alter query accordlingly
   */
  tmp_dir = p4est_virtual_corner_neighbors_search_opts[vid][dir];
  /** find internal neighbor of a virtual quadrant */
  if (tmp_dir < P4EST_CHILDREN) {
#ifndef P4_TO_P8
    neighbor_enc = get_opposite (dir + P4EST_FACES);
#else /* P4_TO_P8 */
    neighbor_enc = get_opposite (dir + P8EST_EDGES + P4EST_FACES);
#endif /* P4_TO_P8 */
    insert_neighbor_elements (n_encs, n_qids, n_vids,
                              neighbor_enc, qid, tmp_dir);
    return 0;
  }

  /* set encoding boundaries to understand the encodings obtained in neighbor
   * search */
#ifndef P4_TO_P8
  l_half_size_face = -8;
  u_half_size_face = l_same_size_face = 0;
  u_same_size_face = l_double_size_face = 8;
  u_double_size_face = 24;
#else /* P4_TO_P8 */
  l_half_size_face = -24;
  u_half_size_face = l_same_size_face = 0;
  u_same_size_face = l_double_size_face = 24;
  u_double_size_face = 120;
  l_half_size_edge = -24;
  u_half_size_edge = l_same_size_edge = 0;
  u_same_size_edge = l_double_size_edge = 24;
  u_double_size_edge = 72;
#endif /* P4_TO_P8 */
  /* external query */
  /* search respective face neighbor */
  if (P4EST_CHILDREN <= tmp_dir
      && tmp_dir <
      (P4EST_CHILDREN + (P4EST_FACES * (1 << (P4EST_DIM - 1))))) {
    /* decode tmp_dir */
    tmp_dir -= P4EST_CHILDREN;
    tmp_subquad = tmp_dir / P4EST_FACES;
    tmp_dir = tmp_dir % P4EST_FACES;
    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                              tmp_dir, NULL, n_encs, n_qids);
    P4EST_ASSERT (0 <= tmp_subquad && tmp_subquad < P4EST_HALF);
    /* if we obtained a set of hanging quadrants: pass the correct quadrant and
     * calculate its encoding */
    if (n_qids->elem_count == P4EST_HALF) {
      neighbor_idx = *(int *) sc_array_index (n_qids, tmp_subquad);
      P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < (lq + gq));
      neighbor_vid = -1;
      neighbor_enc = *(int *) sc_array_index (n_encs, tmp_subquad);
      P4EST_ASSERT (l_half_size_face <= neighbor_enc
                    && neighbor_enc < u_half_size_face);

      p4est_mesh_decode_encoding (neighbor_enc, P4EST_FACES, l_same_size_face,
                                  u_same_size_face, l_double_size_face,
                                  u_double_size_face, l_half_size_face,
                                  u_half_size_face, &tmp_subindex, &tmp_ori,
                                  &tmp_entity);

      set_xor_constants_corner_inside (tmp_dir, &offset);
      c_corner = dir ^ offset;
      neighbor_enc =
        p4est_connectivity_face_neighbor_corner (c_corner, tmp_dir,
                                                 tmp_entity, tmp_ori);
      sc_array_truncate (n_encs);
      sc_array_truncate (n_qids);
      insert_neighbor_elements (n_encs, n_qids, n_vids,
                                neighbor_enc, neighbor_idx, neighbor_vid);
      return 0;
    }

    /* if we obtained a single quadrant it is only relevant iff it is a same
     * sized neighbor of the current host quadrant hosting virtual quadrants.
     */
    else if (n_qids->elem_count == 1) {
      neighbor_enc = *(int *) sc_array_index (n_encs, 0);
      neighbor_idx = *(int *) sc_array_index (n_qids, 0);
      P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < (lq + gq));
      if ((l_same_size_face <= neighbor_enc
           && neighbor_enc < u_same_size_face)
          && (((0 <= neighbor_idx && neighbor_idx < lq)
               && (-1 != virtual_quads->virtual_qflags[neighbor_idx]))
              || ((lq <= neighbor_idx && neighbor_idx < (lq + gq))
                  && (-1 !=
                      virtual_quads->virtual_gflags[neighbor_idx - lq])))) {

        p4est_mesh_decode_encoding (neighbor_enc, P4EST_FACES,
                                    l_same_size_face, u_same_size_face,
                                    l_double_size_face, u_double_size_face,
                                    l_half_size_face, u_half_size_face,
                                    &tmp_subindex, &tmp_ori, &tmp_entity);
        P4EST_ASSERT (tmp_subindex == -1);

        set_xor_constants_corner_inside (tmp_dir, &offset);

        c_corner = dir ^ offset;
        neighbor_enc =
          p4est_connectivity_face_neighbor_corner (c_corner, tmp_dir,
                                                   tmp_entity, tmp_ori);
        neighbor_vid =
          p4est_connectivity_face_neighbor_face_corner (tmp_subquad, tmp_dir,
                                                        tmp_entity, tmp_ori);
        P4EST_ASSERT (0 <= neighbor_vid && neighbor_vid < P4EST_HALF);
        neighbor_vid = p4est_face_corners[tmp_entity][neighbor_vid];
        P4EST_ASSERT (0 <= neighbor_vid && neighbor_vid < P4EST_CHILDREN);

        sc_array_truncate (n_encs);
        sc_array_truncate (n_qids);
        insert_neighbor_elements (n_encs, n_qids, n_vids,
                                  neighbor_enc, neighbor_idx, neighbor_vid);
        return 0;
      }
      else {
        sc_array_truncate (n_encs);
        sc_array_truncate (n_qids);
        return 0;
      }
    }
  }

#ifdef P4_TO_P8
  /* search respective edge neighbor (hard-wire values, no need for dimension
   * independent expression as this is 3D anyways) */
  else if (32 <= tmp_dir && tmp_dir < 56) {
    /* decode tmp_dir */
    tmp_dir -= 32;
    tmp_subquad = tmp_dir / P8EST_EDGES;
    tmp_dir = tmp_dir % P8EST_EDGES;
    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                              P4EST_FACES + tmp_dir, NULL, n_encs, n_qids);
    P4EST_ASSERT (0 <= tmp_subquad && tmp_subquad < 2);
    /* if we obtained a set of hanging quadrants: pass the correct quadrant and
     * calculate its encoding */
    if (n_qids->elem_count == 2) {
      neighbor_idx = *(int *) sc_array_index (n_qids, tmp_subquad);
      P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < (lq + gq));
      neighbor_enc = *(int *) sc_array_index (n_encs, tmp_subquad);

      P4EST_ASSERT (l_half_size_edge <= neighbor_enc
                    && neighbor_enc < u_half_size_edge);

      p4est_mesh_decode_encoding (neighbor_enc, P8EST_EDGES, l_same_size_edge,
                                  u_same_size_edge, l_double_size_edge,
                                  u_double_size_edge, l_half_size_edge,
                                  u_half_size_edge, &tmp_subindex, &tmp_ori,
                                  &tmp_entity);

      /* We see the current quadrant through the opposite corner of the current
       * edge.  We can obtain that corner e.g. by virtually flipping the edge.
       */
      neighbor_enc =
        p8est_connectivity_edge_neighbor_edge_corner (tmp_subquad,
                                                      (tmp_ori + 1) % 2);
      neighbor_enc = p8est_edge_corners[tmp_entity][neighbor_enc];

      neighbor_vid = -1;
      sc_array_truncate (n_encs);
      sc_array_truncate (n_qids);
      insert_neighbor_elements (n_encs, n_qids, n_vids,
                                neighbor_enc, neighbor_idx, neighbor_vid);
      return 0;
    }

    /* if we obtained a single quadrant it is only relevant iff it is a same
     * sized neighbor of the current host quadrant hosting virtual quadrants.
     */
    else if (n_qids->elem_count == 1) {
      neighbor_enc = *(int *) sc_array_index (n_encs, 0);
      neighbor_idx = *(int *) sc_array_index (n_qids, 0);
      P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < (lq + gq));
      if ((l_same_size_edge <= neighbor_enc
           && neighbor_enc < u_same_size_edge)
          && (((0 <= neighbor_idx && neighbor_idx < lq)
               && (-1 != virtual_quads->virtual_qflags[neighbor_idx]))
              || ((lq <= neighbor_idx && neighbor_idx < (lq + gq))
                  && (-1 !=
                      virtual_quads->virtual_gflags[neighbor_idx - lq])))) {

        p4est_mesh_decode_encoding (neighbor_enc, P8EST_EDGES,
                                    l_same_size_face, u_same_size_face,
                                    l_double_size_face, u_double_size_face,
                                    l_half_size_face, u_half_size_face,
                                    &tmp_subindex, &tmp_ori, &tmp_entity);
        /* We see the current quadrant through the opposite corner of the current
         * edge.  We can obtain that corner e.g. by virtually flipping the edge.
         */
        neighbor_enc =
          p8est_connectivity_edge_neighbor_edge_corner (tmp_subquad,
                                                        (tmp_ori + 1) % 2);
        neighbor_enc = p8est_edge_corners[tmp_entity][neighbor_enc];

        neighbor_vid =
          p8est_connectivity_edge_neighbor_edge_corner (tmp_subquad, tmp_ori);
        neighbor_vid = p8est_edge_corners[tmp_entity][neighbor_vid];
        sc_array_truncate (n_encs);
        sc_array_truncate (n_qids);
        insert_neighbor_elements (n_encs, n_qids, n_vids,
                                  neighbor_enc, neighbor_idx, neighbor_vid);
        return 0;
      }
      else {
        sc_array_truncate (n_encs);
        sc_array_truncate (n_qids);

        return 0;
      }
    }
  }
#endif /* P4_TO_P8 */

  /* No need to change search direction: perform neighbor lookup across corner.
   */
  neighbor_idx = mesh->quad_to_corner[P4EST_CHILDREN * qid + dir];
  /* no neighbor present or hanging neighbor: irrelevant for virtual quadrants.
   */
  if (neighbor_idx < 0) {
    return 0;
  }

  P4EST_ASSERT (0 <= neighbor_idx
                && neighbor_idx < (lq + gq + mesh->local_num_corners));
  /* single diagonally-opposite neighbor, no rotation involved, i.e. implicit
   * encoding is used.  Again, we are interested in smaller neighbors and same
   * size neighbors hosting virtual quadrants.*/
  if (0 <= neighbor_idx && neighbor_idx < (lq + gq)) {
    if (0 <= neighbor_idx && neighbor_idx < lq) {
      quad = p4est_mesh_get_quadrant (p4est, mesh, neighbor_idx);
    }
    else if (lq <= neighbor_idx && neighbor_idx < (lq + gq)) {
      quad = p4est_quadrant_array_index (&ghost->ghosts, neighbor_idx - lq);
    }
    else {
      SC_ABORT_NOT_REACHED ();
    }
    if ((curr_quad->level == quad->level)
        && (((0 <= neighbor_idx && neighbor_idx < lq)
             && (virtual_quads->virtual_qflags[neighbor_idx] != -1))
            || ((lq <= neighbor_idx && neighbor_idx < (lq + gq))
                && (virtual_quads->virtual_gflags[neighbor_idx - lq] != -1)))) {
      neighbor_vid = neighbor_enc = dir ^ (P4EST_CHILDREN - 1);
      insert_neighbor_elements (n_encs, n_qids, n_vids,
                                neighbor_enc, neighbor_idx, neighbor_vid);
    }
    else if ((curr_quad->level + 1) == quad->level) {
      neighbor_vid = -1;
      neighbor_enc = dir ^ (P4EST_CHILDREN - 1);
      insert_neighbor_elements (n_encs, n_qids, n_vids,
                                neighbor_enc, neighbor_idx, neighbor_vid);
    }
    return 0;
  }
  /* Anything else.  This must be a tree-boundary. */
  else
    if ((lq + gq) <= neighbor_idx
        && neighbor_idx < (lq + gq + mesh->local_num_corners)) {
    offset = P4EST_FACES;
#ifdef P4_TO_P8
    offset += P8EST_EDGES;
#endif /* P4_TO_P8 */
    p4est_mesh_get_neighbors (p4est, ghost, mesh, qid,
                              dir + offset, NULL, n_encs, n_qids);
    for (i = 0; i < n_encs->elem_count; ++i) {
      /* FIXME: Remove this assertion and handle this case in code */
      P4EST_ASSERT (1 == n_encs->elem_count);
      neighbor_idx = *(int *) sc_array_index (n_qids, i);
      neighbor_enc = *(int *) sc_array_index (n_encs, i);
      P4EST_ASSERT (0 <= neighbor_idx && neighbor_idx < (lq + gq));
      quad =
        neighbor_idx < lq ? p4est_mesh_get_quadrant (p4est,
                                                     mesh, neighbor_idx)
        : p4est_quadrant_array_index (&ghost->ghosts, neighbor_idx - lq);
      if ((curr_quad->level == quad->level)
          && (((0 <= neighbor_idx && neighbor_idx < lq)
               && (virtual_quads->virtual_qflags[neighbor_idx] != -1))
              || ((lq <= neighbor_idx && neighbor_idx < (lq + gq))
                  && (virtual_quads->virtual_gflags[neighbor_idx - lq] !=
                      -1)))) {
        neighbor_vid = neighbor_enc;
        int_ins = (int *) sc_array_push (n_vids);
        *int_ins = neighbor_vid;
      }
      else if ((curr_quad->level + 1) == quad->level) {
        neighbor_vid = -1;
        int_ins = (int *) sc_array_push (n_vids);
        *int_ins = neighbor_vid;
      }
      else {
        sc_array_truncate (n_qids);
        sc_array_truncate (n_encs);
      }
    }
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  return 0;
}

/** Neighbor lookup for virtual quadrants that will only return same sized
 * quadrants as neighbors, real or virtual. That means the comparison from the
 * regular neighbor lookup \ref p4est_mesh_get_neighbors w.r.t. the quadrant's
 * level takes place between the neighbor and the host quadrant, i.e. a
 * same-size neighbor is actually a double-sized neighbor w.r.t. a virtual
 * quadrant.
 * \param[in]      p4est          The forest.
 * \param[in]      ghost          Ghost layer
 * \param[in]      mesh           Neighbor information of real quadrants.
 * \param[in]      virtual_quads  Information which quadrants host virtual
 *                                quadrants.
 * \param[in]      qid            Local quadrant-id of the quadrant whose
 *                                neighbors we are searching.
 * \param[in]      dir            Direction we search for neighbors.
 *                                2D:
 *                                  0 ..  3 neighbor(-s) across face i,
 *                                  4 ..  7 neighbor(-s) across corner i-4.
 *                                3D:
 *                                  0 ..  5 neighbor(-s) across face i,
 *                                  6 .. 17 neighbor(-s) across edge i-6
 *                                 18 .. 25 neighbor(-s) across corner i-18.
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
static int
get_neighbor_virtual (p4est_t * p4est,
                      p4est_ghost_t * ghost,
                      p4est_mesh_t * mesh,
                      p4est_virtual_t *
                      virtual_quads,
                      p4est_locidx_t qid, int vid,
                      int dir, sc_array_t * n_encs,
                      sc_array_t * n_qids, sc_array_t * n_vids)
{
  int                 temp_dir = 0;
  int                 limit;
  p4est_locidx_t      lq = mesh->local_num_quadrants;
#ifdef P4EST_ENABLE_DEBUG
  /* Integrity checks: */
  /* result arrays should be empty, */
  P4EST_ASSERT (n_encs == NULL || n_encs->elem_count == 0);
  P4EST_ASSERT (n_qids == NULL || n_qids->elem_count == 0);
  P4EST_ASSERT (n_vids == NULL || n_vids->elem_count == 0);
  /* result arrays need to have matching element sizes */
  P4EST_ASSERT (n_encs == NULL || n_encs->elem_size == sizeof (int));
  P4EST_ASSERT (n_qids == NULL || n_qids->elem_size == sizeof (int));
  P4EST_ASSERT (n_vids == NULL || n_vids->elem_size == sizeof (int));
  /* mesh and virtual quadrants have to be created, i.e. not NULL, */
  P4EST_ASSERT (mesh != NULL);
  P4EST_ASSERT (virtual_quads != NULL);
  /* qid must identify a local quadrant .. */
  P4EST_ASSERT (0 <= qid && qid < lq);
  /* and dir must be within the allowed range */
  switch (virtual_quads->btype) {
  case P4EST_CONNECT_FACE:
    P4EST_ASSERT (0 <= dir && dir < P4EST_FACES);
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    limit = P4EST_FACES + P8EST_EDGES;
    P4EST_ASSERT (0 <= dir && dir < limit);
    break;
#endif /* P4_TO_P8 */
  case P4EST_CONNECT_CORNER:
#ifdef P4_TO_P8
    limit = P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN;
#else /* P4_TO_P8 */
    limit = P4EST_FACES + P4EST_CHILDREN;
#endif /* P4_TO_P8 */
    P4EST_ASSERT (0 <= dir && dir < limit);
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }
#endif /* P4EST_ENABLE_DEBUG */

#ifdef P4_TO_P8
  temp_dir = P8EST_EDGES;
#endif /* P4_TO_P8 */
  /* inspect what to query */
  if (0 <= dir && dir < P4EST_FACES) {
    get_virtual_face_neighbors (p4est, ghost, mesh, virtual_quads, qid, vid,
                                dir, n_encs, n_qids, n_vids);
  }
#ifdef P4_TO_P8
  else if (P4EST_FACES <= dir && dir < P4EST_FACES + P8EST_EDGES) {
    get_virtual_edge_neighbors (p4est, ghost, mesh, virtual_quads, qid, vid,
                                dir - P4EST_FACES, n_encs, n_qids, n_vids);
  }
  else if (P4EST_FACES + P8EST_EDGES <= dir
           && dir < P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)
#else /* P4_TO_P8 */
  else if (P4EST_FACES <= dir && dir < P4EST_FACES + P4EST_CHILDREN)
#endif /* P4_TO_P8 */
  {
    get_virtual_corner_neighbors (p4est, ghost, mesh, virtual_quads, qid, vid,
                                  dir - P4EST_FACES - temp_dir, n_encs,
                                  n_qids, n_vids);
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
  return 0;
}

int
p4est_virtual_get_neighbor (p4est_t * p4est, p4est_ghost_t * ghost,
                            p4est_mesh_t * mesh,
                            p4est_virtual_t * virtual_quads,
                            p4est_locidx_t qid, int vid, int dir,
                            sc_array_t * n_encs, sc_array_t * n_qids,
                            sc_array_t * n_vids)
{
  int                 dir_max;
  switch (virtual_quads->btype) {
  case P4EST_CONNECT_FACE:
    dir_max = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    dir_max = P4EST_FACES + P8EST_EDGES;
    break;
#endif /* P4_TO_P8 */
  /* *INDENT-OFF* */
  case P4EST_CONNECT_FULL:
    dir_max = P4EST_FACES +
#ifdef P4_TO_P8
              P8EST_EDGES +
#endif /* P4_TO_P8 */
              P4EST_CHILDREN;
  /* *INDENT-ON* */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /* check input conditions */
  P4EST_ASSERT (0 <= qid && qid < p4est->local_num_quadrants);
  P4EST_ASSERT (-1 <= vid && vid < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= dir && dir < dir_max);
  P4EST_ASSERT (0 == n_encs->elem_count && n_encs->elem_size == sizeof (int));
  P4EST_ASSERT (0 == n_qids->elem_count
                && n_qids->elem_size == sizeof (p4est_locidx_t));
  P4EST_ASSERT (0 == n_vids->elem_count && n_vids->elem_size == sizeof (int));
  /* delegate neighbor lookup depending on type of quadrant */
  if (-1 == vid) {
    get_neighbor_real (p4est, ghost, mesh, virtual_quads,
                       qid, dir, n_encs, n_qids, n_vids);
  }
  else {
    get_neighbor_virtual (p4est, ghost, mesh,
                          virtual_quads, qid, vid, dir,
                          n_encs, n_qids, n_vids);
  }

  return 0;
}
