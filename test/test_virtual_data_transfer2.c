/*
 *  This file is part of p4est.
 *  p4est is a C library to manage a collection (a forest) of multiple
 *  connected adaptive quadtrees or octrees in parallel.
 *
 *  Copyright (C) 2010 The University of Texas System
 *  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac
 *
 *  p4est is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  p4est is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with p4est; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_virtual.h>
#else
#include <p8est_extended.h>
#include "p8est_virtual.h"
#endif

void
setup_p4est_structs (sc_MPI_Comm mpicomm, const int level,
                     p4est_connectivity_t ** conn, p4est_t ** p4est,
                     p4est_ghost_t ** ghost, p4est_mesh_t ** mesh,
                     p4est_virtual_t ** vq, p4est_virtual_ghost_t ** vg)
{
  const int           min_quadrants = 0;
  const int           fill_uniform = 0;
  const size_t        data_size = 0;
  const p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  const int           mesh_tree_indices = 1;
  const int           level_lists = 1;
  const int           mesh_par_boundary = 1;
  *conn = p4est_connectivity_new_periodic ();
  *p4est = p4est_new_ext (mpicomm, *conn, min_quadrants, level, fill_uniform,
                          data_size, NULL, NULL);
  p4est_balance (*p4est, btype, NULL);
  *ghost = p4est_ghost_new (*p4est, btype);
  *mesh = p4est_mesh_new_ext (*p4est, *ghost, mesh_tree_indices, level_lists,
                              mesh_par_boundary, btype);
  *vq = p4est_virtual_new_ext (*p4est, *ghost, *mesh, btype, level_lists);
  *vg = p4est_virtual_ghost_new (*p4est, *ghost, *mesh, *vq, btype);
}

void
free_p4est_structs (p4est_connectivity_t * conn, p4est_t * p4est,
                    p4est_ghost_t * ghost, p4est_mesh_t * mesh,
                    p4est_virtual_t * vq, p4est_virtual_ghost_t * vg)
{
  p4est_virtual_ghost_destroy (vg);
  p4est_virtual_destroy (vq);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

typedef struct message_content
{
  int                 foo;
  double              bar;
} message_content_t;

void
perform_test_flat (p4est_connectivity_t * conn, p4est_t * p4est,
                   p4est_ghost_t * ghost, p4est_mesh_t * mesh,
                   p4est_virtual_t * vq, p4est_virtual_ghost_t * vg)
{
  int                 qid, vid;
  int                 ctr = 0;
  const p4est_locidx_t lq = vq->local_num_quadrants;
  const p4est_locidx_t lv = vq->local_num_virtuals;
  const p4est_locidx_t gq = vq->ghost_num_quadrants;
  const p4est_locidx_t gv = vq->ghost_num_virtuals;

  message_content_t  *flat_data_local =
    P4EST_ALLOC (message_content_t, lq + lv);
  message_content_t  *flat_data_ghost =
    P4EST_ALLOC (message_content_t, gq + gv);

  for (qid = 0; qid < lq; ++qid) {
    flat_data_local[ctr].foo = qid;
    flat_data_local[ctr].bar = 0.02;
    ++ctr;
    if (vq->virtual_qflags[qid] != -1) {
      for (vid = 0; vid < P4EST_CHILDREN; ++vid) {
        flat_data_local[ctr].foo = -qid;
        flat_data_local[ctr].bar = 0.06;
        ++ctr;
      }
    }
  }

  p4est_virtual_ghost_exchange_data (p4est, ghost, mesh, vq, vg,
                                     sizeof (message_content_t),
                                     flat_data_local, flat_data_ghost);

  P4EST_FREE (flat_data_ghost);
  P4EST_FREE (flat_data_local);
}

void
populate_per_level_data (p4est_t * p4est, p4est_mesh_t * mesh,
                         p4est_virtual_t * vq,
                         message_content_t ** level_data_local,
                         message_content_t ** level_data_ghost)
{
  int                 qid, vid;
  p4est_quadrant_t   *q;
  const p4est_locidx_t lq = vq->local_num_quadrants;

  for (qid = 0; qid < lq; ++qid) {
    q = p4est_mesh_get_quadrant (p4est, mesh, qid);
    level_data_local[q->level][vq->quad_qreal_offset[qid]].foo = qid;
    level_data_local[q->level][vq->quad_qreal_offset[qid]].bar = 0.42;
    if (vq->virtual_qflags[qid] != -1) {
      for (vid = 0; vid < P4EST_CHILDREN; ++vid) {
        level_data_local[q->level + 1][vq->quad_qvirtual_offset[qid] + vid]
          .foo = -qid;
        level_data_local[q->level + 1][vq->quad_qvirtual_offset[qid] + vid]
          .bar = 0.17;
      }
    }
  }
}

void
perform_test_level (p4est_connectivity_t * conn, p4est_t * p4est,
                    p4est_ghost_t * ghost, p4est_mesh_t * mesh,
                    p4est_virtual_t * vq, p4est_virtual_ghost_t * vg)
{
  int                 level;
  size_t              quads_per_level;

  message_content_t **level_data_local =
    P4EST_ALLOC (message_content_t *, P4EST_QMAXLEVEL);
  message_content_t **level_data_ghost =
    P4EST_ALLOC (message_content_t *, P4EST_QMAXLEVEL);

  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    quads_per_level =
      (mesh->quad_level + level)->elem_count +
      P4EST_CHILDREN * (vq->virtual_qlevels + level)->elem_count;
    level_data_local[level] = P4EST_ALLOC (message_content_t, quads_per_level);
    quads_per_level =
      (mesh->ghost_level + level)->elem_count +
      P4EST_CHILDREN * (vq->virtual_glevels + level)->elem_count;
    level_data_ghost[level] = P4EST_ALLOC (message_content_t, quads_per_level);
  }

  populate_per_level_data (p4est, mesh, vq, level_data_local,
                           level_data_ghost);

  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    p4est_virtual_ghost_exchange_data_level (p4est, ghost, mesh, vq, vg,
                                             level,
                                             sizeof (message_content_t),
                                             (void **) level_data_local,
                                             (void **) level_data_ghost);
  }

  for (level = 0; level < P4EST_QMAXLEVEL; ++level) {
    P4EST_FREE (level_data_local[level]);
    P4EST_FREE (level_data_ghost[level]);
  }
  P4EST_FREE (level_data_local);
  P4EST_FREE (level_data_ghost);
}

void
execute_tests (sc_MPI_Comm mpicomm)
{
  p4est_connectivity_t *conn;
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virt_quads;
  p4est_virtual_ghost_t *virt_ghost;
  const int           level = 6;

  setup_p4est_structs (mpicomm, level, &conn, &p4est, &ghost, &mesh,
                       &virt_quads, &virt_ghost);

  perform_test_flat (conn, p4est, ghost, mesh, virt_quads, virt_ghost);
  perform_test_level (conn, p4est, ghost, mesh, virt_quads, virt_ghost);

  free_p4est_structs (conn, p4est, ghost, mesh, virt_quads, virt_ghost);
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  MPI_Comm_set_errhandler (mpicomm, MPI_ERRORS_RETURN);

  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  execute_tests (mpicomm);

  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
