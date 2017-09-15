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
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#include <p4est_virtual.h>
#else /* !P4_TO_P8 */
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_virtual.h>
#endif /* !P4_TO_P8 */

int
check_flat_exchange (p4est_t * p4est, p4est_ghost_t * ghost,
                     p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                     p4est_virtual_ghost_t * virtual_ghost)
{
  int                 lq, lv, gq, gv;
  int                 offset;
  int                 i, j;
  int                *local_buffer;
  int                *ghost_buffer;
  p4est_quadrant_t   *ghost_quad;

  lq = virtual_quads->local_num_quadrants;
  lv = virtual_quads->local_num_virtuals;
  gq = virtual_quads->ghost_num_quadrants;
  gv = virtual_quads->ghost_num_virtuals;

  local_buffer = P4EST_ALLOC_ZERO (int, lq + lv);
  ghost_buffer = P4EST_ALLOC_ZERO (int, gq + gv);

  /* populate local buffer:
   * All regular quadrants obtain their qid as payload.  Virtual quadrants'
   * payload is set to -qid.
   */
  for (offset = 0, i = 0; i < lq; ++i) {
    local_buffer[i + offset] = i;
    if (-1 < virtual_quads->virtual_qflags[i]) {
      for (j = 1; j <= P4EST_CHILDREN; ++j) {
        local_buffer[i + offset + j] = -i;
      }
      offset += P4EST_CHILDREN;
    }
  }

  /* perform ghost-exchange */
  p4est_virtual_ghost_exchange_data (p4est, ghost, mesh, virtual_quads,
                                     virtual_ghost, sizeof (int),
                                     local_buffer, ghost_buffer);

  for (offset = 0, i = 0; i < gq; ++i) {
    ghost_quad = p4est_quadrant_array_index (&ghost->ghosts, i);
    P4EST_ASSERT (ghost_buffer[i + offset] == ghost_quad->p.piggy3.local_num);
    if (-1 < virtual_quads->virtual_gflags[i]) {
      for (j = 1; j <= P4EST_CHILDREN; ++j) {
        P4EST_ASSERT (-1 * ghost_buffer[i + offset + j] ==
                      ghost_quad->p.piggy3.local_num);
      }
      offset += P4EST_CHILDREN;
    }
  }
  P4EST_ASSERT (offset == gv);

  P4EST_FREE (local_buffer);
  P4EST_FREE (ghost_buffer);

  return 0;
}

int
check_level_exchange (p4est_t * p4est, p4est_ghost_t * ghost,
                      p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads,
                      p4est_virtual_ghost_t * virtual_ghost)
{
  int                 offset;
  int                 level;
  int                 i, i_real, i_virt, j;
  p4est_locidx_t      qid_real, qid_virt;
  int               **local_buffer;
  int               **ghost_buffer;
  p4est_quadrant_t   *ghost_quad;

  local_buffer = P4EST_ALLOC (int *, P4EST_QMAXLEVEL + 1);
  ghost_buffer = P4EST_ALLOC (int *, P4EST_QMAXLEVEL + 1);

  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    local_buffer[level] =
      P4EST_ALLOC_ZERO (int,
                        (mesh->quad_level + level)->elem_count +
                        P4EST_CHILDREN * (virtual_quads->virtual_qlevels +
                                          level)->elem_count);
    ghost_buffer[level] =
      P4EST_ALLOC_ZERO (int,
                        (mesh->ghost_level + level)->elem_count +
                        P4EST_CHILDREN * (virtual_quads->virtual_glevels +
                                          level)->elem_count);
  }

  /* populate local buffer:
   * All regular quadrants obtain their qid as payload.  Virtual quadrants'
   * payload is set to -qid.
   */
  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    for (i_real = 0, i_virt = 0, i = 0;
         i < ((mesh->quad_level + level)->elem_count +
              (virtual_quads->virtual_qlevels + level)->elem_count); ++i) {
      if (i_real < (mesh->quad_level + level)->elem_count) {
        qid_real =
          *(p4est_locidx_t *) sc_array_index (mesh->quad_level + level,
                                              i_real);
      }
      else {
        qid_real = INT_MAX;
      }
      if (i_virt < (virtual_quads->virtual_qlevels + level)->elem_count) {
        qid_virt =
          *(p4est_locidx_t *) sc_array_index (virtual_quads->virtual_qlevels +
                                              level, i_virt);
      }
      else {
        qid_virt = INT_MAX;
      }
      if (qid_real < qid_virt) {
        offset = virtual_quads->quad_qreal_offset[qid_real];
        P4EST_ASSERT (local_buffer[level][offset] == 0);
        local_buffer[level][offset] = qid_real;
        ++i_real;
      }
      else if (qid_virt < qid_real) {
        offset = virtual_quads->quad_qvirtual_offset[qid_virt];
        P4EST_ASSERT (-1 != offset);
        for (j = 0; j < P4EST_CHILDREN; ++j) {
          P4EST_ASSERT (local_buffer[level][offset + j] == 0);
          local_buffer[level][offset + j] = -qid_virt;
        }
        ++i_virt;
      }
      else {
        P4EST_ASSERT (qid_real == qid_virt);
        P4EST_ASSERT (qid_real == INT_MAX);
      }
    }
  }

  /* perform ghost-exchange */
  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    p4est_virtual_ghost_exchange_data_level (p4est, ghost, mesh,
                                             virtual_quads, virtual_ghost,
                                             level, sizeof (int),
                                             (void **) local_buffer,
                                             (void **) ghost_buffer);
  }

  /* check that the correct data was received */
  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    for (i = 0, i_real = 0, i_virt = 0;
         i < ((mesh->ghost_level + level)->elem_count +
              (virtual_quads->virtual_glevels + level)->elem_count); ++i) {
      if (i_real < (mesh->ghost_level + level)->elem_count) {
        qid_real =
          *(p4est_locidx_t *) sc_array_index (mesh->ghost_level + level,
                                              i_real);
        ++i_real;
      }
      else {
        qid_real = INT_MAX;
      }
      if (i_virt < (virtual_quads->virtual_glevels + level)->elem_count) {
        qid_virt =
          *(p4est_locidx_t *) sc_array_index (virtual_quads->virtual_glevels +
                                              level, i_virt);
        ++i_virt;
      }
      else {
        qid_virt = INT_MAX;
      }
      if (qid_real < qid_virt) {
        ghost_quad = p4est_quadrant_array_index (&ghost->ghosts, qid_real);
        offset = virtual_quads->quad_greal_offset[qid_real];
        P4EST_ASSERT (ghost_quad->level == level);
        P4EST_ASSERT (ghost_buffer[level][offset] ==
                      ghost_quad->p.piggy3.local_num);
      }
      else if (qid_virt < qid_real) {
        ghost_quad = p4est_quadrant_array_index (&ghost->ghosts, qid_virt);
        offset = virtual_quads->quad_gvirtual_offset[qid_virt];
        P4EST_ASSERT (ghost_quad->level + 1 == level);
        for (j = 0; j < P4EST_CHILDREN; ++j) {
          P4EST_ASSERT (ghost_buffer[level][offset + j] ==
                        -ghost_quad->p.piggy3.local_num);
        }
      }
      else {
        P4EST_ASSERT (qid_real == qid_virt);
        P4EST_ASSERT (qid_real == INT_MAX);
      }
    }
  }

  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    P4EST_FREE (local_buffer[level]);
    P4EST_FREE (ghost_buffer[level]);
  }
  P4EST_FREE (local_buffer);
  P4EST_FREE (ghost_buffer);

  return 0;
}

int
check_virtual_ghost (sc_MPI_Comm mpicomm)
{
  p4est_locidx_t     *received_mirror_flags;
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  p4est_virtual_ghost_t *virtual_ghost;
  int                 min_level = 3;
  sc_MPI_Request     *r;
  int                 i, ng, n_proc, mpiret;
  int                 offset_old, offset_new;
  sc_array_t         *requests;

  /* setup p4est */
  conn = p4est_connectivity_new_periodic ();
  p4est = p4est_new_ext (mpicomm, conn, 0, min_level, 0, 0, NULL, NULL);
  p4est_balance (p4est, btype, NULL);

  n_proc = p4est->mpisize;

  /* setup anything else */
  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  virtual_quads = p4est_virtual_new_ext (p4est, ghost, mesh, btype, 1);
  virtual_ghost =
    p4est_virtual_ghost_new (p4est, ghost, mesh, virtual_quads, btype);

  /** Exchange information if ghost quadrants contain virtual quadrants on
   * local processor. */
  received_mirror_flags =
    P4EST_ALLOC (p4est_locidx_t, ghost->mirror_proc_offsets[n_proc]);
  requests = sc_array_new (sizeof (sc_MPI_Request));

  /** receive virtual_gflags from neighboring processes */
  offset_old = 0;
  for (i = 0; i < n_proc; ++i) {
    offset_new = ghost->mirror_proc_offsets[i + 1];
    ng = offset_new - offset_old;
    P4EST_ASSERT (0 <= ng);
    if (0 < ng) {
      r = (sc_MPI_Request *) sc_array_push (requests);
      mpiret = sc_MPI_Irecv ((void *) (received_mirror_flags + offset_old),
                             ng * sizeof (p4est_locidx_t), sc_MPI_BYTE, i,
                             P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      offset_old = offset_new;
    }
  }

  /** send virtual_gflags to neighboring processes:
   * Iterate over mpisize and check proc_offsets. If the current offset differs
   * from the offset before this is a processor from which we will receive
   * messages. Send to that processor if the current process expects virtual
   * quadrants during ghost exchange.
   */
  offset_old = 0;
  for (i = 0; i < n_proc; ++i) {
    offset_new = ghost->proc_offsets[i + 1];
    ng = offset_new - offset_old;
    P4EST_ASSERT (0 <= ng);
    if (0 < ng) {
      r = (sc_MPI_Request *) sc_array_push (requests);
      mpiret =
        sc_MPI_Isend ((void *) (virtual_quads->virtual_gflags + offset_old),
                      ng * sizeof (p4est_locidx_t), sc_MPI_BYTE, i,
                      P4EST_COMM_GHOST_EXCHANGE, p4est->mpicomm, r);
      SC_CHECK_MPI (mpiret);
      offset_old = offset_new;
    }
  }

  /** Wait for communication to finish */
  mpiret =
    sc_MPI_Waitall (requests->elem_count, (sc_MPI_Request *) requests->array,
                    sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);
  sc_array_destroy (requests);

  /* compare obtained array with the locally taken decision */
  for (i = 0; i < ghost->mirror_proc_offsets[n_proc]; ++i) {
    if (received_mirror_flags[i] == -1) {
      P4EST_ASSERT (virtual_ghost->mirror_proc_virtuals[i] == 0);
    }
    else {
      /* include a small test for virtual_gflags. At most all ghosts can embed
       * virtual quadrants, such that we have an upper bound here */
      P4EST_ASSERT ((virtual_ghost->mirror_proc_virtuals[i] == 1) &&
                    (0 <= received_mirror_flags[i]
                     && received_mirror_flags[i] < ghost->ghosts.elem_count));
    }
  }

  check_flat_exchange (p4est, ghost, mesh, virtual_quads, virtual_ghost);
  check_level_exchange (p4est, ghost, mesh, virtual_quads, virtual_ghost);

  /* cleanup */
  P4EST_FREE (received_mirror_flags);
  p4est_virtual_ghost_destroy (virtual_ghost);
  p4est_virtual_destroy (virtual_quads);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  return 0;
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
  MPI_Comm_set_errhandler(mpicomm, MPI_ERRORS_RETURN);
  
  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  check_virtual_ghost (mpicomm);

  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
