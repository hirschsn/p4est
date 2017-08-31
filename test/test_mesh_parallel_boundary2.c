/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
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

#include <inttypes.h>
#include <unistd.h>
#ifndef P4_TO_P8
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#else /* !P4_TO_P8 */
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#endif /* !P4_TO_P8 */

/** Check that indices in parallel boundary and mirrors are consistent.
 *
 * \param[in] p4est     The forest.
 * \param[in] ghost     Ghost layer.
 * \param[in] mesh      Neighbor information about forest and ghost layer.
 *                      Must contain information about parallel boundary.
 */
int
check_consistency_with_mirrors (p4est_t * p4est, p4est_ghost_t * ghost,
                                p4est_mesh_t * mesh)
{
  size_t              i;
  int                 count_mirrors = 0;
  int                 qid;
  p4est_quadrant_t   *q, *m;

  P4EST_ASSERT (mesh->parallel_boundary != 0);

  for (i = 0; i < ghost->mirrors.elem_count; ++i) {
    qid = mesh->mirror_qid[i];
    P4EST_ASSERT (0 <= qid && qid < mesh->local_num_quadrants);
    P4EST_ASSERT ((p4est_locidx_t) i == mesh->parallel_boundary[qid]);
  }
  for (i = 0; i < (size_t) mesh->local_num_quadrants; ++i) {
    if (-1 < mesh->parallel_boundary[i]) {
      q = p4est_mesh_get_quadrant (p4est, mesh, i);
      m =
        (p4est_quadrant_t *) sc_array_index (&ghost->mirrors, count_mirrors);
      // check for equality: position and level has to match
      P4EST_ASSERT (p4est_quadrant_is_equal (q, m));
      ++count_mirrors;
    }
  }
  return 0;
}

/** Check that all mirrors are detected correctly, i.e. if a quadrant has
 * neighbors in the ghost layer, parallel boundary must be set.  If all
 * neighbors are owned by the same rank, the flag must not be set.
 *
 * \param[in] p4est     The forest.
 * \param[in] ghost     Ghost layer.
 * \param[in] mesh      Neighbor information about forest and ghost layer.
 *                      Must contain information about parallel boundary.
 */
int
check_parallel_boundary_flag_is_valid (p4est_t * p4est, p4est_ghost_t * ghost,
                                       p4est_mesh_t * mesh)
{
  int                 quad, norm_quad;
  int                 i, imax, j;
  int                 lq = mesh->local_num_quadrants;
  int                 gq = mesh->ghost_num_quadrants;
  sc_array_t         *neighboring_qids;
  int                 neighbor_qid;
  int                 n_ghost_neighbors;

  P4EST_ASSERT (mesh->parallel_boundary != 0);

  switch (ghost->btype) {
  case P4EST_CONNECT_FACE:
    imax = P4EST_FACES;
    break;
#ifdef P4_TO_P8
  case P8EST_CONNECT_EDGE:
    imax = P4EST_FACES + P8EST_EDGES;
    break;
#endif /* P4_TO_P8 */
  case P4EST_CONNECT_FULL:
#ifdef P4_TO_P8
    imax = P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN;
#else /* P4_TO_P8 */
    imax = P4EST_FACES + P4EST_CHILDREN;
#endif /* P4_TO_P8 */
    break;
  default:
    SC_ABORT_NOT_REACHED ();
  }

  /** allocate containers */
  neighboring_qids = sc_array_new (sizeof (int));

  for (quad = 0; quad < p4est->global_num_quadrants; ++quad) {
    sc_MPI_Barrier (p4est->mpicomm);
    /** loop over all quads, verify only on the processor owning the respective
        quadrant. */
    if (p4est->global_first_quadrant[p4est->mpirank] <= quad &&
        quad < p4est->global_first_quadrant[p4est->mpirank + 1]) {
      /** norm global quad index to local index */
      norm_quad = quad - p4est->global_first_quadrant[p4est->mpirank];
      /**(quadrant must be local by design) */
      P4EST_ASSERT (0 <= norm_quad && norm_quad < lq);
      n_ghost_neighbors = 0;

      for (i = 0; i < imax; ++i) {
        /** empty containers */
        sc_array_truncate (neighboring_qids);

        /** search neighbors */
        p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad, i, NULL,
                                  NULL, neighboring_qids);

        /** inspect obtained neighbors */
        for (j = 0; j < (int) neighboring_qids->elem_count; ++j) {
          neighbor_qid = *(int *) sc_array_index (neighboring_qids, j);

          P4EST_ASSERT (0 <= neighbor_qid && neighbor_qid < (lq + gq));
          if (lq <= neighbor_qid && neighbor_qid < (lq + gq)) {
            ++n_ghost_neighbors;
          }
        }
      }
      if (n_ghost_neighbors) {
        P4EST_ASSERT (-1 < mesh->parallel_boundary[norm_quad]
                      && mesh->parallel_boundary[norm_quad] <
                      (p4est_locidx_t) ghost->mirrors.elem_count);
      }
      else {
        P4EST_ASSERT (-1 == mesh->parallel_boundary[norm_quad]);
      }
    }
  }
  /** de-allocate containers */
  sc_array_destroy (neighboring_qids);

  return 0;
}

/* Function for testing p4est_mesh for a single tree scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \param [in] mpicomm   MPI communicator
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_one_tree (p4est_t * p4est, p4est_connectivity_t * conn,
                    int8_t periodic, sc_MPI_Comm mpicomm)
{
  int                 minLevel = 3;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  /* create connectivity */
#ifndef P4_TO_P8
  conn = periodic == 1 ? p4est_connectivity_new_periodic ()
    : p4est_connectivity_new_unitsquare ();
#else /* !P4_TO_P8 */
  conn = periodic == 1 ? p8est_connectivity_new_periodic ()
    : p8est_connectivity_new_unitcube ();
#endif /* !P4_TO_P8 */

  /* setup p4est */
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);
  p4est_balance (p4est, btype, NULL);

  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  check_consistency_with_mirrors (p4est, ghost, mesh);
  check_parallel_boundary_flag_is_valid (p4est, ghost, mesh);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;

  return 0;
}

/** Function for testing p4est_mesh for all existing two_trees scenarios
 * \param [in] p4est     The forest, NULL.
 * \param [in] conn      The connectivity structure, NULL.
 * \param [in] mpicomm   MPI communicator
 */
int
test_mesh_two_trees (p4est_t * p4est, p4est_connectivity_t * conn,
                     sc_MPI_Comm mpicomm)
{
  int                 conn_face_tree1, conn_face_tree2, orientation;
  int                 minLevel = 3;
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);
  for (conn_face_tree1 = 0; conn_face_tree1 < P4EST_FACES; ++conn_face_tree1) {
    for (conn_face_tree2 = 0; conn_face_tree2 < P4EST_FACES;
         ++conn_face_tree2) {
      for (orientation = 0; orientation < P4EST_HALF; ++orientation) {
        P4EST_VERBOSEF ("Check if parallel boundary is detected properly for"
                        " two trees, left face: %i, right face: %i,"
                        " orientation: %i\n", conn_face_tree1,
                        conn_face_tree2, orientation);
        /* create connectivity */
        conn =
          p4est_connectivity_new_twotrees (conn_face_tree1, conn_face_tree2,
                                           orientation);
        /* setup p4est */
        p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);
        p4est_balance (p4est, btype, NULL);

        ghost = p4est_ghost_new (p4est, btype);
        mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
        check_consistency_with_mirrors (p4est, ghost, mesh);
        check_parallel_boundary_flag_is_valid (p4est, ghost, mesh);

        /* cleanup */
        p4est_ghost_destroy (ghost);
        p4est_mesh_destroy (mesh);
        p4est_destroy (p4est);
        p4est_connectivity_destroy (conn);

        conn = 0;
        p4est = 0;
      }
    }
  }

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \param [in] mpicomm   MPI communicator
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_brick (p4est_t * p4est, p4est_connectivity_t * conn,
                                int8_t periodic, sc_MPI_Comm mpicomm)
{
  int                 minLevel = 3;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  P4EST_VERBOSEF
    ("Check if parallel boundary is detected properly for brick of"
     " trees, periodic: %i\n", periodic);

  /* create connectivity */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (2, 1, periodic, periodic);
#else /* !P4_TO_P8 */
  conn = p8est_connectivity_new_brick (2, 1, 1, periodic, periodic, periodic);
#endif /* !P4_TO_P8 */

  /* setup p4est */
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);
  p4est_balance (p4est, btype, NULL);

  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  check_consistency_with_mirrors (p4est, ghost, mesh);
  check_parallel_boundary_flag_is_valid (p4est, ghost, mesh);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a non-brick scenario
 *
 * \param [in] p4est     The forest.
 * \param [in] conn      The connectivity structure
 * \param [in] periodic  Flag for checking if we have periodic boundaries
 * \returns 0 for success, -1 for failure
 */
int
test_mesh_multiple_trees_nonbrick (p4est_t * p4est,
                                   p4est_connectivity_t * conn,
                                   int8_t periodic, sc_MPI_Comm mpicomm)
{
  return 0;
}

int
test_saved_tree (p4est_t * p4est, p4est_connectivity_t * conn,
                 sc_MPI_Comm mpicomm)
{
#ifdef P4_TO_P8
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;

  p4est = p4est_load ("broken_forest.p4est", mpicomm, 0, 0, 0, &conn);
  P4EST_ASSERT (p4est_is_balanced (p4est, btype));

  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  check_consistency_with_mirrors (p4est, ghost, mesh);
  check_parallel_boundary_flag_is_valid (p4est, ghost, mesh);

  /* cleanup */
  p4est_ghost_destroy (ghost);
  p4est_mesh_destroy (mesh);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;
#endif /* P4_TO_P8 */

  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  int8_t              periodic_boundaries;
  int                 test_single, test_two_trees;
  int                 test_multi_brick, test_multi_non_brick;
  int                 test_saved_config;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;
  mpiret = sc_MPI_Comm_size (mpicomm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &mpirank);
  SC_CHECK_MPI (mpiret);
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);
  p4est = 0;
  conn = 0;

  test_single = 1;
  test_two_trees = 1;
  test_multi_brick = 1;
  test_multi_non_brick = 0;
  test_saved_config = 0;

  /* test both periodic and non-periodic boundaries */
  if (test_single) {
    /* test one tree */
    periodic_boundaries = 0;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);

    periodic_boundaries = 1;
    test_mesh_one_tree (p4est, conn, periodic_boundaries, mpicomm);
  }

  if (test_two_trees) {
    /* test two trees
     * (using all combinations of faces connected and all orientations) */
    test_mesh_two_trees (p4est, conn, mpicomm);
  }

  if (test_multi_brick) {
    /* test multiple trees; brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_brick (p4est, conn, periodic_boundaries,
                                    mpicomm);
  }

  if (test_multi_non_brick) {
    /* test multiple trees; non-brick */
    periodic_boundaries = 0;
    test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                       mpicomm);
    periodic_boundaries = 1;
    test_mesh_multiple_trees_nonbrick (p4est, conn, periodic_boundaries,
                                       mpicomm);
  }

  if (test_saved_config) {
    test_saved_tree (p4est, conn, mpicomm);
  }
  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
