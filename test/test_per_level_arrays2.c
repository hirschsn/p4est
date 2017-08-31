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
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_mesh.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#endif /* !P4_TO_P8 */

/** Check that all quadrants are found in a per-level traversal.
 *
 * \param[in] p4est     The forest.
 * \param[in] ghost     Ghost layer.
 * \param[in] mesh      Neighbor information about forest and ghost layer.
 *                      Must contain level information
 */
int
check_consistency_of_level_array (p4est_t * p4est, p4est_ghost_t * ghost,
                                  p4est_mesh_t * mesh)
{
  p4est_locidx_t      lq, gq;
  int                 level, i;
  p4est_locidx_t      qid;
  p4est_locidx_t     *found_local, *found_ghost;
  p4est_quadrant_t   *q;

  P4EST_ASSERT (mesh->quad_level != 0 && mesh->ghost_level != 0);

  lq = mesh->local_num_quadrants;
  gq = mesh->ghost_num_quadrants;

  found_local = P4EST_ALLOC_ZERO (p4est_locidx_t, lq);
  found_ghost = P4EST_ALLOC_ZERO (p4est_locidx_t, gq);

  for (level = 0; level < P4EST_QMAXLEVEL + 1; ++level) {
    for (i = 0; i < (mesh->quad_level + level)->elem_count; ++i) {
      qid = *(p4est_locidx_t *) sc_array_index (mesh->quad_level + level, i);
      P4EST_ASSERT (found_local[qid] == 0);
      found_local[qid] = 1;
      q = p4est_mesh_get_quadrant (p4est, mesh, qid);
      P4EST_ASSERT (level == q->level);
    }
    for (i = 0; i < (mesh->ghost_level + level)->elem_count; ++i) {
      qid = *(p4est_locidx_t *) sc_array_index (mesh->ghost_level + level, i);
      P4EST_ASSERT (found_ghost[qid] == 0);
      found_ghost[qid] = 1;
      q = p4est_quadrant_array_index (&ghost->ghosts, qid);
      P4EST_ASSERT (level == q->level);
    }
  }

  for (qid = 0; qid < lq; ++qid) {
    P4EST_ASSERT (1 == found_local[qid]);
  }
  for (qid = 0; qid < gq; ++qid) {
    P4EST_ASSERT (1 == found_ghost[qid]);
  }

  P4EST_FREE (found_local);
  P4EST_FREE (found_ghost);

  return 0;
}

/* Function for testing p4est-mesh for a single tree scenario
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
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
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
  check_consistency_of_level_array (p4est, ghost, mesh);

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
        check_consistency_of_level_array (p4est, ghost, mesh);

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
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

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
  check_consistency_of_level_array (p4est, ghost, mesh);

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

  /* create mesh */
  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  check_consistency_of_level_array (p4est, ghost, mesh);

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

  int                 test_single, test_two_trees;
  int                 test_multi_brick, test_multi_non_brick;
  int                 test_saved_config;
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
