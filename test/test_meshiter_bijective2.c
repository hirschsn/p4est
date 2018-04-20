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
#include <p4est_meshiter.h>
#include <p4est_virtual.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_meshiter.h>
#include <p8est_virtual.h>
#endif /* !P4_TO_P8 */

int
check_bijectivity (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads)
{
  int                 level;
  int                 i, imax, iinv;
  int                 success;
  int                 status_main, status_help;
  p4est_meshiter_t   *main_iterator, *help_iterator;

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

  for (level = 0; level < P4EST_MAXLEVEL; ++level) {
    status_main = 0;
    main_iterator =
      p4est_meshiter_new_ext (p4est, ghost, mesh, virtual_quads, level,
                              ghost->btype, P4EST_TRAVERSE_LOCAL,
                              P4EST_TRAVERSE_REALVIRTUAL,
                              P4EST_TRAVERSE_PARBOUNDINNER);

    while (status_main != P4EST_MESHITER_DONE) {
      status_main = p4est_meshiter_next (main_iterator);
      if (status_main != P4EST_MESHITER_DONE) {
        for (i = 0; i < imax; ++i) {
          success = 0;

          p4est_meshiter_set_neighbor_quad_info (main_iterator, i);

          if ((main_iterator->neighbor_qid == -1)
              || (main_iterator->neighbor_is_ghost)) {
            success = 1;
          }
          else {
            iinv = main_iterator->neighbor_entity_index;
            help_iterator =
              p4est_meshiter_new_ext (p4est, ghost, mesh, virtual_quads,
                                      level, ghost->btype,
                                      P4EST_TRAVERSE_LOCAL,
                                      P4EST_TRAVERSE_REALVIRTUAL,
                                      P4EST_TRAVERSE_PARBOUNDINNER);

            while (status_help != P4EST_MESHITER_DONE) {
              status_help = p4est_meshiter_next (help_iterator);
              if ((help_iterator->current_qid == main_iterator->neighbor_qid)
                  && (help_iterator->current_vid ==
                      main_iterator->neighbor_vid)) {
                p4est_meshiter_set_neighbor_quad_info (help_iterator, iinv);
                P4EST_ASSERT (help_iterator->neighbor_qid ==
                              main_iterator->current_qid);
                P4EST_ASSERT (help_iterator->neighbor_vid ==
                              main_iterator->current_vid);
                P4EST_ASSERT (p4est_meshiter_get_current_storage_id
                              (main_iterator) ==
                              p4est_meshiter_get_neighbor_storage_id
                              (help_iterator));
                P4EST_ASSERT (p4est_meshiter_get_current_storage_id
                              (help_iterator) ==
                              p4est_meshiter_get_neighbor_storage_id
                              (main_iterator));
                success = 1;
                break;
              }
            }
            p4est_meshiter_destroy (help_iterator);
          }
          P4EST_ASSERT (1 == success);
        }
      }
    }

    p4est_meshiter_destroy (main_iterator);
  }
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
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  int                 minLevel = 3;

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
  virtual_quads = p4est_virtual_new_ext (p4est, ghost, mesh, btype, 1);

  check_bijectivity (p4est, ghost, mesh, virtual_quads);

  /* cleanup */
  p4est_virtual_destroy (virtual_quads);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
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
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  int                 minLevel = 3;
  int                 conn_face_tree1, conn_face_tree2, orientation;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  for (conn_face_tree1 = 0; conn_face_tree1 < P4EST_FACES; ++conn_face_tree1) {
    for (conn_face_tree2 = 0; conn_face_tree2 < P4EST_FACES;
         ++conn_face_tree2) {
      for (orientation = 0; orientation < P4EST_HALF; ++orientation) {
        /* create connectivity */
        conn =
          p4est_connectivity_new_twotrees (conn_face_tree1, conn_face_tree2,
                                           orientation);
        /* setup p4est */
        p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);
        p4est_balance (p4est, btype, NULL);

        /* create mesh */
        ghost = p4est_ghost_new (p4est, btype);
        mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
        virtual_quads = p4est_virtual_new_ext (p4est, ghost, mesh, btype, 1);

        check_bijectivity (p4est, ghost, mesh, virtual_quads);

        /* cleanup */
        p4est_virtual_destroy (virtual_quads);
        p4est_mesh_destroy (mesh);
        p4est_ghost_destroy (ghost);
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
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  int                 minLevel = 3;

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

  /* create mesh */
  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  virtual_quads = p4est_virtual_new_ext (p4est, ghost, mesh, btype, 1);

  check_bijectivity (p4est, ghost, mesh, virtual_quads);

  /* cleanup */
  p4est_virtual_destroy (virtual_quads);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
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
  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;
}
