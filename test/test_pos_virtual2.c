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
#include <p4est_virtual.h>
#else /* !P4_TO_P8 */
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>
#include <p8est_virtual.h>
#endif /* !P4_TO_P8 */

/** Set constants needed to decode an encoding obtained from \ref
 * p4est_mesh_get_neighbors
 * \param [in]      direction             Encoded direction for neighbor search
 *                                        analog to direction used in
 *                                        \ref p4est_mesh_get_neighbors
 * \param [in]      btype                 Connectivity type (used for asserting
 *                                        that direction is valid).
 * \param     [out] l_same_size           Minimum value of an encoding
 *                                        indicating a same-sized neighbor.
 * \param     [out] u_same_size           Maximum value of an encoding
 *                                        indicating a same-sized neighbor.
 * \param     [out] l_double_size         Minimum value of an encoding
 *                                        indicating a double-sized neighbor.
 * \param     [out] u_double_size         Maximum value of an encoding
 *                                        indicating a double-sized neighbor.
 * \param     [out] l_half_size           Minimum value of an encoding
 *                                        indicating a half-sized neighbor.
 * \param     [out] u_half_size           Maximum value of an encoding
 *                                        indicating a half-sized neighbor.
 * \param     [out] n_neighbor_entities   Number of different entities of the
 *                                        respective neighbor entity, i.e.
 *                                        number of faces, (edges), or corners.
 * \param     [out] n_hanging_quads       Number of hanging quads across face
 *                                        (or edge).
 */
static int
set_limits (int direction, p4est_connect_type_t btype,
            int *l_same_size, int *u_same_size, int *l_double_size,
            int *u_double_size, int *l_half_size, int *u_half_size)
{
#ifndef P4_TO_P8
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = 2 * P4EST_FACES;
    *l_double_size = *u_same_size;
    *u_double_size = 24;
    *l_half_size = -8;
    *u_half_size = *l_same_size;
  }
  else if (P4EST_FACES <= direction &&
           direction < (P4EST_FACES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_half_size = *u_half_size = *l_same_size = 0;
    *u_same_size = *l_double_size = *u_double_size = P4EST_CHILDREN;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#else /* !P4_TO_P8 */
  if (0 <= direction && direction < P4EST_FACES) {
    *l_same_size = 0;
    *u_same_size = 4 * P4EST_FACES;
    *l_double_size = *u_same_size;
    *u_double_size = 120;
    *l_half_size = -24;
    *u_half_size = *l_same_size;
  }
  else if (P4EST_FACES <= direction &&
           direction < (P4EST_FACES + P8EST_EDGES)) {
    P4EST_ASSERT (btype >= P8EST_CONNECT_EDGE);
    *l_same_size = 0;
    *u_same_size = 2 * P8EST_EDGES;
    *l_double_size = *u_same_size;
    *u_double_size = 72;
    *l_half_size = -24;
    *u_half_size = *l_same_size;
  }
  else if ((P4EST_FACES + P8EST_EDGES) <= direction &&
           direction < (P4EST_FACES + P8EST_EDGES + P4EST_CHILDREN)) {
    P4EST_ASSERT (btype == P4EST_CONNECT_CORNER);
    *l_half_size = *u_half_size = *l_same_size = 0;
    *u_same_size = *l_double_size = *u_double_size = P4EST_CHILDREN;
  }
  else {
    SC_ABORT_NOT_REACHED ();
  }
#endif /* !P4_TO_P8 */

  return 0;
}

void
check_pos_virtual (p4est_t * p4est, p4est_ghost_t * ghost,
                   p4est_mesh_t * mesh, p4est_virtual_t * virtual_quads)
{
  int                 quad, norm_quad;
  int                 i, imax, j;
  p4est_locidx_t      lq = mesh->local_num_quadrants;
  p4est_locidx_t      gq = mesh->ghost_num_quadrants;
  p4est_quadrant_t   *current_quad;
  int                 l_same_size, u_same_size;
  int                 l_double_size, u_double_size;
  int                 l_half_size, u_half_size;
  int                 needs_virtuals;
  sc_array_t         *neighboring_quads;
  sc_array_t         *neighboring_encs;
  sc_array_t         *neighboring_qids;
  p4est_quadrant_t   *neighbor_quad;
  int                 neighbor_qid, neighbor_enc;

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

  /** create containers */
  neighboring_quads = sc_array_new (sizeof (p4est_quadrant_t *));
  neighboring_encs = sc_array_new (sizeof (int));
  neighboring_qids = sc_array_new (sizeof (int));

  for (quad = 0; quad < p4est->global_num_quadrants; ++quad) {
    sc_MPI_Barrier (p4est->mpicomm);
    /** loop over all quads, verify only on the processor owning the respective
        quadrant. */
    if (p4est->global_first_quadrant[p4est->mpirank] <= quad &&
        quad < p4est->global_first_quadrant[p4est->mpirank + 1]) {
      /** norm global quad index to local index */
      norm_quad = quad - p4est->global_first_quadrant[p4est->mpirank];
      P4EST_ASSERT (0 <= norm_quad && norm_quad < lq);
      current_quad = p4est_mesh_get_quadrant (p4est, mesh, norm_quad);

      needs_virtuals = 0;

      for (i = 0; i < imax; ++i) {
        /** empty containers */
        sc_array_truncate (neighboring_quads);
        sc_array_truncate (neighboring_encs);
        sc_array_truncate (neighboring_qids);

        /** set constants for decoding */
        set_limits (i, ghost->btype, &l_same_size, &u_same_size,
                    &l_double_size, &u_double_size, &l_half_size,
                    &u_half_size);

        /** search neighbors */
        p4est_mesh_get_neighbors (p4est, ghost, mesh, norm_quad, i,
                                  neighboring_quads, neighboring_encs,
                                  neighboring_qids);
        P4EST_ASSERT (neighboring_quads->elem_count ==
                      neighboring_qids->elem_count);
        P4EST_ASSERT (neighboring_encs->elem_count ==
                      neighboring_qids->elem_count);

        /** inspect obtained neighbors */
        for (j = 0; j < (int) neighboring_qids->elem_count; ++j) {
          neighbor_quad =
            *(p4est_quadrant_t **) sc_array_index (neighboring_quads, j);
          neighbor_qid = *(int *) sc_array_index (neighboring_qids, j);
          neighbor_enc = *(int *) sc_array_index (neighboring_encs, j);

          if ((l_double_size <= neighbor_enc && neighbor_enc < u_double_size)
              || (current_quad->level == neighbor_quad->level + 1)) {
            /* check if ghost quadrant has virtual quadrants */
            if (lq <= neighbor_qid && neighbor_qid < (lq + gq)) {
              P4EST_ASSERT (virtual_quads->
                            virtual_gflags[(neighbor_qid - lq)] != -1);
            }
          }
          else if ((l_half_size <= neighbor_enc && neighbor_enc < u_half_size)
                   || (current_quad->level + 1 == neighbor_quad->level)) {
            ++needs_virtuals;
          }
        }
      }
      /* if we need virtuals virtuals need to be build and if we do not need
         them they must not be present. */
      if (needs_virtuals) {
        P4EST_ASSERT (-1 != virtual_quads->virtual_qflags[norm_quad]);
      }
      else {
        P4EST_ASSERT (-1 == virtual_quads->virtual_qflags[norm_quad]);
      }
    }
  }
  /** de-allocate containers */
  sc_array_destroy (neighboring_quads);
  sc_array_destroy (neighboring_encs);
  sc_array_destroy (neighboring_qids);
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
  p4est_virtual_t    *virtual_quads;

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

  /* create mesh */

  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  virtual_quads = p4est_virtual_new (p4est, ghost, mesh, btype);

  /* check mesh */
  check_pos_virtual (p4est, ghost, mesh, virtual_quads);

  /* cleanup */
  p4est_virtual_destroy (virtual_quads);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

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
  p4est_virtual_t    *virtual_quads;

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
        virtual_quads = p4est_virtual_new (p4est, ghost, mesh, btype);

        /* check mesh */
        check_pos_virtual (p4est, ghost, mesh, virtual_quads);

        /* cleanup */
        p4est_virtual_destroy (virtual_quads);
        p4est_mesh_destroy (mesh);
        p4est_ghost_destroy (ghost);
        p4est_destroy (p4est);
        p4est_connectivity_destroy (conn);

        conn = 0;
        p4est = 0;
        P4EST_ASSERT (p4est == NULL);
        P4EST_ASSERT (conn == NULL);
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
  p4est_virtual_t    *virtual_quads;

  /* ensure that we have null pointers at beginning and end of function */
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  /* create connectivity */
#ifndef P4_TO_P8
  conn = p4est_connectivity_new_brick (2, 2, periodic, periodic);
#else /* !P4_TO_P8 */
  conn = p8est_connectivity_new_brick (2, 2, 2, periodic, periodic, periodic);
#endif /* !P4_TO_P8 */

  /* setup p4est */
  p4est = p4est_new_ext (mpicomm, conn, 0, minLevel, 0, 0, NULL, NULL);
  p4est_balance (p4est, btype, NULL);

  /* create mesh */
  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 1, 1, btype);
  virtual_quads = p4est_virtual_new (p4est, ghost, mesh, btype);

  /* check mesh */
  check_pos_virtual (p4est, ghost, mesh, virtual_quads);

  /* cleanup */
  p4est_virtual_destroy (virtual_quads);
  p4est_mesh_destroy (mesh);
  p4est_ghost_destroy (ghost);
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  conn = 0;
  p4est = 0;
  P4EST_ASSERT (p4est == NULL);
  P4EST_ASSERT (conn == NULL);

  return 0;
}

/* Function for testing p4est-mesh for multiple trees in a non-brick
 * scenario. Not used yet.
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
  test_single = 1;
  test_two_trees = 1;
  test_multi_brick = 1;
  test_multi_non_brick = 0;

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
