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
#include <p4est_virtual.h>
#else /* !P4_TO_P8 */
#include <p8est_connectivity.h>
#include <p8est_virtual.h>
#endif /* !P4_TO_P8 */

/** Investigate p4est_face_virtual_neighbors_inside matrix.  Count number of
 * occurances of each entry and verify that neighborhood within a quadrant is
 * detected correctly.
 */
static int
check_face_neighbor_matrix ()
{
  int                 space_to_alloc;
  int                *count_occurances;
  int                 val, i, j, k;
  int                 l_internal, u_internal;
  int                 l_face, u_face;
  int                 internal_face_neighbors[P4EST_DIM];
  int                 must_be_internal;

#ifndef P4_TO_P8
  space_to_alloc = 12;
#else /* !P4_TO_P8 */
  space_to_alloc = 32;
#endif /* !P4_TO_P8 */
  l_internal = 0;
  u_internal = l_face = P4EST_CHILDREN;
  u_face = space_to_alloc;

  count_occurances = P4EST_ALLOC_ZERO (int, space_to_alloc);

  for (i = 0; i < P4EST_CHILDREN; ++i) {
    for (j = 0; j < P4EST_DIM; ++j) {
      internal_face_neighbors[j] = p4est_face_dual[p4est_corner_faces[i][j]];
    }
    for (j = 0; j < P4EST_FACES; ++j) {
      must_be_internal = 0;
      for (k = 0; k < P4EST_DIM; ++k) {
        if (j == internal_face_neighbors[k]) {
          must_be_internal = 1;
          break;
        }
      }
      val = p4est_face_virtual_neighbors_inside[i][j];
      if (must_be_internal) {
        P4EST_ASSERT (0 <= val && val < P4EST_CHILDREN);
      }
      else {
        P4EST_ASSERT (P4EST_CHILDREN <= val && val < space_to_alloc);
      }
      P4EST_ASSERT (0 <= val && val < space_to_alloc);
      ++count_occurances[val];
    }
  }

  val = 0;
  for (i = l_internal; i < u_internal; ++i) {
    P4EST_ASSERT (count_occurances[i] == P4EST_DIM);
    ++val;
  }
  for (i = l_face; i < u_face; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }

  P4EST_ASSERT (val == space_to_alloc);
  P4EST_FREE (count_occurances);

  return 0;
}

#ifdef P4_TO_P8
/** Investigate p8est_edge_virtual_neighbors_inside matrix.  Count number of
 * occurances of each entry and verify that neighborhood within a quadrant is
 * detected correctly.
 */
static int
check_edge_neighbor_matrix ()
{
  int                 space_to_alloc;
  int                *count_occurances;
  int                 val, i, j, k;
  int                 must_be_internal;
  int                 must_be_edge_query;
  int                 l_internal, u_internal;
  int                 l_face, u_face;
  int                 l_edge, u_edge;
  int                 internal_edge_neighbors[P4EST_DIM];

  space_to_alloc = 56;
  l_internal = 0;
  u_internal = l_face = P4EST_CHILDREN;
  u_face = l_edge = 32;
  u_edge = space_to_alloc;

  count_occurances = P4EST_ALLOC_ZERO (int, space_to_alloc);

  for (i = 0; i < P4EST_CHILDREN; ++i) {
    for (j = 0; j < P4EST_DIM; ++j) {
      internal_edge_neighbors[j] = p8est_corner_edges[i][j] ^ 3;
    }
    for (j = 0; j < P8EST_EDGES; ++j) {
      must_be_internal = 0;
      must_be_edge_query = 0;
      val = p8est_edge_virtual_neighbors_inside[i][j];
      for (k = 0; k < P4EST_DIM; ++k) {
        if (j == internal_edge_neighbors[k]) {
          must_be_internal = 1;
          break;
        }
        if (j == p8est_corner_edges[i][k]) {
          must_be_edge_query = 1;
          break;
        }
      }
      P4EST_ASSERT (!(must_be_internal && must_be_edge_query));
      if (must_be_internal) {
        P4EST_ASSERT (0 <= val && val < P4EST_CHILDREN);
      }
      else if (must_be_edge_query) {
        P4EST_ASSERT (l_edge <= val && val < u_edge);
      }
      else {
        P4EST_ASSERT (l_face <= val && val < u_face);
      }
      ++count_occurances[val];
    }
  }

  val = 0;
  for (i = l_internal; i < u_internal; ++i) {
    P4EST_ASSERT (P4EST_DIM == count_occurances[i]);
    ++val;
  }
  for (i = l_face; i < u_face; ++i) {
    P4EST_ASSERT (2 == count_occurances[i]);
    ++val;
  }
  for (i = l_edge; i < u_edge; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }

  P4EST_ASSERT (val == space_to_alloc);
  P4EST_FREE (count_occurances);

  return 0;
}
#endif /* P4_TO_P8 */

/** Investigate p4est_corner_virtual_neighbors_inside matrix.  Count number of
 * occurances of each entry and verify that neighborhood within a quadrant is
 * detected correctly.
 */
static int
check_corner_neighbor_matrix ()
{
  int                 space_to_alloc;
  int                *count_occurances;
  int                 val, i, j;
  int                 l_internal, u_internal;
  int                 l_face, u_face;
#ifdef P4_TO_P8
  int                 k, l;
  int                 l_edge, u_edge;
  int                 must_be_edge_query;
#endif /* P4_TO_P8 */
  int                 l_corner, u_corner;
  int                 internal_corner_neighbor;

#ifndef P4_TO_P8
  space_to_alloc = 16;
#else /* !P4_TO_P8 */
  space_to_alloc = 64;
#endif /* !P4_TO_P8 */
  l_internal = 0;
  u_internal = l_face = P4EST_CHILDREN;
#ifndef P4_TO_P8
  u_face = l_corner = 12;
#else /* !P4_TO_P8 */
  u_face = l_edge = 32;
  u_edge = l_corner = 56;
#endif /* !P4_TO_P8 */
  u_corner = space_to_alloc;
  count_occurances = P4EST_ALLOC_ZERO (int, space_to_alloc);

  for (i = 0; i < P4EST_CHILDREN; ++i) {
    /* types of queries:
     * a) internal:  an internal virtual quadrant needs to be queried for the
     *               diagonally opposite corner of the corner the current
     *               virtual quadrant is associated with.
     * b) corner:    a corner neighbor needs to be queried for the corner that
     *               the current virtual quadrant is associated with.
     * c) edge (3D): an edge neighbor needs to be queried if the corner index is
     *               among the three corner indices that are opposite across one
     *               of the edges that are adjacent to the corner the current
     *               virtual quadrant is associated with.
     * d) face:      Anything else.
     */
    internal_corner_neighbor = i ^ (P4EST_CHILDREN - 1);
    for (j = 0; j < P4EST_CHILDREN; ++j) {
#ifdef P4_TO_P8
      must_be_edge_query = -1;
      for (k = 0; k < P4EST_DIM; ++k) {
        if (must_be_edge_query == j) {
          break;
        }
        for (l = 0; l < 2; ++l) {
          must_be_edge_query =
            p8est_edge_corners[p8est_corner_edges[i][k]][l];
          if (must_be_edge_query == i) {
            must_be_edge_query = -1;
            continue;
          }
          else {
            break;
          }
        }
      }
      if (must_be_edge_query != j) {
        must_be_edge_query = 0;
      }
      else {
        must_be_edge_query = 1;
      }
#endif /* P4_TO_P8 */
      val = p4est_corner_virtual_neighbors_inside[i][j];
      if (j == internal_corner_neighbor) {
        P4EST_ASSERT (0 <= val && val < P4EST_CHILDREN);
      }
#ifdef P4_TO_P8
      else if (must_be_edge_query) {
        P4EST_ASSERT (l_edge <= val && val < u_edge);
      }
#endif /* P4_TO_P8 */
      else if (j == i) {
        P4EST_ASSERT (l_corner <= val && val < u_corner);
      }
      else {
        P4EST_ASSERT (l_face <= val && val < u_face);
      }
      ++count_occurances[val];
    }
  }

  val = 0;
  for (i = l_internal; i < u_internal; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }
  for (i = l_face; i < u_face; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }
#ifdef P4_TO_P8
  for (i = l_edge; i < u_edge; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }
#endif /* P4_TO_P8 */
  for (i = l_corner; i < u_corner; ++i) {
    P4EST_ASSERT (1 == count_occurances[i]);
    ++val;
  }

  P4EST_ASSERT (val == space_to_alloc);
  P4EST_FREE (count_occurances);

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

  /* initialize libsc and p4est */
  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  check_face_neighbor_matrix ();
#ifdef P4_TO_P8
  check_edge_neighbor_matrix ();
#endif /* P4_TO_P8 */
  check_corner_neighbor_matrix ();

  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
