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

#include <math.h>
#include <stdlib.h>
#include <time.h>

#ifndef P4_TO_P8
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_search.h>
#else /* P4_TO_P8 */
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <p8est_search.h>
#endif /* P4_TO_P8 */

/** Generate a random double between 0 and 1
 */
double
random_double ()
{
  double              r = (double) rand () / (double) RAND_MAX;
  return r;
}

/** Calculate a 2D Morton Index by bit-interleaving
 *
 * \param[in] x, y    Spatial coordinate
 */
int64_t
morton_idx_2d (int x, int y)
{
  int64_t             idx = 0;
  int64_t             pos = 1;

  for (int i = 0; i < 32; ++i) {
    if ((x & 1))
      idx += pos;
    x >>= 1;
    pos <<= 1;
    if ((y & 1))
      idx += pos;
    y >>= 1;
    pos <<= 1;
  }

  return idx;
}

/** Calculate a 3D Morton Index by bit-interleaving
 *
 * \param[in] x, y, z  Spatial coordinate
 */
int64_t
morton_idx_3d (int x, int y, int z)
{
  int64_t             idx = 0;
  int64_t             pos = 1;

  for (int i = 0; i < 21; ++i) {
    if ((x & 1))
      idx += pos;
    x >>= 1;
    pos <<= 1;
    if ((y & 1))
      idx += pos;
    y >>= 1;
    pos <<= 1;
    if ((z & 1))
      idx += pos;
    z >>= 1;
    pos <<= 1;
  }

  return idx;
}

/** Search a quadrant that contains a spatial coordinate
 *
 * \param [in]    p4est     The forest.
 */
int
test_pos_to_quad (p4est_t * p4est)
{
  p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, 0);
  p4est_quadrant_t    q;
  const int           forest_level = tree->maxlevel;

  /** number of quadrants per spatial dimension given a regular grid */
  const int           nq = 1 << forest_level;

  /** number of positions to draw */
  const int           n_samples = 5000;

  /** init random number generator */
  time_t              t;
  srand ((unsigned) time (&t));

  double              pos[P4EST_DIM];
  int                 qpos[P4EST_DIM];

  for (int i = 0; i < n_samples; ++i) {
    sc_MPI_Barrier (p4est->mpicomm);
    if (p4est->mpirank == 0) {
      for (int j = 0; j < P4EST_DIM; ++j) {
        pos[j] = random_double ();
        P4EST_ASSERT (0 <= pos[j] && pos[j] < 1);
      }
    }
    sc_MPI_Bcast (pos, P4EST_DIM, sc_MPI_DOUBLE, 0, p4est->mpicomm);
    for (int j = 0; j < P4EST_DIM; ++j) {
      qpos[j] = nq * pos[j];
    }
    /** generate Morton-index */
#ifdef P4_TO_P8
    const int           qid = morton_idx_3d (qpos[0], qpos[1], qpos[2]);
#else /* P4_TO_P8 */
    const int           qid = morton_idx_2d (qpos[0], qpos[1]);
#endif /* P4_TO_P8 */

    /** initialize dummy quadrant (level and spatial position) */
    P4EST_QUADRANT_INIT (&q);
    q.level = forest_level;
    p4est_quadrant_set_morton (&q, forest_level, qid);

    /** search first quadrant q* where Morton-Index(q*) <= Morton-Index(q) */
    int                 idx_lb =
      p4est_find_lower_bound_overlap (&tree->quadrants, &q,
                                      0.5 * tree->quadrants.elem_count);
    int                 idx_hb =
      p4est_find_higher_bound_overlap (&tree->quadrants, &q,
                                       0.5 * tree->quadrants.elem_count);

    p4est_quadrant_t   *q_found;
    if (p4est->global_first_quadrant[p4est->mpirank] <= idx_lb &&
        idx_lb < p4est->global_first_quadrant[p4est->mpirank + 1]) {
      /** Verification: quadrants q and q* must overlap */
      q_found = p4est_quadrant_array_index (&tree->quadrants, idx_lb);
      P4EST_ASSERT (p4est_quadrant_overlaps (&q, q_found));
    }

    if (p4est->global_first_quadrant[p4est->mpirank] <= idx_hb &&
        idx_hb < p4est->global_first_quadrant[p4est->mpirank + 1]) {
      /** Verification: quadrants q and q* must overlap */
      q_found = p4est_quadrant_array_index (&tree->quadrants, idx_hb);
      P4EST_ASSERT (p4est_quadrant_overlaps (&q, q_found));
    }

    if (idx_hb != idx_lb) {
#ifdef P4_TO_P8
      printf ("[p4est %i] idx_hb: %i idx_lb: %i pos %lf %lf %lf\n",
              p4est->mpirank, idx_hb, idx_lb, pos[0], pos[1], pos[2]);
#else
      printf ("[p4est %i] idx_hb: %i idx_lb: %i pos %lf %lf\n",
              p4est->mpirank, idx_hb, idx_lb, pos[0], pos[1]);
#endif
    }
    P4EST_ASSERT (idx_hb == idx_lb);
  }

  return 0;
}

/** Test upper bound functionality by constructing artificial quadrants
 *
 * \param[in]   p4est     The forest.
 */
int
test_upper_bound (p4est_t * p4est)
{
  p4est_tree_t       *tree = p4est_tree_array_index (p4est->trees, 0);
  int                 forest_level = tree->maxlevel;

  /** artificially refine quadrants compared to max level of forest */
  const int           quad_level = forest_level + 2;

  /** create 3 quadrants that are contained in the first quadrant (aligned with
   * lower left corner, somewhere inbetween, and aligned with the upper right
   * corner). Additionally, create another quadrant aligned with the lower left
   * corner of the second quadrant of the forest.*/
  p4est_quadrant_t    lower_left, inbetween, upper_right, next_ll;

  P4EST_QUADRANT_INIT (&lower_left);
  P4EST_QUADRANT_INIT (&inbetween);
  P4EST_QUADRANT_INIT (&upper_right);
  P4EST_QUADRANT_INIT (&next_ll);

  /** initialize level and position */
  lower_left.level = quad_level;
  inbetween.level = quad_level;
  upper_right.level = quad_level;
  next_ll.level = quad_level;

  /** Calculate index of upper right quadrant depending on level difference
   * between forest and quadrants */
  int                 upper_right_index =
    pow (P4EST_CHILDREN, quad_level - forest_level) - 1;
  p4est_quadrant_set_morton (&lower_left, quad_level, 0);
  p4est_quadrant_set_morton (&inbetween, quad_level, 0.5 * upper_right_index);
  p4est_quadrant_set_morton (&upper_right, quad_level, upper_right_index);
  p4est_quadrant_set_morton (&next_ll, quad_level, upper_right_index + 1);

  int                 ub_ll, ub_ib, ub_ur, ub_nll;
  ub_ll = p4est_find_higher_bound_overlap (&tree->quadrants, &lower_left,
                                           0.5 * tree->quadrants.elem_count);
  ub_ib = p4est_find_higher_bound_overlap (&tree->quadrants, &inbetween,
                                           0.5 * tree->quadrants.elem_count);
  ub_ur = p4est_find_higher_bound_overlap (&tree->quadrants, &upper_right,
                                           0.5 * tree->quadrants.elem_count);
  ub_nll = p4est_find_higher_bound_overlap (&tree->quadrants, &next_ll,
                                            0.5 * tree->quadrants.elem_count);

  if (ub_ll != -1 && ub_nll != -1) {
    P4EST_ASSERT (0 == ub_ll);
    P4EST_ASSERT (0 == ub_ib);
    P4EST_ASSERT (0 == ub_ur);
    P4EST_ASSERT (1 == ub_nll);
  }

  return 0;
}

int
main (int argc, char **argv)
{
  sc_MPI_Comm         mpicomm;
  int                 mpiret;
  int                 mpisize, mpirank;
  const int           forest_level = 4;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;

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

  /* build example forest */
#ifdef P4_TO_P8
  conn = p8est_connectivity_new_periodic ();
#else /* P4_TO_P8 */
  conn = p4est_connectivity_new_periodic ();
#endif /* P4_TO_P8 */
  p4est = p4est_new_ext (mpicomm, conn, 0, forest_level, 1, 0, 0, 0);

  /** perform tests */
  test_upper_bound (p4est);
  p4est_destroy (p4est);
  p4est = p4est_new_ext (mpicomm, conn, 0, forest_level, 0, 0, 0, 0);
  test_pos_to_quad (p4est);

  /* cleanup */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
