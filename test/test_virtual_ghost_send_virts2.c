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
check_virtual_ghost (sc_MPI_Comm mpicomm)
{
  int8_t             *received_mirror_flags;
  p4est_connect_type_t btype = P4EST_CONNECT_FULL;
  p4est_t            *p4est;
  p4est_connectivity_t *conn;
  p4est_ghost_t      *ghost;
  p4est_mesh_t       *mesh;
  p4est_virtual_t    *virtual_quads;
  p4est_virtual_ghost_t *virtual_ghost;
  int                 min_level = 3;

  /* setup p4est */
  conn = p4est_connectivity_new_periodic ();
  p4est = p4est_new_ext (mpicomm, conn, 0, min_level, 0, 0, NULL, NULL);
  p4est_balance (p4est, btype, NULL);

  /* setup anything else */
  ghost = p4est_ghost_new (p4est, btype);
  mesh = p4est_mesh_new_ext (p4est, ghost, 1, 0, 1, btype);
  virtual_quads = p4est_virtual_new (p4est, ghost, mesh, btype);
  virtual_ghost =
    p4est_virtual_ghost_new (p4est, ghost, mesh, virtual_quads, btype);
  received_mirror_flags =
    P4EST_ALLOC (int8_t, ghost->mirror_proc_offsets[p4est->mpisize]);

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
