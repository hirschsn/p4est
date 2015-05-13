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

#ifdef P4_TO_P8
#include <p8est_vtk.h>
#include <p8est_nodes.h>
#define P4EST_VTK_CELL_TYPE     11      /* VTK_VOXEL */
#else
#include <p4est_vtk.h>
#include <p4est_nodes.h>
#define P4EST_VTK_CELL_TYPE      8      /* VTK_PIXEL */
#endif /* !P4_TO_P8 */

#include <sc_io.h>

static const double p4est_vtk_scale = 0.95;
static const int    p4est_vtk_write_tree = 1;
static const int    p4est_vtk_write_level = 1;
static const int    p4est_vtk_write_rank = 1;
static const int    p4est_vtk_wrap_rank = 0;

#ifndef P4EST_VTK_DOUBLES
#define P4EST_VTK_FLOAT_NAME "Float32"
#define P4EST_VTK_FLOAT_TYPE float
#else
#define P4EST_VTK_FLOAT_NAME "Float64"
#define P4EST_VTK_FLOAT_TYPE double
#endif

#ifndef P4EST_VTK_BINARY
#define P4EST_VTK_ASCII 1
#define P4EST_VTK_FORMAT_STRING "ascii"
#else
#define P4EST_VTK_FORMAT_STRING "binary"

static int
p4est_vtk_write_binary (FILE * vtkfile, char *numeric_data,
                        size_t byte_length)
{
#ifndef P4EST_VTK_COMPRESSION
  return sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
#else
  return sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
#endif /* P4EST_VTK_COMPRESSION */
}

#endif /* P4EST_VTK_BINARY */

/** Opaque context type for writing VTK output with multiple function calls.
 *
 * This structure holds all the information needed for the p4est vtk context.
 * It is used to relay necessary vtk information to the \b p4est_vtk_write_*
 * functions. This structure is initialized by \b p4est_vtk_write_header and
 * destroyed by \b p4est_vtk_write_footer; it can also be destroyed manually
 * using the \b p4est_vtk_context_destroy function if necessary.
 *
 * The \a p4est member is a pointer to the local p4est.
 * The \a geom member is a pointer to the geometry used to create the p4est.
 * The \a num_nodes member holds the number of nodes present in the vtk output;
 * this is determined in \b p4est_vtk_write_header using the \a scale parameter
 * and is used to assure the proper number of point variables are provided.
 * The \a filename member holds the vtk file basename: for error reporting.
 * The \a vtufilename, \a pvtufilename, and \a visitfilename members are the
 * vtk file names.
 * The \a vtufile, \a pvtufile, and \a visitfile members are the vtk file
 * pointers; opened by \b p4est_vtk_write_header and closed by \b
 * p4est_vtk_write_footer.
 *
 */
struct p4est_vtk_context
{
  p4est_t            *p4est;       /**< The p4est structure must be alive. */
  p4est_geometry_t   *geom;        /**< The geometry may be NULL. */
  p4est_locidx_t      num_nodes;   /**< Number of Q1 dofs passed. */
  const char         *filename;    /**< Original filename provided. */
  char                vtufilename[BUFSIZ];   /**< Each process writes one. */
  char                pvtufilename[BUFSIZ];  /**< Only root writes this one. */
  char                visitfilename[BUFSIZ]; /**< Only root writes this one. */
  FILE               *vtufile;     /**< File pointer for the VTU file. */
  FILE               *pvtufile;    /**< Paraview meta file. */
  FILE               *visitfile;   /**< Visit meta file. */
};

void
p4est_vtk_context_destroy (p4est_vtk_context_t * context)
{
  /* Close all file pointers. */
  P4EST_ASSERT (context->vtufile != NULL);
  if (fclose (context->vtufile))
    P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                   context->vtufilename);
  context->vtufile = NULL;

  /* Only the root process opens/closes these files. */
  if (context->p4est->mpirank == 0) {

    /* Close paraview master file */
    P4EST_ASSERT (context->pvtufile != NULL);
    if (fclose (context->pvtufile))
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->pvtufilename);
    context->pvtufile = NULL;

    /* Close visit master file */
    P4EST_ASSERT (context->visitfile != NULL);
    if (fclose (context->visitfile))
      P4EST_LERRORF (P4EST_STRING "_vtk: Error closing <%s>.\n",
                     context->visitfilename);
    context->visitfile = NULL;
  }

  /* Free context structure. */
  P4EST_FREE (context);
}

void
p4est_vtk_write_file (p4est_t * p4est, p4est_geometry_t * geom,
                      const char *filename)
{
  int                 retval;
  p4est_vtk_context_t *cont = NULL;

  cont = p4est_vtk_write_header (p4est, geom, p4est_vtk_scale, filename);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing header");

  /* Write the tree/level/rank data. */
  cont =
    p4est_vtk_write_cell_data (cont, p4est_vtk_write_tree,
                               p4est_vtk_write_level, p4est_vtk_write_rank,
                               p4est_vtk_wrap_rank, 0, 0, cont);
  SC_CHECK_ABORT (cont != NULL, P4EST_STRING "_vtk: Error writing cell data");

  retval = p4est_vtk_write_footer (cont);
  SC_CHECK_ABORT (!retval, P4EST_STRING "_vtk: Error writing footer");
}

p4est_vtk_context_t *
p4est_vtk_write_header (p4est_t * p4est, p4est_geometry_t * geom,
                        double scale, const char *filename)
{
  /* Allocate, initialize the vtk context. */
  p4est_vtk_context_t *cont = P4EST_ALLOC_ZERO (p4est_vtk_context_t, 1);

  cont->filename = filename;
  cont->p4est = p4est;
  cont->geom = geom;
  /*
   * TODO: [ehereth 2014-12-25 02:18:43]
   *       Move all context file initialization here for clarity.
   */

  p4est_connectivity_t *connectivity = p4est->connectivity;
  sc_array_t         *trees = p4est->trees;
  const int           mpirank = p4est->mpirank;
  const double        intsize = 1.0 / P4EST_ROOT_LEN;
  const double       *v = connectivity->vertices;
  const p4est_topidx_t first_local_tree = p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = p4est->last_local_tree;
  const p4est_topidx_t *tree_to_vertex = connectivity->tree_to_vertex;
  const p4est_locidx_t Ncells = p4est->local_num_quadrants;
  const p4est_locidx_t Ncorners = P4EST_CHILDREN * Ncells;
#ifdef P4EST_VTK_ASCII
  double              wx, wy, wz;
  p4est_locidx_t      sk;
#else
  int                 retval;
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  int                 xi, yi, j, k;
#ifdef P4_TO_P8
  int                 zi;
#endif
  double              h2, eta_x, eta_y, eta_z = 0.;
  double              xyz[3], XYZ[3];   /* 3 not P4EST_DIM */
  size_t              num_quads, zz;
  p4est_topidx_t      jt;
  p4est_topidx_t      vt[P4EST_CHILDREN];
  p4est_locidx_t      quad_count, Ntotal;
  p4est_locidx_t      il;
  P4EST_VTK_FLOAT_TYPE *float_data;
  sc_array_t         *quadrants, *indeps;
  p4est_tree_t       *tree;
  p4est_quadrant_t   *quad;
  p4est_nodes_t      *nodes;
  p4est_indep_t      *in;

  SC_CHECK_ABORT (p4est->connectivity->num_vertices > 0,
                  "Must provide connectivity with vertex information");

  P4EST_ASSERT (0. <= scale && scale <= 1.);
  P4EST_ASSERT (v != NULL && tree_to_vertex != NULL);

  if (scale < 1.) {
    /* when we scale the quadrants we need each corner separately */
    nodes = NULL;
    indeps = NULL;
    Ntotal = Ncorners;
  }
  else {
    /* when scale == 1. we can reuse shared quadrant corners */
    nodes = p4est_nodes_new (p4est, NULL);
    indeps = &nodes->indep_nodes;
    Ntotal = nodes->num_owned_indeps;
    P4EST_ASSERT ((size_t) Ntotal == indeps->elem_count);
  }

  /* Store the number of nodes in the vtk context. */
  cont->num_nodes = Ntotal;

  /* Have each proc write to its own file */
  snprintf (cont->vtufilename, BUFSIZ, "%s_%04d.vtu", filename, mpirank);
  /* Use "w" for writing the initial part of the file.
   * For further parts, use "r+" and fseek so write_compressed succeeds.
   */
  cont->vtufile = fopen (cont->vtufilename, "wb");
  if (cont->vtufile == NULL) {
    P4EST_LERRORF ("Could not open %s for output\n", cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    if (nodes != NULL)
      p4est_nodes_destroy (nodes);

    return NULL;
  }

  fprintf (cont->vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (cont->vtufile,
           "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_VTK_BINARY && defined P4EST_VTK_COMPRESSION
  fprintf (cont->vtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
  fprintf (cont->vtufile, " byte_order=\"BigEndian\">\n");
#else
  fprintf (cont->vtufile, " byte_order=\"LittleEndian\">\n");
#endif
  fprintf (cont->vtufile, "  <UnstructuredGrid>\n");
  fprintf (cont->vtufile,
           "    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n",
           (long long) Ntotal, (long long) Ncells);
  fprintf (cont->vtufile, "      <Points>\n");

  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, 3 * Ntotal);

  /* write point position data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"Position\""
           " NumberOfComponents=\"3\" format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);

  if (nodes == NULL) {
    /* loop over the trees */
    for (jt = first_local_tree, quad_count = 0; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;

      /* retrieve corners of the tree */
      for (k = 0; k < P4EST_CHILDREN; ++k)
        vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];

      /* loop over the elements in tree and calculate vertex coordinates */
      for (zz = 0; zz < num_quads; ++zz, ++quad_count) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        h2 = .5 * intsize * P4EST_QUADRANT_LEN (quad->level);
        k = 0;
#ifdef P4_TO_P8
        for (zi = 0; zi < 2; ++zi) {
          eta_z = intsize * quad->z + h2 * (1. + (zi * 2 - 1) * scale);
#endif
          for (yi = 0; yi < 2; ++yi) {
            eta_y = intsize * quad->y + h2 * (1. + (yi * 2 - 1) * scale);
            for (xi = 0; xi < 2; ++xi) {
              P4EST_ASSERT (0 <= k && k < P4EST_CHILDREN);
              eta_x = intsize * quad->x + h2 * (1. + (xi * 2 - 1) * scale);
              if (geom != NULL) {
                xyz[0] = eta_x;
                xyz[1] = eta_y;
                xyz[2] = eta_z;
                geom->X (geom, jt, xyz, XYZ);
                for (j = 0; j < 3; ++j) {
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
                    (P4EST_VTK_FLOAT_TYPE) XYZ[j];
                }
              }
              else {
                for (j = 0; j < 3; ++j) {
                  /* *INDENT-OFF* */
                xyz[j] =
            ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                   eta_x  * v[3 * vt[1] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                   eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
             +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                   eta_x  * v[3 * vt[5] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                   eta_x  * v[3 * vt[7] + j]))
#endif
            );
                  /* *INDENT-ON* */
                  float_data[3 * (P4EST_CHILDREN * quad_count + k) + j] =
                    (P4EST_VTK_FLOAT_TYPE) xyz[j];
                }
              }
              ++k;
            }
          }
#ifdef P4_TO_P8
        }
#endif
        P4EST_ASSERT (k == P4EST_CHILDREN);
      }
    }
    P4EST_ASSERT (P4EST_CHILDREN * quad_count == Ntotal);
  }
  else {
    for (zz = 0; zz < indeps->elem_count; ++zz) {
      in = (p4est_indep_t *) sc_array_index (indeps, zz);

      /* retrieve corners of the tree */
      jt = in->p.which_tree;
      for (k = 0; k < P4EST_CHILDREN; ++k)
        vt[k] = tree_to_vertex[jt * P4EST_CHILDREN + k];

      /* calculate vertex coordinates */
      eta_x = intsize * in->x;
      eta_y = intsize * in->y;
#ifdef P4_TO_P8
      eta_z = intsize * in->z;
#endif
      if (geom != NULL) {
        xyz[0] = eta_x;
        xyz[1] = eta_y;
        xyz[2] = eta_z;
        geom->X (geom, jt, xyz, XYZ);
        for (j = 0; j < 3; ++j) {
          float_data[3 * zz + j] = (P4EST_VTK_FLOAT_TYPE) XYZ[j];
        }
      }
      else {
        for (j = 0; j < 3; ++j) {
          /* *INDENT-OFF* */
          xyz[j] =
            ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                   eta_x  * v[3 * vt[1] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                   eta_x  * v[3 * vt[3] + j]))
  #ifdef P4_TO_P8
             +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                   eta_x  * v[3 * vt[5] + j]) +
                                   eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                   eta_x  * v[3 * vt[7] + j]))
  #endif
            );
          /* *INDENT-ON* */

          float_data[3 * zz + j] = (P4EST_VTK_FLOAT_TYPE) xyz[j];
        }
      }
    }
  }

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ntotal; ++il) {
    wx = float_data[3 * il + 0];
    wy = float_data[3 * il + 1];
    wz = float_data[3 * il + 2];

#ifdef P4EST_VTK_DOUBLES
    fprintf (cont->vtufile, "     %24.16e %24.16e %24.16e\n", wx, wy, wz);
#else
    fprintf (cont->vtufile, "          %16.8e %16.8e %16.8e\n", wx, wy, wz);
#endif
  }
#else
  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * 3 * Ntotal);
  fprintf (cont->vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (float_data);
    if (nodes != NULL)
      p4est_nodes_destroy (nodes);

    return NULL;
  }
#endif
  P4EST_FREE (float_data);

  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Points>\n");
  fprintf (cont->vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (cont->vtufile,
           "        <DataArray type=\"%s\" Name=\"connectivity\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  for (sk = 0, il = 0; il < Ncells; ++il) {
    fprintf (cont->vtufile, "         ");
    for (k = 0; k < P4EST_CHILDREN; ++sk, ++k) {
      fprintf (cont->vtufile, " %lld", nodes == NULL ?
               (long long) sk : (long long) nodes->local_nodes[sk]);
    }
    fprintf (cont->vtufile, "\n");
  }
#else
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncorners);
  fprintf (cont->vtufile, "          ");
  if (nodes == NULL) {
    for (il = 0; il < Ncorners; ++il) {
      locidx_data[il] = il;
    }
    retval =
      p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                              sizeof (*locidx_data) * Ncorners);
  }
  else {
    retval =
      p4est_vtk_write_binary (cont->vtufile, (char *) nodes->local_nodes,
                              sizeof (*locidx_data) * Ncorners);
  }
  fprintf (cont->vtufile, "\n");
  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding connectivity\n");
    p4est_vtk_context_destroy (cont);

    if (nodes != NULL)
      p4est_nodes_destroy (nodes);
    P4EST_FREE (locidx_data);

    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"offsets\""
           " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %lld", (long long) (P4EST_CHILDREN * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  for (il = 1; il <= Ncells; ++il)
    locidx_data[il - 1] = P4EST_CHILDREN * il;  /* same type */

  fprintf (cont->vtufile, "          ");
  retval = p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                   sizeof (*locidx_data) * Ncells);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (locidx_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding offsets\n");
    p4est_vtk_context_destroy (cont);

    if (nodes != NULL)
      p4est_nodes_destroy (nodes);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (cont->vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
  fprintf (cont->vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (cont->vtufile, " %d", P4EST_VTK_CELL_TYPE);
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (cont->vtufile, "\n         ");
  }
  fprintf (cont->vtufile, "\n");
#else
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
  for (il = 0; il < Ncells; ++il)
    uint8_data[il] = P4EST_VTK_CELL_TYPE;

  fprintf (cont->vtufile, "          ");
  retval = p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                                   sizeof (*uint8_data) * Ncells);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (uint8_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
    p4est_vtk_context_destroy (cont);

    if (nodes != NULL)
      p4est_nodes_destroy (nodes);

    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");
  fprintf (cont->vtufile, "      </Cells>\n");

  if (nodes != NULL) {
    p4est_nodes_destroy (nodes);
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing header\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    snprintf (cont->pvtufilename, BUFSIZ, "%s.pvtu", filename);

    cont->pvtufile = fopen (cont->pvtufilename, "wb");
    if (!cont->pvtufile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->pvtufilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    fprintf (cont->pvtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (cont->pvtufile,
             "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined P4EST_VTK_BINARY && defined P4EST_VTK_COMPRESSION
    fprintf (cont->pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_IS_BIGENDIAN
    fprintf (cont->pvtufile, " byte_order=\"BigEndian\">\n");
#else
    fprintf (cont->pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

    fprintf (cont->pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf (cont->pvtufile, "    <PPoints>\n");
    fprintf (cont->pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
             " NumberOfComponents=\"3\" format=\"%s\"/>\n",
             P4EST_VTK_FLOAT_NAME, P4EST_VTK_FORMAT_STRING);
    fprintf (cont->pvtufile, "    </PPoints>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);
      return NULL;
    }

    /* Create a master file for visualization in Visit; this will be used
     * only in p4est_vtk_write_footer().
     */
    snprintf (cont->visitfilename, BUFSIZ, "%s.visit", filename);
    cont->visitfile = fopen (cont->visitfilename, "wb");
    if (!cont->visitfile) {
      P4EST_LERRORF ("Could not open %s for output\n", cont->visitfilename);
      p4est_vtk_context_destroy (cont);
      return NULL;
    }
  }

  return cont;
}

/** Write VTK point data.
 *
 * This function exports custom point data to the vtk file; it is functionally
 * the same as \b p4est_vtk_write_point_data with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b p4est_vtk_write_point_data
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static p4est_vtk_context_t *
p4est_vtk_write_point_datav (p4est_vtk_context_t * cont,
                             int num_point_scalars,
                             int num_point_vectors, va_list ap)
{
  /* This function needs to do nothing if there is no data. */
  if (!(num_point_scalars || num_point_vectors))
    return cont;

  int                 retval;
  int                 i, all = 0;
  const int           mpirank = cont->p4est->mpirank;
  int                 scalar_strlen, vector_strlen;
  char                point_scalars[BUFSIZ], point_vectors[BUFSIZ];
  const char         *name, **names;
  sc_array_t        **values;

  values = P4EST_ALLOC (sc_array_t *, num_point_scalars + num_point_vectors);
  names = P4EST_ALLOC (const char *, num_point_scalars + num_point_vectors);

  /* Gather point data. */
  scalar_strlen = 0;
  point_scalars[0] = '\0';
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data type; scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == cont->num_nodes,
                    P4EST_STRING
                    "_vtk: Error: incorrect point scalar data count; see "
                    P4EST_STRING ".h for more details.");
  }

  vector_strlen = 0;
  point_vectors[0] = '\0';
  for (i = 0; i < num_point_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (point_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting point vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data type; vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count == (3 * cont->num_nodes),
                    P4EST_STRING
                    "_vtk: Error: incorrect point vector data count; see "
                    P4EST_STRING ".h for more details.");
  }

  /* Check for pointer variable marking the end of variable data input. */
  p4est_vtk_context_t *end = va_arg (ap, p4est_vtk_context_t *);
  SC_CHECK_ABORT (end == cont, P4EST_STRING "_vtk Error: the end of variable "
                  "data must be specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t struct. See " P4EST_STRING "_vtk.h for more information.");

  fprintf (cont->vtufile, "      <PointData");
  if (point_scalars != NULL)
    fprintf (cont->vtufile, " Scalars=\"%s\"", point_scalars);
  if (point_vectors != NULL)
    fprintf (cont->vtufile, " Vectors=\"%s\"", point_vectors);
  fprintf (cont->vtufile, ">\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  all = 0;
  for (i = 0; i < num_point_scalars; ++all, ++i) {
    cont = p4est_vtk_write_point_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point scalars");
  }

  for (i = 0; i < num_point_vectors; ++all, ++i) {
    cont = p4est_vtk_write_point_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing point vectors");
  }

  fprintf (cont->vtufile, "      </PointData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PPointData>\n");

    all = 0;
    for (i = 0; i < num_point_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_point_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PPointData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_point_data (p4est_vtk_context_t * cont,
                            int num_point_scalars, int num_point_vectors, ...)
{
  va_list             ap;
  va_start (ap, num_point_vectors);

  cont = p4est_vtk_write_point_datav (cont, num_point_scalars,
                                      num_point_vectors, ap);
  va_end (ap);

  return cont;
}

/** Write VTK cell data.
 *
 * This function exports custom cell data to the vtk file; it is functionally
 * the same as \b p4est_vtk_write_cell_data with the only difference being
 * that instead of a variable argument list, an initialized \a va_list is
 * passed as the last argument. The \a va_list is initialized from the variable
 * argument list of the calling function.
 *
 * \note This function is actually called from \b p4est_vtk_write_cell_data
 * and does all of the work.
 *
 * \param [in,out] cont    A vtk context created by p4est_vtk_write_header.
 * \param [in] num_point_scalars Number of point scalar datasets to output.
 * \param [in] num_point_vectors Number of point vector datasets to output.
 * \param [in,out] ap      An initialized va_list used to access the
 *                         scalar/vector data.
 *
 * \return          On success, the context that has been passed in.
 *                  On failure, returns NULL and deallocates the context.
 */
static p4est_vtk_context_t *
p4est_vtk_write_cell_datav (p4est_vtk_context_t * cont,
                            int write_tree, int write_level,
                            int write_rank, int wrap_rank,
                            int num_cell_scalars,
                            int num_cell_vectors, va_list ap)
{
  /* This function needs to do nothing if there is no data. */
  if (!
      (write_tree || write_level || write_rank || wrap_rank
       || num_cell_vectors || num_cell_vectors))
    return cont;

  const int           mpirank = cont->p4est->mpirank;
  int                 retval;
  int                 i, all = 0;
  int                 scalar_strlen, vector_strlen;
  sc_array_t         *trees = cont->p4est->trees;
  p4est_tree_t       *tree;
  const p4est_topidx_t first_local_tree = cont->p4est->first_local_tree;
  const p4est_topidx_t last_local_tree = cont->p4est->last_local_tree;
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  char                cell_scalars[BUFSIZ], cell_vectors[BUFSIZ];
  const char         *name, **names;
  sc_array_t        **values;
  size_t              num_quads, zz;
  sc_array_t         *quadrants;
  p4est_quadrant_t   *quad;
#ifdef P4EST_VTK_ASCII
  p4est_locidx_t      sk;
#else
  uint8_t            *uint8_data;
  p4est_locidx_t     *locidx_data;
#endif
  p4est_topidx_t      jt;
  p4est_locidx_t      il;

  P4EST_ASSERT (wrap_rank >= 0);

  values = P4EST_ALLOC (sc_array_t *, num_cell_scalars + num_cell_vectors);
  names = P4EST_ALLOC (const char *, num_cell_scalars + num_cell_vectors);

  /* Gather cell data. */
  scalar_strlen = 0;
  cell_scalars[0] = '\0';
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_scalars + scalar_strlen, BUFSIZ - scalar_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell scalars");
    scalar_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data type; scalar data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell scalar data count; scalar data must contain exactly p4est->local_num_quadrants doubles.");
  }

  vector_strlen = 0;
  cell_vectors[0] = '\0';
  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    name = names[all] = va_arg (ap, const char *);
    retval = snprintf (cell_vectors + vector_strlen, BUFSIZ - vector_strlen,
                       "%s%s", i == 0 ? "" : ",", name);
    SC_CHECK_ABORT (retval > 0,
                    P4EST_STRING "_vtk: Error collecting cell vectors");
    vector_strlen += retval;
    values[all] = va_arg (ap, sc_array_t *);

    /* Validate input. */
    SC_CHECK_ABORT (values[all]->elem_size == sizeof (double),
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data type; vector data must contain doubles.");
    SC_CHECK_ABORT (values[all]->elem_count ==
                    3 * cont->p4est->local_num_quadrants,
                    P4EST_STRING
                    "_vtk: Error: incorrect cell vector data count; vector data must contain exactly 3*p4est->local_num_quadrants doubles.");
  }

  /* Check for pointer variable marking the end of variable data input. */
  p4est_vtk_context_t *end = va_arg (ap, p4est_vtk_context_t *);
  SC_CHECK_ABORT (end == cont, P4EST_STRING "_vtk Error: the end of variable "
                  "data must be specified by passing, as the last argument, the current "
                  P4EST_STRING "_vtk_context_t struct. See " P4EST_STRING "_vtk.h for more information.");

  char                vtkCellDataString[BUFSIZ] = "";
  int                 printed = 0;

  if (write_tree)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed, "treeid");

  if (write_level)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",level" : "level");

  if (write_rank)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",mpirank" : "mpirank");

  if (num_cell_scalars)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_scalars);

  if (num_cell_vectors)
    printed +=
      snprintf (vtkCellDataString + printed, BUFSIZ - printed,
                printed > 0 ? ",%s" : "%s", cell_vectors);

  fprintf (cont->vtufile, "      <CellData Scalars=\"%s\">\n",
           vtkCellDataString);

#ifndef P4EST_VTK_ASCII
  locidx_data = P4EST_ALLOC (p4est_locidx_t, Ncells);
  uint8_data = P4EST_ALLOC (uint8_t, Ncells);
#endif

  if (write_tree) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"treeid\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        fprintf (cont->vtufile, " %lld", (long long) jt);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      num_quads = tree->quadrants.elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        locidx_data[il] = (p4est_locidx_t) jt;
      }
    }
    fprintf (cont->vtufile, "          ");
    retval = p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (cont->vtufile, "\n");
    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);
      P4EST_FREE (uint8_data);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
    P4EST_ASSERT (il == Ncells);
  }

  if (write_level) {
    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"level\""
             " format=\"%s\">\n", "UInt8", P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++sk, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        fprintf (cont->vtufile, " %d", (int) quad->level);
        if (!(sk % 20) && il != (Ncells - 1))
          fprintf (cont->vtufile, "\n         ");
      }
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0, jt = first_local_tree; jt <= last_local_tree; ++jt) {
      tree = p4est_tree_array_index (trees, jt);
      quadrants = &tree->quadrants;
      num_quads = quadrants->elem_count;
      for (zz = 0; zz < num_quads; ++zz, ++il) {
        quad = p4est_quadrant_array_index (quadrants, zz);
        uint8_data[il] = (uint8_t) quad->level;
      }
    }

    fprintf (cont->vtufile, "          ");
    retval = p4est_vtk_write_binary (cont->vtufile, (char *) uint8_data,
                                     sizeof (*uint8_data) * Ncells);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (uint8_data);

    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);
      P4EST_FREE (locidx_data);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (write_rank) {
    const int           wrapped_rank =
      wrap_rank > 0 ? mpirank % wrap_rank : mpirank;

    fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"mpirank\""
             " format=\"%s\">\n", P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);
#ifdef P4EST_VTK_ASCII
    fprintf (cont->vtufile, "         ");
    for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
      fprintf (cont->vtufile, " %d", wrapped_rank);
      if (!(sk % 20) && il != (Ncells - 1))
        fprintf (cont->vtufile, "\n         ");
    }
    fprintf (cont->vtufile, "\n");
#else
    for (il = 0; il < Ncells; ++il)
      locidx_data[il] = (p4est_locidx_t) wrapped_rank;

    fprintf (cont->vtufile, "          ");
    retval = p4est_vtk_write_binary (cont->vtufile, (char *) locidx_data,
                                     sizeof (*locidx_data) * Ncells);
    fprintf (cont->vtufile, "\n");

    P4EST_FREE (locidx_data);

    if (retval) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error encoding types\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (values);
      P4EST_FREE (names);

      return NULL;
    }
#endif
    fprintf (cont->vtufile, "        </DataArray>\n");
  }

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (values);
    P4EST_FREE (names);

    return NULL;
  }

  all = 0;
  for (i = 0; i < num_cell_scalars; ++all, ++i) {
    cont = p4est_vtk_write_cell_scalar (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell scalars");
  }

  for (i = 0; i < num_cell_vectors; ++all, ++i) {
    cont = p4est_vtk_write_cell_vector (cont, names[all], values[all]);
    SC_CHECK_ABORT (cont != NULL,
                    P4EST_STRING "_vtk: Error writing cell vectors");
  }

  fprintf (cont->vtufile, "      </CellData>\n");

  P4EST_FREE (values);

  if (ferror (cont->vtufile)) {
    P4EST_LERRORF (P4EST_STRING "_vtk: Error writing %s\n",
                   cont->vtufilename);
    p4est_vtk_context_destroy (cont);

    P4EST_FREE (names);

    return NULL;
  }

  /* Only have the root write to the parallel vtk file */
  if (mpirank == 0) {
    fprintf (cont->pvtufile, "    <PCellData Scalars=\"%s\">\n",
             vtkCellDataString);

    if (write_tree)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"treeid\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);

    if (write_level)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"level\" format=\"%s\"/>\n",
               "UInt8", P4EST_VTK_FORMAT_STRING);

    if (write_rank)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"mpirank\" format=\"%s\"/>\n",
               P4EST_VTK_LOCIDX, P4EST_VTK_FORMAT_STRING);

    all = 0;
    for (i = 0; i < num_cell_scalars; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    for (i = 0; i < num_cell_vectors; ++all, i++)
      fprintf (cont->pvtufile, "      "
               "<PDataArray type=\"%s\" Name=\"%s\" format=\"%s\"/>\n",
               P4EST_VTK_FLOAT_NAME, names[all], P4EST_VTK_FORMAT_STRING);

    fprintf (cont->pvtufile, "    </PCellData>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel header\n");
      p4est_vtk_context_destroy (cont);

      P4EST_FREE (names);

      return NULL;
    }
  }

  P4EST_FREE (names);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_data (p4est_vtk_context_t * cont,
                           int write_tree, int write_level,
                           int write_rank, int wrap_rank,
                           int num_cell_scalars, int num_cell_vectors, ...)
{
  va_list             ap;

  va_start (ap, num_cell_vectors);
  cont = p4est_vtk_write_cell_datav (cont,
                                     write_tree, write_level,
                                     write_rank, wrap_rank,
                                     num_cell_scalars, num_cell_vectors, ap);
  va_end (ap);

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_point_scalar (p4est_vtk_context_t * cont,
                              const char *scalar_name,
                              sc_array_t * values)
{
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  const p4est_locidx_t Ncorners = P4EST_CHILDREN * Ncells;      /* type ok */
  p4est_locidx_t      il;
#ifndef P4EST_VTK_ASCII
  int                 retval;
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif

  /* write point data */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, scalar_name, P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ncorners; ++il) {
#ifdef P4EST_VTK_DOUBLES
    fprintf (cont->vtufile, "     %24.16e\n",
             *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il)));
#else
    fprintf (cont->vtufile, "          %16.8e\n",
             *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il)));
#endif
  }
#else
  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Ncorners);
  for (il = 0; il < Ncorners; ++il) {
    float_data[il] = *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il));
  }

  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * Ncorners);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (float_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding points\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing point scalar\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_point_vector (p4est_vtk_context_t * cont,
                              const char *vector_name,
                              sc_array_t * values)
{
  SC_ABORT (P4EST_STRING "_vtk_write_point_vector not implemented");
}

p4est_vtk_context_t *
p4est_vtk_write_cell_scalar (p4est_vtk_context_t * cont,
                             const char *scalar_name,
                             sc_array_t * values)
{
  const p4est_locidx_t Ncells = cont->p4est->local_num_quadrants;
  p4est_locidx_t      il;
#ifndef P4EST_VTK_ASCII
  int                 retval;
  P4EST_VTK_FLOAT_TYPE *float_data;
#endif

  /* Write cell data. */
  fprintf (cont->vtufile, "        <DataArray type=\"%s\" Name=\"%s\""
           " format=\"%s\">\n",
           P4EST_VTK_FLOAT_NAME, scalar_name, P4EST_VTK_FORMAT_STRING);

#ifdef P4EST_VTK_ASCII
  for (il = 0; il < Ncells; ++il) {
#ifdef P4EST_VTK_DOUBLES
    fprintf (cont->vtufile, "     %24.16e\n",
             *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il)));
#else
    fprintf (cont->vtufile, "          %16.8e\n",
             *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il)));
#endif
  }
#else
  float_data = P4EST_ALLOC (P4EST_VTK_FLOAT_TYPE, Ncells);
  for (il = 0; il < Ncells; ++il) {
    float_data[il] = *((P4EST_VTK_FLOAT_TYPE *) sc_array_index (values, il));
  }

  fprintf (cont->vtufile, "          ");
  /* TODO: Don't allocate the full size of the array, only allocate
   * the chunk that will be passed to zlib and do this a chunk
   * at a time.
   */
  retval = p4est_vtk_write_binary (cont->vtufile, (char *) float_data,
                                   sizeof (*float_data) * Ncells);
  fprintf (cont->vtufile, "\n");

  P4EST_FREE (float_data);

  if (retval) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error encoding scalar cell data\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }
#endif
  fprintf (cont->vtufile, "        </DataArray>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing cell scalar file\n");
    p4est_vtk_context_destroy (cont);
    return NULL;
  }

  return cont;
}

p4est_vtk_context_t *
p4est_vtk_write_cell_vector (p4est_vtk_context_t * cont,
                             const char *vector_name,
                             sc_array_t * values)
{
  SC_ABORT (P4EST_STRING "_vtk_write_cell_vector not implemented");
}

int
p4est_vtk_write_footer (p4est_vtk_context_t * cont)
{
  int                 p;
  int                 procRank = cont->p4est->mpirank;
  int                 numProcs = cont->p4est->mpisize;

  fprintf (cont->vtufile, "    </Piece>\n");
  fprintf (cont->vtufile, "  </UnstructuredGrid>\n");
  fprintf (cont->vtufile, "</VTKFile>\n");

  if (ferror (cont->vtufile)) {
    P4EST_LERROR (P4EST_STRING "_vtk: Error writing footer\n");
    p4est_vtk_context_destroy (cont);
    return -1;
  }

  /* Only have the root write to the parallel vtk file */
  if (procRank == 0) {
    fprintf (cont->visitfile, "!NBLOCKS %d\n", numProcs);

    /* Write data about the parallel pieces into both files */
    for (p = 0; p < numProcs; ++p) {
      fprintf (cont->pvtufile,
               "    <Piece Source=\"%s_%04d.vtu\"/>\n", cont->filename, p);
      fprintf (cont->visitfile, "%s_%04d.vtu\n", cont->filename, p);
    }
    fprintf (cont->pvtufile, "  </PUnstructuredGrid>\n");
    fprintf (cont->pvtufile, "</VTKFile>\n");

    if (ferror (cont->pvtufile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      p4est_vtk_context_destroy (cont);
      return -1;
    }

    if (ferror (cont->visitfile)) {
      P4EST_LERROR (P4EST_STRING "_vtk: Error writing parallel footer\n");
      p4est_vtk_context_destroy (cont);
      return -1;
    }
  }

  /* Destroy context structure. */
  p4est_vtk_context_destroy (cont);

  return 0;
}
