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

#include <p4est_virtual.h>

/* *INDENT-OFF* */
const int           p4est_face_virtual_neighbors_inside[P4EST_CHILDREN]
                                                       [P4EST_FACES] =
{{  4,  1,  6,  2 },
 {  0,  5, 10,  3 },
 {  8,  3,  0,  7 },
 {  2,  9,  1, 11 }};

const int           p4est_corner_virtual_neighbors_inside[P4EST_CHILDREN]
                                                         [P4EST_CHILDREN] =
{{ 12, 10,  8,  3 },
 {  6, 13,  2,  9 },
 {  4,  1, 14, 11 },
 {  0,  5,  7, 15 }};
/* *INDENT-ON* */
