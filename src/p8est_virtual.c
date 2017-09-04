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

#include <p4est_to_p8est.h>
#include <p8est_virtual.h>
#include <p8est_connectivity.h>

/* *INDENT-OFF* */
const int           p8est_face_virtual_neighbors_inside[P8EST_CHILDREN]
                                                       [P8EST_FACES] =
{{  8,  1, 10,  2, 12,  4 },
 {  0,  9, 16,  3, 18,  5 },
 { 14,  3,  0, 11, 24,  6 },
 {  2, 15,  1, 17, 30,  7 },
 { 20,  5, 22,  6,  0, 13 },
 {  4, 21, 28,  7,  1, 19 },
 { 26,  7,  4, 23,  2, 25 },
 {  6, 27,  5, 29,  3, 31 }};

const int           p8est_edge_virtual_neighbors_inside[P8EST_CHILDREN]
                                                       [P8EST_EDGES] =
{{ 32, 24, 22,  6, 36, 18, 20,  5, 40, 16, 14,  3 },
 { 44, 30, 28,  7, 12, 37,  4, 21, 10, 41,  2, 15 },
 { 12, 33,  4, 23, 48, 30, 26,  7,  8,  1, 42, 17 },
 { 18, 45,  5, 29, 24, 49,  6, 27,  0,  9, 11, 43 },
 { 10,  2, 34, 25,  8,  1, 38, 19, 52, 28, 26,  7 },
 { 16,  3, 46, 31,  0,  9, 13, 39, 22, 53,  6, 27 },
 {  0, 11, 13, 35, 14,  3, 50, 31, 20,  5, 54, 29 },
 {  1, 17, 19, 47,  2, 15, 25, 51,  4, 21, 23, 55 }};

const int           p8est_corner_virtual_neighbors_inside[P8EST_CHILDREN]
                                                         [P8EST_CHILDREN] =
{{ 56, 44, 48, 30, 52, 28, 26,  7 },
 { 32, 57, 24, 49, 22, 53,  6, 27 },
 { 36, 18, 58, 45, 20,  5, 54, 29 },
 { 12, 37, 33, 59,  4, 21, 23, 55 },
 { 40, 16, 14,  3, 60, 46, 50, 31 },
 { 10, 41,  2, 15, 34, 61, 25, 51 },
 {  8,  1, 42, 17, 38, 19, 62, 47 },
 {  0,  9, 11, 43, 13, 39, 35, 63 }};
/* *INDENT-ON* */

#include "p4est_virtual.c"
