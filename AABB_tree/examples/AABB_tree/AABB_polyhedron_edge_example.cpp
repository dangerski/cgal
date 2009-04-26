// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_polyhedron_edge_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_edge_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Polyhedron_tree;

int main(void)
{
        Point p(1.0, 0.0, 0.0);
        Point q(0.0, 1.0, 0.0);
        Point r(0.0, 0.0, 1.0);
        Point s(0.0, 0.0, 0.0);
        Polyhedron polyhedron;
        polyhedron.make_tetrahedron( p, q, r, s);

        Polyhedron_tree tree(polyhedron.edges_begin(),polyhedron.edges_end());
        Point base(0.2, 0.2, 0.2);
        Vector normal(0.1, 0.2, 0.3);
        Plane plane(base,normal);
        std::cout << tree.number_of_intersections(plane) 
                << " intersections(s) with plane" << std::endl;
        return 0;
}