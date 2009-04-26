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
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

template <class Tree, class Polyhedron, class K>
void test_speed(Tree& tree, 
                Polyhedron& polyhedron)
{
        typedef K::FT FT;
        typedef K::Ray_3 Ray;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;

        CGAL::Timer timer;
        unsigned int nb = 0;
        timer.start();
        Point source((FT)0.0, (FT)0.0, (FT)0.0);
        Vector vec((FT)0.1, (FT)0.2, (FT)0.3);
        Ray ray(source, vec);
        while(timer.time() < 1.0)
        {
                tree.do_intersect(ray); 
                nb++;
        }
        double speed = (double)nb / timer.time();
        std::cout << speed << " intersections/s" << std::endl;
        timer.stop();
}

template <class K>
void test(const char *filename)
{
        typedef K::FT FT;
        typedef K::Ray_3 Ray;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;
        typedef CGAL::Polyhedron_3<K> Polyhedron;
        typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
        typedef CGAL::AABB_traits<K, Primitive> Traits;
        typedef CGAL::AABB_tree<Traits> Polyhedron_tree;

        Polyhedron polyhedron;
        std::ifstream ifs(filename);
        ifs >> polyhedron;

        // construct tree
        Polyhedron_tree tree(polyhedron.facets_begin(),polyhedron.facets_end());

        // TODO
        // - compare tree tests with exhaustive ones
        // - query with ray/line/segment

        Point source((FT)0.5, (FT)0.5, (FT)0.5);
        Ray ray(source, Vector((FT)0.1, (FT)0.2, (FT)0.3));
        std::cout << tree.number_of_intersections(ray) 
                << " intersections(s) with ray" << std::endl;

        test_speed<Polyhedron_tree,Polyhedron,K>(tree,polyhedron);
}

void test_kernels(const char *filename)
{
        std::cout << std::endl;
        std::cout << "Polyhedron " << filename << std::endl;

        std::cout << "Simple cartesian float kernel" << std::endl;
        test<CGAL::Simple_cartesian<float> >(filename);

        std::cout << "Simple cartesian double kernel" << std::endl;
        test<CGAL::Simple_cartesian<double> >(filename);

        std::cout << "Epic kernel" << std::endl;
        test<CGAL::Exact_predicates_inexact_constructions_kernel>(filename);
}

int main(void)
{
        std::cout << "AABB intersection tests" << std::endl;
        test_kernels("../data/cube.off");
        test_kernels("../data/coverrear.off");
        test_kernels("../data/nested_spheres.off");
        return 0;
}