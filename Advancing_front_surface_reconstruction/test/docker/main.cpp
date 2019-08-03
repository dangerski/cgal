#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Surface_mesh.h>
#include <boost/function_output_iterator.hpp>
#include <CGAL/squared_distance_3.h>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point3 = Kernel::Point_3;
using CGSurfaceMesh = CGAL::Surface_mesh<Kernel::Point_3>;
using Facet = CGAL::cpp11::array<std::size_t, 3>;

namespace fs = boost::filesystem;
using std::cout;
using std::endl;


struct Construct {
	CGSurfaceMesh& mesh;
	template < typename PointIterator>
	Construct(CGSurfaceMesh& mesh,PointIterator b, PointIterator e) : mesh(mesh)
	{
		for(; b!=e; ++b){
			boost::graph_traits<CGSurfaceMesh>::vertex_descriptor v;
			v = add_vertex(mesh);
			mesh.point(v) = *b;
		}
	}
	Construct& operator=(const Facet f)
	{
		typedef boost::graph_traits<CGSurfaceMesh>::vertex_descriptor vertex_descriptor;
		typedef boost::graph_traits<CGSurfaceMesh>::vertices_size_type size_type;
		mesh.add_face(vertex_descriptor(static_cast<size_type>(f[0])),
		              vertex_descriptor(static_cast<size_type>(f[1])),
		              vertex_descriptor(static_cast<size_type>(f[2])));
		return *this;
	}
	Construct& operator*() { return *this; }
	Construct& operator++() { return *this; }
	Construct operator++(int) { return *this; }
};

struct Perimeter {
	double bound;

	Perimeter(double bound) : bound(bound){}

	template <typename AdvancingFront, typename Cell_handle>
	double operator() (const AdvancingFront& adv, Cell_handle& c, const int& index) const
	{
		if (bound != 0) {
			double d = 0;
			d = sqrt(squared_distance(c->vertex((index + 1) % 4)->point(),
			                          c->vertex((index + 2) % 4)->point()));
			if (d > bound) return adv.infinity();
			d += sqrt(squared_distance(c->vertex((index + 2) % 4)->point(),
			                           c->vertex((index + 3) % 4)->point()));
			if (d > bound) return adv.infinity();
			d += sqrt(squared_distance(c->vertex((index + 1) % 4)->point(),
			                           c->vertex((index + 3) % 4)->point()));
			if (d > bound) return adv.infinity();
		}
		// Otherwise, return usual priority value: smallest radius of
		// delaunay sphere
		return adv.smallest_radius_delaunay_sphere (c, index);
	}
};

std::vector<std::string> list_files(const std::string& directory, const std::string& regex) {
	namespace ba = boost::adaptors;

	std::vector<std::string> results;

	if (fs::is_directory(directory)) {
		const boost::regex file_matcher(regex);
		boost::smatch what;
		for (auto &entry: boost::make_iterator_range(fs::directory_iterator(directory), {})
		                  | ba::filtered(static_cast<bool (*)(const fs::path &)>(&fs::is_regular_file))
		                  | ba::filtered([&](const fs::path &path){ return boost::regex_match(path.string(), what, file_matcher); })
				)
		{
			results.push_back(entry.path().string());
		}
	} else {
		cout << " This is not a directory: " << directory << endl;
	}

	std::sort(results.begin(), results.end());
	return results;
}

int main() {

	std::string base_path{"/home/afsr/"};
	std::vector<Kernel::Point_3> points_list;
	CGSurfaceMesh mesh;

	std::string regex_xyzb{".*\\.xyzb"};
	std::vector<std::string> file_names = list_files(base_path, regex_xyzb);
	for (const auto& filename : file_names) {
		points_list.clear();
		mesh.clear();

		cout << "reading file: " << filename << endl;

		std::ifstream fin;
		fin.open(filename, (std::ios::in | std::ios::binary));
		if (!fin.is_open()) {
			cout << "Error opening file " << filename << endl;
			continue;
		}

		{
			std::size_t num_points = 0;
			fin.read((char*)&num_points, sizeof(std::size_t));
			std::cout << "number of points: " << num_points << std::endl;
			points_list.reserve(num_points);
			for(std::size_t i=0; i < num_points; ++i) {
				Point3::FT x, y, z;
				fin.read((char*)&x, sizeof(Point3::FT));
				fin.read((char*)&y, sizeof(Point3::FT));
				fin.read((char*)&z, sizeof(Point3::FT));
				points_list.emplace_back(Point3(x,y,z));
			}
		}

		fin.close();
		
		Construct construct(mesh, points_list.begin(), points_list.end());
		Perimeter perimeter(3.5);
		std::cout << "Advancing front..." << std::endl;
		CGAL::advancing_front_surface_reconstruction(points_list.begin(), points_list.end(), construct, perimeter);
		std::cout << "done." << std::endl;
	}

}