/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */


#include "photon.h"
#include "colour.h"
#include "ray.h"
#include "object.h"
#include "light.h"
#include "hit.h"
#include "vector.h"
#include "vertex.h"
#include "nanoflann-master/include/nanoflann.hpp"
//#include "nanoflann-master\\utils.h"

using namespace nanoflann;

// Scene is a class that is used to build a scene database of objects
// and lights and then trace a ray through it.


struct PointCloud
{
	struct Point
	{            
		float  x,y,z;
		Vector incident_direction;
		Colour power; 
	};

	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline float kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

};


class Scene {
public:
	
	typedef KDTreeSingleIndexAdaptor<
	L2_Simple_Adaptor<float, PointCloud > ,
	PointCloud,
	3 /* dim */
	> my_kd_tree_t;
	
	my_kd_tree_t *tree;
	my_kd_tree_t *caustic_tree;
	Object *object_list;
	Light *light_list;
	PointCloud map;
	PointCloud caustic_map;
	Scene();

	// Trace a ray through the scene and find the closest if any object
	// intersection in front of the ray.
	void trace(Ray ray, Object *objects, Hit &hit);
	
	// Trace a ray through the scene and return its colour. This function
	// is the one that should recurse down the reflection/refraction tree.
	void raytrace(Ray ray, int level, Object *objects, Light *lights, Colour &colour);
	
	// photon
    void photon_trace(Photon photon, int level, int diffuse_level, int specular_level,bool caustic);
    void emit_photon_diffuse_point(Vertex light_location,int num_photons);
    
    void hemisphere_diffuse(Vector normal,Vector &new_direction);
	void store_photon(Photon photon,bool caustic);
	void photon_raytrace(Ray ray, Vertex eye,Colour &out_colour,Vertex light_pos, int level,bool caustic,Colour &caustic_output, Colour &direct_output, Colour &global_map_output);
	void emit_photon_caustic(Vertex light_location,int num_photons,float power_ratio);
};
