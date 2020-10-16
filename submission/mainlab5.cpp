#include "nanoflann-master/include/nanoflann.hpp"
//#include "nanoflann-master\\utils.h"
#include <ctime>
#include <cstdlib>
#include <iostream>

#include "math.h"
#include "vector.h"
#include "vertex.h"
#include "photon.h"
#include "scene.h"
#include "hit.h"
#include "material.h"
#include "polymesh.h"
#include "phong.h"
#include "sphere.h"
#include "framebuffer.h"
#include <algorithm>

using namespace std;
using namespace nanoflann;


// photon mapping lab

// this creates structure for storing photon to use with nanoflann.


int main()
{

    

    int width = 300;
    int height = 300;
    
    FrameBuffer *fb = new FrameBuffer(width,height);
    FrameBuffer *fb1 = new FrameBuffer(width,height);
    FrameBuffer *fb2 = new FrameBuffer(width,height);
    FrameBuffer *fb3 = new FrameBuffer(width,height);

    Sphere *sphere  = new Sphere(Vertex(2.5,1.2,8.0), 1.0f);
    Sphere *sphere2 = new Sphere(Vertex(0.0,1.2,8.0),1.0f);
    Sphere *sphere3 = new Sphere(Vertex(-2.5,1.2,8.0),1.0f); 
    Sphere *sphere4 = new Sphere(Vertex(0.0,-1.5,4),0.6f);// this sphere is inside teapot
   

    
    
    // TEAPOT STUFF
    
    // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
    Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
                                        0.0f, 0.0f, 1.0f, -2.7f,
                                        0.0f, 1.0f, 0.0f, 5.0f,
                                        0.0f, 0.0f, 0.0f, 1.0f); 

    //  Read in the teapot model.
    PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);
    PolyMesh *plane = new PolyMesh((char *)"plane.ply");
    PolyMesh *plane2 = new PolyMesh((char *)"plane2.ply");
    
    sphere->next = sphere2;
    sphere2->next = sphere3;
    sphere3->next = sphere4;
    sphere4->next = pm;
    
    pm->next = plane;
    //plane->next = plane2;
    //splane2->next = plane;
    //pm->next = plane;
    //sphere4->next = plane;
    

    Phong red;
    // Don't need ambient for photon mapping
 
    red.diffuse_photon = Colour(1,0,0,0);
    red.diffuse = red.diffuse_photon;
    red.specular_photon = Colour(1,0,0,0);
    red.specular = red.specular_photon;
    red.power = 20;
    red.index_refrac = 1.01;
    red.reflect = 1.0f;
    red.refract = 0.0f;


    Phong blue;
    // Don't need ambient for photon mapping

    blue.diffuse_photon = Colour(0,0,1,0);
    blue.diffuse = blue.diffuse_photon;
    blue.specular_photon = Colour(0,0,1,0);
    blue.specular = blue.specular_photon;
    blue.power =20;
    blue.index_refrac = 1.01f;
    blue.reflect = 1.0f;
    blue.refract = 0.0f;
    Phong green;
    
    green.diffuse_photon = Colour(0,1,0,0);
    green.diffuse = green.diffuse_photon;
    green.specular_photon = Colour(0,1,0,0);
    green.specular = green.specular_photon;
    green.power = 20;
    green.index_refrac = 1.01f;
    green.reflect = 1.0f;
    green.refract = 0.0f;
    Phong white;

    white.diffuse_photon = Colour(0.0f,0.0f,0.0f,0.0f); // should be 0 if totally transparent.
    white.diffuse = white.diffuse_photon;
    white.specular_photon = Colour(1.0f,1.0f,1.0f,0.0f);
    white.specular = white.specular_photon;
    white.power = 40;
    white.index_refrac = 1.1f;
    white.reflect = 0.5f; // change back to 1.0
    white.refract = 0.5f; // should both be 1 and work out actual kt and kr values from index of refraction.

    Phong white_reflect;
    white_reflect.diffuse_photon = Colour(0.6f,0.6f,0.6f,0.0f);
    white_reflect.diffuse = white_reflect.diffuse_photon;
    white_reflect.specular_photon = Colour(1.0f,1.0f,1.0f,1.0f);
    white_reflect.specular = white_reflect.specular_photon;
    white_reflect.power = 40;
    white_reflect.index_refrac = 1.1f;
    white_reflect.reflect = 1.0f;
    white_reflect.refract = 0.0;

    Phong no_trans_white;
    no_trans_white.diffuse_photon = Colour(1.0f,1.0f,1.0f,0.0f);
    no_trans_white.diffuse = no_trans_white.diffuse_photon;
    no_trans_white.specular_photon = Colour(0.0f,0.0f,0.0f,0.0f);
    no_trans_white.specular = no_trans_white.specular_photon;
    no_trans_white.power = 40;
    no_trans_white.index_refrac = 1.1f;
    no_trans_white.reflect = 0.0f;
    no_trans_white.refract = 0.0f;



    /*
    left->material = &green;
    ceiling->material = &red;
    floor->material = &blue;
    back->material = &red;
    right->material = &blue;
    */

    sphere->material=&blue;
    sphere2->material=&green;
    sphere3->material = &red;
    sphere4->material = &red;
    pm->material = &white;
    plane->material = &green;
    plane2->material= &no_trans_white;
    plane->material = &no_trans_white;
    // light
    // 
    Vertex light_location = Vertex(0.0,1.2f,4.0f);//(0.0f,5,-1); 


    // create full scene (need to add objects):
    Scene *my_scene = new Scene();
    my_scene->object_list = sphere;


    // my emit_photon function updates photon_map
    // With that photon map we build a tree. 

    int num_photons_caustic=100000;
    int num_photons = 100000;
    float power_ratio = (float)num_photons/(float)num_photons_caustic;
    my_scene->emit_photon_caustic(light_location,num_photons_caustic);
    std::cout << "caustic map built with size: " << endl;
    std::cout << my_scene->caustic_map.pts.size() << endl;
   
    
    my_scene->emit_photon_diffuse_point(light_location,num_photons); 
    std::cout << "photon map built with size: " << endl;
    std::cout << my_scene->map.pts.size() << endl;

    // building trees
	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<float, PointCloud > ,
		PointCloud,
		3 /* dim */
		> my_kd_tree_t;
    
    my_kd_tree_t caustic_index(3 /*dim*/, my_scene->caustic_map, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    caustic_index.buildIndex();
    my_scene->caustic_tree = &caustic_index;
	my_kd_tree_t  index(3 /*dim*/, my_scene->map, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	index.buildIndex();
    my_scene->tree = &index;


    // We have a tree now called index.
    // This tree stores photons.
    // query_point should be hit.position when raytracing, number closest determines how many photons we want to include.
    // out_indices gives indices for map.


    // Creating a camera.
    Vector up,look,eye;
    up = Vector(0,-1,0);
    look = Vector(0.01f,0.01f,8.0f);
    eye = Vector(0.01f,1,-10); // eye is a vector to allow vector operations.
    float distance = width*2; // 2000 for 1024
    float scale = 1;

    Ray ray;
    ray.position.x = eye.x;
    ray.position.y = eye.y;
    ray.position.z = eye.z;
    Vector u,my_v;
    Vector w = eye-look;
    w.normalise();
    up.cross(w,u);
    u.normalise();
    w.cross(u,my_v);
    // level
    int level = 5;

    for(int x=0;x<width;x++)
    {
        std::cout << x << " ";
        for(int y=0;y<height;y++)
        {
            // camera stuff
            float fx = scale*((float)x-width/2.0f);
            float fy = scale*((float)y-height/2.0f);
        


            ray.direction.x = fx*u.x+fy*my_v.x-distance*w.x; 
            ray.direction.y = fx*u.y+fy*my_v.y-distance*w.y; 
            ray.direction.z = fx*u.z+fy*my_v.z-distance*w.z;
            ray.direction.normalise();         
            
            Colour colour,caustic_col,global_col,direct_col;
            bool caustic = true;
            my_scene->photon_raytrace(ray,Vertex(eye.x,eye.y,eye.z),colour,light_location,level,caustic,caustic_col,direct_col,global_col);
            
            fb->plotPixel(x,y,colour.r,colour.g,colour.b);
            fb1->plotPixel(x,y,caustic_col.r,caustic_col.g,caustic_col.b);
            fb2->plotPixel(x,y,global_col.r,global_col.g,global_col.b);
            fb3->plotPixel(x,y,direct_col.r,direct_col.g,direct_col.b);
        }
    }
    std::cout << "outputting photon_test.ppm" << endl;
    fb->writeRGBFile((char *)"photon_test.ppm"); // full image
    fb1->writeRGBFile((char*)"caustic_test.ppm"); // caustic and global contribution
    fb2->writeRGBFile((char*)"global_map_test.ppm"); // global photon map contribution
    fb3->writeRGBFile((char*)"direct_test.ppm"); // direct lighting contribution
    return 0;



}

