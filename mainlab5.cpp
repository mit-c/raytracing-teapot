#include "nanoflann-master/include/nanoflann.hpp"
//#include "nanoflann-master\\utils.h"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <algorithm>
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


using namespace std;
using namespace nanoflann;


// photon mapping lab

// this creates structure for storing photon to use with nanoflann.
// to fix code: divide by 2 whenever you use phong and also divide photons by total emitted.

int main()
{
    int pixel_power = 3; // 3 for desktop resolution

    int width = 240*pow(2,pixel_power); // 1920 x 1080 for my background size (pixel_power = 3).
    int height = 135*pow(2,pixel_power);
    std::cout << "creating image with width: " << width << " and height: " << height << endl;
    FrameBuffer *fb = new FrameBuffer(width,height);
    FrameBuffer *fb1 = new FrameBuffer(width,height);
    FrameBuffer *fb2 = new FrameBuffer(width,height);
    FrameBuffer *fb3 = new FrameBuffer(width,height);

    Sphere *sphere  = new Sphere(Vertex(2.5,1.8,5.0), 1.0f);
    Sphere *sphere2 = new Sphere(Vertex(0.0,1.8,5.0),1.0f);
    Sphere *sphere3 = new Sphere(Vertex(-2.5,1.8,5.0),1.0f); 
    Sphere *sphere4 = new Sphere(Vertex(0.0,-1.5,4.6),0.5f);// this sphere is inside teapot
    // (x,y,depth)
    Sphere *sphere5 = new Sphere(Vertex(0.0,-1.5,4.6),1.0f);
    
    
    // TEAPOT STUFF
    
    // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
    Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
                                        0.0f, 0.0f, 1.0f, -2.7f,
                                        0.0f, 1.0f, 0.0f, 5.0f,
                                        0.0f, 0.0f, 0.0f, 1.0f); 
    Transform *f_16_transform = new Transform(0.5f, 0.0f, 0.0f,  -3.5f,
                                            0.0f, 0.5f, 0.0f, 0.0f,
                                            0.0f, 0.0f, 0.5f, 0.0f,
                                            0.0f, 0.0f, 0.0f, 1.0f); 
    //  Read in the teapot model.
    PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);
    PolyMesh *plane = new PolyMesh((char *)"plane.ply");
    PolyMesh *plane2 = new PolyMesh((char *)"plane2.ply");
    PolyMesh *plane3 = new PolyMesh((char *)"plane3.ply");
    PolyMesh *plane4 = new PolyMesh((char *)"plane4.ply");
    PolyMesh *plane5 = new PolyMesh((char *)"plane5.ply");
    PolyMesh *plane6 = new PolyMesh((char *)"plane6.ply");
    PolyMesh *f16 = new PolyMesh((char *)"f16.ply", f_16_transform);
    sphere->next = sphere2;
    sphere2->next = sphere3;
    sphere3->next = sphere4;
    sphere4->next = pm;
    f16->next = plane;
    pm->next = plane;
    plane->next = plane2; // I think the problem with the planes is the direct lighting is messing up.
    //plane2->next = plane3; // plane 3 is roof I think.
    //plane3 -> next = plane4;
    //plane4 -> next = plane5;
    //plane5 -> next = plane6;
 
    

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

    Phong gold;
    // Don't need ambient for photon mapping
 
    gold.diffuse_photon = Colour(1,0.843,0,0);
    gold.diffuse = gold.diffuse_photon;
    gold.specular_photon = Colour(1,0.843,0,0);
    gold.specular = gold.specular_photon;
    gold.power = 200;
    gold.index_refrac = 1.01;
    gold.reflect = 1.0f;
    gold.refract = 0.0f;


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
    no_trans_white.diffuse_photon = Colour(1,1,1,0.0f);
    no_trans_white.diffuse = no_trans_white.diffuse_photon;
    no_trans_white.specular_photon = Colour(0.0f,0.0f,0.0f,0.0f);
    no_trans_white.specular = no_trans_white.specular_photon;
    no_trans_white.power = 40;
    no_trans_white.index_refrac = 1.0f;
    no_trans_white.reflect = 0.4f; // 0.9 a bit too bright
    no_trans_white.refract = 0.0f;

    Phong sky_blue;
    sky_blue.diffuse_photon = Colour(0.59f,0.808f,0.922,0.0f);
    sky_blue.diffuse = sky_blue.diffuse_photon;
    sky_blue.specular_photon =  Colour(0.59f,0.808f,0.922,0.0f);
    sky_blue.specular =  Colour(0.59f,0.808f,0.922,0.0f);;
    sky_blue.power = 40.0f;
    sky_blue.reflect = 0.0f;
    sky_blue.refract = 0.0f;


    /*
    left->material = &green;
    ceiling->material = &red;
    floor->material = &blue;
    back->material = &red;
    right->material = &blue;
    */
    f16->material=&white;
    sphere->material=&blue;
    sphere2->material=&green;
    sphere3->material = &red;
    sphere4->material = &red;
    pm->material = &white;// &white;
    plane->material = &green;
    plane2->material= &no_trans_white;
    plane->material = &no_trans_white;
    plane3->material = &no_trans_white;
    plane4->material = &no_trans_white;
    plane5->material = &no_trans_white;
    plane6->material = &no_trans_white;
    sphere5->material = &white;

    sphere5->next=plane;
    // light
    // 
    Vertex light_location = Vertex(0,5,-1);//(0.0f,4,-1); (0,1.2,4) 0.0,-1.5,4.6)


    // create full scene (need to add objects):
    Scene *my_scene = new Scene();
    my_scene->object_list = sphere;


    // my emit_photon function updates photon_map
    // With that photon map we build a tree. 

    int num_photons_caustic=300000;
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

    // In order to take out need to build photon map every time need to store 



    // We have a tree now called index.
    // This tree stores photons.
    // query_point should be hit.position when raytracing, number closest determines how many photons we want to include.
    // out_indices gives indices for map.


    // Creating a camera. I think y=height z=depth for eye = Vector(x,y,z).
    Vector up,look,eye;
    up = Vector(0,-1,0);
    look = Vector(0.0,-1.5,7.0f);//(0.01f,0.01f,8.0f);
    eye = Vector(0.01,-2,-12);//(0.01f,-2,-12); // (0.01,1,-10). location of sphere inside teapot 0.0,-1.5,4.6 
    float distance = height*1.9; // was 2
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
    fb1->writeRGBFile((char*)"caustic_test.ppm"); // caustic contribution
    fb2->writeRGBFile((char*)"global_map_test.ppm"); // global photon map contribution
    fb3->writeRGBFile((char*)"direct_test.ppm"); // direct lighting contribution reduce by a factor of 5 or 4
    system("magick convert photon_test.ppm photon_test4.png");
    system("magick convert caustic_test.ppm caustic_test4.png");
    system("magick convert global_map_test.ppm global_map_test4.png");
    system("magick convert direct_test.ppm direct_test4.png");
    return 0;



}

