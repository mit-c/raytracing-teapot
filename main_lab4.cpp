/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

/* This is the entry point function for the program you need to create for lab two.
 * You should not need to modify this code.
 * It creates a framebuffer, loads an triangle mesh object, calls the drawing function to render the object and then outputs the framebuffer as a ppm file.
 *
 * On linux.bath.ac.uk:
 *
 * Compile the code using g++ -o lab2executable main_lab2.cpp framebuffer.cpp linedrawer.cpp polymesh.cpp -lm
 *
 * Execute the code using ./lab2executable
 *
 * This will produce an image file called test.ppm. You can convert this a png file for viewing using
 *
 * pbmropng test.ppm > test.png
 *
 * You are expected to fill in the missing code in polymesh.cpp.
 */
#include "scene.h"
#include "global.h"
#include "framebuffer.h"
#include "ray.h"
#include "hit.h"
#include "polymesh.h"
#include "sphere.h"
#include "light.h"
#include "directional_light.h"
#include "material.h"
#include "phong.h"
#include "point_light.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using namespace std;
/*
void object_test(Ray ray, Object *objects, Hit &best_hit)
{
  Object *obj = objects;

  best_hit.flag = false;


  while(obj != 0)
  {
    Hit obj_hit;
    obj_hit.flag=false;
	  
    obj->intersection(ray, obj_hit,outside);

    
    if (obj_hit.flag)
    {
      if (obj_hit.t > 0.0f)
      {
        if (best_hit.flag == false)
      {
          best_hit = obj_hit;
      } else if (obj_hit.t < best_hit.t)
      {
        best_hit = obj_hit;
      }
      }
    }
    
    obj = obj->next;
  }


  return;
}


void raytrace(Ray ray, Object *objects, Light *lights, Colour &colour, float &depth)
{
  // first step, find the closest primitive

  Hit shadow_hit;
  Hit best_hit;
  object_test(ray, objects, best_hit);

  
  // if we found a primitive then compute the colour we should see
  if(best_hit.flag)
  {
    best_hit.what->material->compute_base_colour(colour);
    depth = best_hit.t;
    Light *light = lights;

    while (light != (Light *)0)
    {
      Vector viewer;
      Vector ldir;

      viewer.x = -best_hit.position.x;
      viewer.y = -best_hit.position.y;
      viewer.z = -best_hit.position.z;
      viewer.normalise();

      bool lit;
      lit = light->get_direction(best_hit.position, ldir);

      if(ldir.dot(best_hit.normal)>0)
      {
	      lit=false;//light is facing wrong way.
      }

      if(lit)
      {
            
        Ray shadow_ray;

        shadow_ray.direction.x = -ldir.x;
        shadow_ray.direction.y = -ldir.y;
        shadow_ray.direction.z = -ldir.z;
        shadow_ray.position.x = best_hit.position.x + (0.0001f * shadow_ray.direction.x);
        shadow_ray.position.y = best_hit.position.y + (0.0001f * shadow_ray.direction.y);
        shadow_ray.position.z = best_hit.position.z + (0.0001f * shadow_ray.direction.z);



        object_test(shadow_ray, objects, shadow_hit);

        if(shadow_hit.flag==true)
        {
          if (shadow_hit.t < 1000000000.0f)
          {
            lit = false; //there's a shadow so no lighting, if realistically close
          }
        }
      }

      if (lit)
      {
        Colour intensity;
        Colour scaling;

        light->get_intensity(best_hit.position, scaling);

        best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

        intensity.scale(scaling);

        colour.add(intensity);
      }

      light = light->next;
    }

    // TODO: compute reflection ray if material supports it.
    if(1) // if best_hit.what->material is global
    {

    }

    // TODO: compute refraction ray if material supports it.
    if(1)
    {
    }

  } else
  {
    depth = 7.0f;
    colour.r = 0.0f;
    colour.g = 0.0f;
    colour.b = 0.0f;
  }	
}
*/

int main(int argc, char *argv[])
{
  int width = 2048;
  int height = 2048;
  // Create a framebuffer
  FrameBuffer *fb = new FrameBuffer(width,height);

  // The following transform allows 4D homogeneous coordinates to be transformed. It moves the supplied teapot model to somewhere visible.
  Transform *transform = new Transform(1.0f, 0.0f, 0.0f,  0.0f,
				                               0.0f, 0.0f, 1.0f, -2.7f,
                                       0.0f, 1.0f, 0.0f, 5.0f,
                                       0.0f, 0.0f, 0.0f, 1.0f); 

  //  Read in the teapot model.
  PolyMesh *pm = new PolyMesh((char *)"teapot_smaller.ply", transform);
  PolyMesh *plane = new PolyMesh((char *)"plane.ply"); 
  //PolyMesh *plane2 = new PolyMesh((char*)"plane2.ply");
  //PolyMesh *plane3 = new PolyMesh((char*)"plane3.ply");
  //PolyMesh *plane4 = new PolyMesh((char*)"plane4.ply");
  //PolyMesh *plane5 = new PolyMesh((char*)"plane5.ply");
  //PolyMesh *plane6 = new PolyMesh((char*)"plane6.ply");
  
  Sphere *sphere  = new Sphere(Vertex(2.5,1.2,8.0), 1.0f);
  Sphere *sphere2 = new Sphere(Vertex(0.0,1.2,8.0),1.0f);
  Sphere *sphere3 = new Sphere(Vertex(-2.5,1.2,8.0),1.0f); // y 1.2
  Sphere *sphere4 = new Sphere(Vertex(0.0,-1.5,4),0.6f);//(Vertex(1.75,-0.5,-2.0),1.0f);
  //sphere->next = pm;

  Ray ray;
  
 

  // light stuff
  DirectionalLight *dl = new DirectionalLight(Vector(1.01f, -1.0f, 1.0f),Colour(1.0f, 1.0f, 1.0f, 0.0f));
  PointLight *pl = new PointLight(Vertex(0.1,20,20),Colour(1.0f, 1.0f, 1.0f, 0.0f));
  Phong bp1;

	bp1.ambient.r = 0.6f;
	bp1.ambient.g = 0.0f;
	bp1.ambient.b = 0.0f;
	bp1.diffuse.r = 0.6f;
	bp1.diffuse.g = 0.0f;
	bp1.diffuse.b = 0.0f;
	bp1.specular.r = 0.6f;
	bp1.specular.g = 0.6f;
	bp1.specular.b = 0.6f;
	bp1.power = 40.0f;
  bp1.reflect = 1.0f;
  bp1.refract = 0.0f;
  bp1.index_refrac = 1.37f;
	
  pm->material = &bp1;

  Phong bp2; // green sphere

  bp2.ambient.r = 0.0f;
	bp2.ambient.g = 0.4f;
	bp2.ambient.b = 0.0f;
	bp2.diffuse.r = 0.0f;
	bp2.diffuse.g = 0.6f;
	bp2.diffuse.b = 0.0f;
	bp2.specular.r = 0.6f;
	bp2.specular.g = 0.6f;
	bp2.specular.b = 0.6f;
	bp2.power = 40.0f;
  bp2.reflect = 1.0f;
  bp2.refract = 0.0f;
  bp2.index_refrac = 1.37f;

  Phong bp3; // white reflective sphere

  bp3.ambient.r = 1.0f;
  bp3.ambient.g = 1.0f;
  bp3.ambient.b = 1.0f;
  bp3.diffuse.r = 1.0f;
  bp3.diffuse.g = 1.0f;
  bp3.diffuse.b = 1.0f;
  bp3.specular.r = 1.0f;
  bp3.specular.g = 1.0f;
  bp3.specular.b = 1.0f;
  bp3.power = 40.0f;
  bp3.reflect=1.0f;
  bp3.refract=0.0f;
  bp3.index_refrac = 1.1f;

  Phong bp4; // blue sphere
  bp4.ambient.r = 0.0f;
	bp4.ambient.g = 0.0f;
	bp4.ambient.b = 0.6f;
	bp4.diffuse.r = 0.0f;
	bp4.diffuse.g = 0.0f;
	bp4.diffuse.b = 0.6f;
	bp4.specular.r = 0.6f;
	bp4.specular.g = 0.6f;
	bp4.specular.b = 0.6f;
	bp4.power = 40.0f;
  bp4.reflect = 1.0f;
  bp4.refract = 1 - bp4.reflect;
  bp4.index_refrac = 2.37f;

  Phong bp5; //  white refractive
  bp5.ambient = Colour(0.4f,0.4f,0.4f,0.0f);
  bp5.diffuse = Colour(0.4f,0.4f,0.4f,0.0f);
  bp5.specular = Colour(0.6f,0.6f,0.6f,0.0f);
  bp5.power = 40.0f;
  bp5.reflect = 1.0f;
  bp5.refract = 1.0f;
  bp5.index_refrac = 1.1;



  pm->material = &bp5;
	sphere->material = &bp4;
  sphere2->material = &bp2;
  sphere3->material = &bp1;
  sphere4->material = &bp1;
  plane->material = &bp3;
  //plane2->material = &bp5;
  /*
  plane3->material = &bp5;
  plane4->material = &bp5;
  plane5->material = &bp5;
  plane6->material = &bp5;
  */
	pm->next = sphere2;
  sphere2->next = sphere;
  sphere->next = sphere3;
  sphere3->next = sphere4;
  sphere4->next = plane;
  //plane->next = plane2;
  //plane2->next = plane3;
  //plane3->next = plane4;
  //plane4->next = plane5;
  //plane5->next = plane6;
  Scene my_scene;
  
   // start outside object
  // OBJECTS AND LIGHTS
  my_scene.object_list = pm; // pm,sphere2,sphere,sphere3 is the list
  my_scene.light_list = dl; // only have 1 light at the moment
 
  // camera stuff  
  Vector up,look,eye;
  up = Vector(0,-1,0);
  look = Vector(0.01f,0.01f,8.0f);
  eye = Vector(0.01f,1,-10); //3,1,-10
  float distance = width*2; // 2000 for 1024
  float scale = 1;
  
  // ignore below
  ray.position.x = eye.x;
  ray.position.y = eye.y;
  ray.position.z = eye.z;
  Vector u,my_v;
  Vector w = eye-look;
  w.normalise();
  up.cross(w,u);
  u.normalise();
  w.cross(u,my_v);
  //

  for (int y = 0; y < height; y += 1)
  {
    for (int x = 0; x < width; x += 1)
    {
      float fx = scale*((float)x-width/2.0f);
      float fy = scale*((float)y-height/2.0f);
 

      //Vector direction;

      ray.direction.x = fx*u.x+fy*my_v.x-distance*w.x; //= (fx-0.5f); // ray points at (x y 0.5) where x,y is eye
      ray.direction.y = fx*u.y+fy*my_v.y-distance*w.y; //= (0.5f-fy);
      ray.direction.z = fx*u.z+fy*my_v.z-distance*w.z;//= 0.5f;
      ray.direction.normalise();

      Colour colour;
      float depth;

      my_scene.raytrace(ray,5, my_scene.object_list, my_scene.light_list, colour);

      fb->plotPixel(x, y, colour.r, colour.g, colour.b);
      fb->plotDepth(x,y, depth);
    }

    cerr << y << " " << flush;
  }
  
  // Output the framebuffer.
  fb->writeRGBFile((char *)"test.ppm");
  fb->writeDepthFile((char *)"depth.ppm");
  return 0;
  
}
