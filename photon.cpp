// #include "nanoflann-master\\include\\nanoflann.hpp"
// #include "math.h"
// #include "vector.h"
// #include "vertex.h"
// #include "photon.h"
// #include "scene.h"
// #include "hit.h"
// #include "material.h"


// #include <algorithm>

// using namespace std;
// using namespace nanoflann;




// //my_kd_tree_t index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /*max leaf*/) );
// // need to actually store stuff in this tree then reindex. see nanonflann.

// //PointCloud<float> cloud;

// // the below function finds random direction of emission in hemisphere -- we will most likely need to do something more complicated for specular.
// void Photon::hemisphere_diffuse(Vector normal, Vector &photon_direction) 
// {
//     // generating hemisphere coordinates
//     float x,y,z;
//     Vector N_t = Vector(normal.z,0.0f,-normal.x);
//     Vector N_b;
//     normal.cross(N_t,N_b);
//     N_t.normalise();
//     N_b.normalise();

//     x = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults N_t
//     y = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults normal
//     z = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults N_b
//     y = abs(y); // we want hemisphere.
//     Vector photon_direction = x*N_t + y*normal + z*N_b;
//     photon_direction.normalise();
    



// }

// void Photon::emit_photon_diffuse_point(Scene *scene,Vertex light_location)
// {
//     float x,y,z;
//     int num_of_photons = 1000;
//     for(int i=0;i<num_of_photons;i++)
//     {
//             x = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
//             y = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
//             z = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
//             int level = 5;
//             Photon photon;
//             photon.photon_ray.position = light_location;
//             photon.photon_ray.direction = Vector(x,y,z);
//             photon.photon_ray.direction.normalise();
//             photon.power = Colour(1.0,1.0f,1.0f,0); // light colour
//             Hit best_photon_hit;
//             photon_trace(photon,best_photon_hit,scene,level); // so russian roulette and shadow photons
//             // extract material from best_photon_hit -- from this we create probabilities for new rays.
//             // not sure if this next bit should be in photon trace.


//     }
// }


// void Photon::photon_trace(Photon photon, Hit &hit, Scene *scene, int level)
// {
//     Ray ray = photon.photon_ray;
//     level = level - 1;
//     if(level < 0)
//     {
//         return;
//     }


//     // Need to trace ray through scene
//     // First need to find first primative
//     scene->trace(ray,scene->object_list,hit); // this gives us closest hit for our photon.
//     // We have closest hit now we need to figure out probabilities.
//     if(ray.direction.dot(hit.normal)>0) //makes sure normal is pointing correct direction.
//     {
//         hit.normal.negate();
//     }

    
//     if(photon.power.r == 0 && photon.power.g == 0 && photon.power.b == 0)
//     {
//         // this is a shadow photon. Not sure if I need to use this but might be useful.
//     }

//     // storage here


//     // storage end.
    
//     Colour diffuse = hit.what->material->diffuse_photon;
//     Colour specular = hit.what->material->specular_photon;
//     float prob_reflect = max({diffuse.r + specular.r, diffuse.g + specular.g,diffuse.b + specular.b});
//     float sum_diffuse = diffuse.r + diffuse.g + diffuse.b;
//     float sum_specular = specular.r + specular.g + specular.b;
//     float prob_diffuse = (sum_diffuse/(sum_diffuse + sum_specular))*prob_reflect;
//     float prob_specular = prob_reflect - prob_diffuse;    
//     float xi = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
//     if(xi>=0 && xi<=prob_diffuse)
//     {
//         // diffuse reflection
//         // store photon
//         Colour reflect_power;
//         Photon reflect_photon;
//         // implementing lambertian diffusion. For this type of diffusion the power is constant.
//         Vector L = ray.direction;
//         L.negate();
//         reflect_power.r = photon.power.r*L.dot(hit.normal)*diffuse.r;
//         reflect_power.g = photon.power.g*L.dot(hit.normal)*diffuse.g;
//         reflect_power.b = photon.power.b*L.dot(hit.normal)*diffuse.b;


//         /*
//         reflect_power.r = (power.r*diffuse.r)/prob_diffuse;
//         reflect_power.g = (power.g*diffuse.g)/prob_diffuse;  
//         reflect_power.b = (power.b*diffuse.b)/prob_diffuse;
//         */   
//         // if diffuse this direction should be random in a hemisphere around normal.
//         // create new photon

//         Vector reflect_ray;
//         hemisphere_diffuse(hit.normal, reflect_ray);

//         reflect_photon.photon_ray.direction = reflect_ray;
//         reflect_photon.photon_ray.position.x = hit.position.x + 0.0001*photon_ray.direction.x;
//         reflect_photon.photon_ray.position.y = hit.position.y + 0.0001*photon_ray.direction.y;
//         reflect_photon.photon_ray.position.z = hit.position.z + 0.0001*photon_ray.direction.z;
//         reflect_photon.photon_ray.direction.normalise();
//         reflect_photon.power = reflect_power; // update photon power.
//         Hit reflect_hit;
//         photon_trace(reflect_photon,reflect_hit,scene,level);


//     }
//     else if(xi>prob_diffuse && xi <= prob_specular + prob_diffuse)
//     {
//         // specular
//         // store photon?
//         Colour reflect_power;
//         Photon reflect_photon;
//         reflect_power.r = (power.r*specular.r)/prob_specular;
//         reflect_power.g = (power.g*specular.g)/prob_specular;
//         reflect_power.b = (power.b*specular.b)/prob_specular;
//         reflect_photon.photon_ray.direction = ray.direction - 2.0f*(hit.normal.dot(ray.direction))*hit.normal;
//         reflect_photon.photon_ray.position.x = hit.position.x + 0.0001*photon_ray.direction.x;
//         reflect_photon.photon_ray.position.y = hit.position.y + 0.0001*photon_ray.direction.y;
//         reflect_photon.photon_ray.position.z = hit.position.z + 0.0001*photon_ray.direction.z;
//         reflect_photon.photon_ray.direction.normalise();
//         reflect_photon.power = reflect_power;
//         Hit reflect_hit;
//         photon_trace(reflect_photon,reflect_hit,scene,level);

//     }
//     else
//     {
//         // absorbtion
//     }
//     // for shadow photon 
//     Photon shadow_photon;
//     Hit shadow_hit;
//     shadow_photon.photon_ray.position = hit.position;
//     shadow_photon.photon_ray.direction = ray.direction;
//     shadow_photon.power = Colour(0,0,0,0); // store 0 power (shadow).
//     photon_trace(shadow_photon,shadow_hit,scene,1); // only need to do one level for shadow photon- hence level is 1.
// }

