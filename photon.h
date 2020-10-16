#pragma once

#include "vector.h"
#include "ray.h"
#include "colour.h"

//#include "nanoflann-master\\include\\nanoflann.hpp"

class Photon {
public:
    Ray photon_ray;
    Colour power;
    //bool kdflag;
    //void emit_photon_diffuse_point(Scene *scene,Vertex light_location);
    //void photon_trace(Photon photon, Hit &hit,Scene *scene,int level);
    //void hemisphere_diffuse(Vector normal,Vector &new_direction);
};