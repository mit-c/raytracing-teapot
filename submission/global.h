/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

// Phong is a child class of Material and implement the simple Phong
// surface illumination model.

#pragma once

#include "material.h"
#include "colour.h"
#include "vector.h"
#include "vertex.h"
#include "ray.h"
#include "hit.h"

class Global : public Material {
public:
    Colour ambient;
	Colour diffuse;
	Colour specular;
    

	float  power; // maybe needed
    

    Ray compute_reflection_ray(Vector &viewer, Vector &normal, Vector &ldir, Hit best_hit);
	void compute_base_colour(Colour &result);
	void compute_light_colour(Vector &viewer, Vector &normal,Vector &ldir, Colour &result);
    void compute_reflection(Hit reflect_hit, Colour &colour);
    bool is_reflective()
    {
        return true;
    }
};
