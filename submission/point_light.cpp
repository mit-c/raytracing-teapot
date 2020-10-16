/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#include "point_light.h"

PointLight::PointLight()
{
	Light();
}

PointLight::PointLight(Vertex pos, Colour col)
{
	Light();

	position = pos;
	intensity = col;
}

bool PointLight::get_direction(Vertex &surface, Vector &dir)
{
	dir.x = - position.x + surface.x;
    dir.y = - position.y + surface.y;
    dir.z = - position.z + surface.z;
    dir.normalise();
	return true;
}

void PointLight::get_intensity(Vertex &surface, Colour &level)
{
	level = intensity;
}
