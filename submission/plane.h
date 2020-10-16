/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "object.h"

class Plane : public Object {
	Vector Normal;
    Vector p0;
public:
	Plane(Vector normal,Vector p0);
	void intersection(Ray ray, Hit &hit);
};
