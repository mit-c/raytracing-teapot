/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2019.
 *
 * Do what you like with this code as long as you retain this comment.
 */

#pragma once

#include "vertex.h"
#include "transform.h"
#include "object.h"

typedef int TriangleIndex[3];

class PolyMesh:public Object{
public:
	int vertex_count;
	int triangle_count;
        Vertex *vertex;
	Vector *face_normal;
	Vector *vertex_normal;
	TriangleIndex *triangle;
	bool flip_normals;
	bool interpolate;
	// add a cheaky thing here saying if you want normals interpolated or not.
	void do_construct(char *file, Transform *transform);
	float test_edge(Vector &normal, Vertex &p, Vertex &v1, Vertex &v0);
	void triangle_intersection(Ray ray, Hit &hit, int which_triangle);
	void intersection(Ray ray, Hit &hit);
	void compute_face_normal(int which_triangle, Vector &normal);
	void compute_vertex_normals(void);
	bool rayTriangleIntersect(const Ray& ray, const Vector &v0, const Vector &v1, const Vector &v2, float &t,float &u, float &v);
	PolyMesh(char *file, bool interpolate_on, bool flip);
	PolyMesh(char *file, Transform *transform, bool interpolate_on, bool  flip);
	~PolyMesh(){}
};
