#include "plane.h"

Plane::Plane(Vector N, Vector point) // ax + by + cz + d = 0, normal is (a,b,c)
{
    Normal = N;
    p0 = point;
}

void Plane::intersection(Ray ray, Hit &hit)
{
    Vector raydir = ray.direction;
    Vector raypos = Vector(ray.position.x,ray.position.y,ray.position.z);
    if(Normal.dot(raydir)==0)
    {
        hit.flag = false;
        return;
    }

    float t = -Normal.dot(p0 - raypos)/(Normal.dot(raydir));

    if(t<0)
    {
        hit.flag = false;
    }

    hit.flag = true;
    Vector hitpos = raypos + t*raydir;
    hit.position = Vertex(hitpos.x,hitpos.y,hitpos.z);
    hit.what = this;
    hit.normal = Normal;

}