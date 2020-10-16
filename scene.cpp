/***************************************************************************
 *
 * krt - Kens Raytracer - Coursework Edition. (C) Copyright 1997-2018.
 *
 * Do what you like with this code as long as you retain this comment.
 */
#define _USE_MATH_DEFINES
#include <math.h>
#include "scene.h"
#include "vector.h"
#include "vertex.h"
#include "global.h"
#include <algorithm>
#include <complex>
#include "photon.h"
#include <cmath>

#include "nanoflann-master/include/nanoflann.hpp"
//#include "nanoflann-master\\utils.h"


using namespace std;	


Scene::Scene()
{
	object_list = 0;
	light_list = 0;
	
	// std::vector<PointCloud::Point> pts = std::vector<PointCloud::Point>();
	//map.pts = pts;
}


void Scene::trace(Ray ray, Object *objects, Hit &hit)
{
	Hit current_hit;

	hit.flag = false;
	hit.t = 0.0f;
	hit.what = 0;

	while (objects != 0)
	{
		Hit hit_current;

		objects->intersection(ray, hit_current);

		if (hit_current.flag == true)
		{
			if (hit.flag == false)
			{
				hit = hit_current;

			} else if (hit_current.t < hit.t)
			{
				hit = hit_current;
			}
		}


		objects = objects->next;
	}
}

void Scene::raytrace(Ray ray, int level, Object *objects, Light *lights, Colour &colour)
{

	// a default colour if we hit nothing.
	float red   = 0.2f;
	float green = 0.2f;
	float blue  = 0.2f;

	// check we've not recursed too far.
	level = level - 1;
	if (level < 0)
	{
		return;
	}

	// we might need to check (if we are inside an object, if we are then we need to only look at negated normals for polymesh -- spheres should be fine
	// so we could add bool inside_object to raytrace.
	// Then add inside_object to trace.
	// If inside object is true we negate all the normals of our objects.
	// This would involve adding inside_object to all our objects.
	// how do we know if we are in an object though? -- if we are outside an object and we refract we are inside
	// if we are inside and refract we are outside so in refraction part do opposite bool of previous
	// start outside.
	// now we know if we are in an object how do we implement.
	// change all intersection functions to take bool inside_object.
	// if false do function as before.
	// if true do N.negate() in polymesh/hit.normal.negate() in sphere;

	// first step, find the closest primitive

	Hit best_hit;
	trace(ray, objects, best_hit);
	


	// if we found a primitive then compute the colour we should see
	if (best_hit.flag)
	{
		// if normals are the wrong way round need to swap before reflections
		float material_eta = best_hit.what->material->index_refrac;
		float outside_eta = 1.0f;
		Vector N = best_hit.normal;
		Vector I = ray.direction;
		float NdotI = N.dot(I);

		
		if(I.dot(N)<0) // this does not support polymesh as polymesh does this already.
		{
			NdotI = - NdotI;
		}
		else
		{
			N.negate();
			std::swap(material_eta,outside_eta);
		}
		float eta = outside_eta/material_eta; // so 1/eta = material/outside
		float costhetai = NdotI;
		float costhetat = sqrt(1.0f-(eta*eta)*(1.0f-pow(costhetai,2)));
		float k = 1 - eta*eta*(1-costhetai*costhetai);


		best_hit.what->material->compute_base_colour(colour);

		Light *light = lights;

		while (light != (Light *)0)
		{
			Vector viewer;
			Vector ldir;

			viewer.x = -best_hit.position.x; // This seems really wrong. Should be eye-besthit
			viewer.y = -best_hit.position.y;
			viewer.z = -best_hit.position.z;
			viewer.normalise();

			bool lit;
			float distance;
			lit = light->get_direction(best_hit.position, ldir, distance);

			if(lit)
			{
				Ray shadow_ray;
				Hit shadow_hit;
				// ldir is direction to light
       			shadow_ray.direction.x = -ldir.x;
        		shadow_ray.direction.y = -ldir.y;
        		shadow_ray.direction.z = -ldir.z;
        		shadow_ray.position.x = best_hit.position.x + (0.001f * shadow_ray.direction.x);
        		shadow_ray.position.y = best_hit.position.y + (0.001f * shadow_ray.direction.y);
        		shadow_ray.position.z = best_hit.position.z + (0.001f * shadow_ray.direction.z);				
				trace(shadow_ray,objects,shadow_hit);

				if(shadow_hit.flag == true) // changed this
				{
					lit = false; // just affects ambient stuff still get light from reflections/refractions
					Vector dir;
					/*
					if(shadow_hit.t > distance) // the case where the light is in front of the blocking object
					{
						std::cout << "shadow_hit.t " << shadow_hit.t << " then distance" << distance << endl;
						lit = true;
					}
					*/
				}
				else
				{


					Colour intensity;
					Colour scaling;

					light->get_intensity(best_hit.position, scaling);

					best_hit.what->material->compute_light_colour(viewer, best_hit.normal, ldir, intensity);

					intensity.scale(scaling);

					colour.add(intensity);
				}
			}

			light = light->next;
		}
		////////////////////////////////////////////////////fresnel -- seems to hard to get good looking materials e.g. reflective materials have complex eta.

		float kr = best_hit.what->material->reflect;
		float kt = best_hit.what->material->refract;
		
		if(k<0.0f) // tir if material reflective
		{
			kt = 0.0f; // no refraction only reflection
			if(kr!=0.0f) // if material reflective make 100% reflective
			{
				kr = 1.0f; 
			}
		}
		else // do fresnel
		{
			if(kt!=0)
			{
				float r_par = (material_eta*costhetai-outside_eta*costhetat)/(material_eta*costhetai+outside_eta*costhetat); //eta2 = materialeta
				float r_per = (outside_eta*costhetai-material_eta*costhetat)/(outside_eta*costhetai+material_eta*costhetat); 
				kr = 0.5f*(r_par*r_par+r_per*r_per);
				kt = 1 - kr;
			}
		}
		
		
		
	
		
		//float kr = 0.5f*(Fr_parallel + Fr_orthog); // kr is almost always 0.
	
		//float kt = 1 - kr; // kt is really high
		
		////////////////////////////////////////////////////end of fresnel start of reflect 
		if(kr!=0) // kr=0 => no reflection
		{
			
			Colour r_colour;
			Ray reflect_ray;

			Vector R = I-2.0f*(N.dot(I))*N; 
			R.normalise();
			reflect_ray.position.x = best_hit.position.x + 0.0001f*N.x;
			reflect_ray.position.y = best_hit.position.y + 0.0001f*N.y;
			reflect_ray.position.z = best_hit.position.z + 0.0001f*N.z;

			
			
			reflect_ray.direction = R;



      		raytrace(reflect_ray, level, objects, lights, r_colour);
			r_colour.r = kr*r_colour.r;
			r_colour.g = kr*r_colour.g;
			r_colour.b = kr*r_colour.b;

			colour.add(r_colour); 
		}

		// TODO: if refract inside material need to flip normal.
		
		
		if(kt!=0)
		{

			Colour refract_colour;

			Vector T;
			T = eta*I-(eta*costhetai-costhetat)*N;//(1/eta)*I-(costhetat-(1/eta)*costhetai)*N;
			T.normalise();
			Ray refract_ray;
			refract_ray.direction = T;
			
			refract_ray.position.x = best_hit.position.x + 0.001*T.x;
			refract_ray.position.y = best_hit.position.y + 0.001*T.y;
			refract_ray.position.z = best_hit.position.z + 0.001*T.z; 
			raytrace(refract_ray, level, objects, lights, refract_colour); // if we are inside object we might need to negate the normals
			
			refract_colour.r = kt*refract_colour.r;
			refract_colour.g = kt*refract_colour.g;
			refract_colour.b = kt*refract_colour.b;

			colour.add(refract_colour);
		}
		
	}
	else
	{
		colour.r = 0.0f;
		colour.g = 0.0f;
		colour.b = 0.0f;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
//START PHOTON MAPPING HERE///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////




// the below function finds random direction of emission in hemisphere -- we will most likely need to do something more complicated for specular.
void Scene::hemisphere_diffuse(Vector normal, Vector &photon_direction) 
{
    // generating hemisphere coordinates


    float x,y,z;
    Vector N_t = Vector(normal.z,0.0f,-normal.x); // normal.z = 0 normal.x = 0 normal.y = 1
    Vector N_b;
    normal.cross(N_t,N_b);
    N_t.normalise();
    N_b.normalise();
	//std::cout << "normal in hemisphere" << normal.x << normal.y << normal.z << endl; 
    x = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults N_t
    y = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults normal
    z = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))); // mults N_b
    y = abs(y); // we want hemisphere.
    photon_direction = x*N_t + y*normal + z*N_b;
    photon_direction.normalise();

}

void Scene::emit_photon_caustic(Vertex light_location, int num_of_photons)
{

    float x,y,z;
    for(int i=0;i<num_of_photons;i++)
    {	
		if((i+1)%1000==0)	
		{
			std::cout << "caustic photon: " << i+1 << endl;
		}
			
		x = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
		y = -abs(-1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)))); // emits light downwards
		z = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
		int level = 6; // level is larger for caustic photons
		Photon photon;
		photon.photon_ray.position = light_location;
		photon.photon_ray.direction = Vector(x,y,z);
		photon.photon_ray.direction.normalise();
		float photon_power =1.0f/450.0f; // in order to get good image you need to have enough power to show the weaker caustics then set a maximum intensity.
		photon.power = Colour(photon_power,photon_power,photon_power,0.0f);
		int diffuse_level = 1;
		int specular_level = 0;
		bool caustic = true;
		photon_trace(photon,level,diffuse_level,specular_level,caustic); // so russian roulette and shadow photons
		// extract material from best_photon_hit -- from this we create probabilities for new rays.
		// not sure if this next bit should be in photon trace.


    }
}

void Scene::emit_photon_diffuse_point(Vertex light_location, int num_of_photons)
{
	
    float x,y,z;
    for(int i=0;i<num_of_photons;i++)
    {	
		if((i+1)%1000==0)	
		{
			std::cout << "global photon: " << i+1 << endl;
		}
		// emit in random direction but downwards.
		x = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
		y = -abs(-1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2))));
		z = -1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(2)));
		int level = 5;
		Photon photon;
		photon.photon_ray.position = light_location;
		photon.photon_ray.direction = Vector(x,y,z);
		photon.photon_ray.direction.normalise();
		float photon_power = 1.0f/450.0f;;//((float)num_of_photons*0.0005); // 200 = 20,000 / 100 = num_of_photons/ 100
		photon.power = Colour(photon_power,photon_power,photon_power,0.0f); // light colour
		int diffuse_level = 1000; // this parameter only affects our photon_trace when it is less than level.
		// diffuse_level is designed to implement caustics.
		bool caustic = false;
		int specular_level = 0; // this doesn't matter as caustic = false
		photon_trace(photon,level,diffuse_level,specular_level, caustic); // so russian roulette and shadow photons
		// extract material from best_photon_hit -- from this we create probabilities for new rays.
		// not sure if this next bit should be in photon trace.


    }
}


void Scene::photon_trace(Photon photon, int level, int diffuse_level,int specular_level, bool caustic)
{


    Ray ray = photon.photon_ray;
	Hit hit;
    level = level - 1;

	
	
    if(level < 0)
    {
        return;
    }
	// for caustics we only want specular
	// having a diffuse photon interaction causes diffuse_level to reduce by 1.
	// for caustics we only want specular reflections to diffuse so we set diffuse_level to be 1.
	// if we have a diffuse interaction our diffuse_level goes to 0 and we stop tracing.
	if(diffuse_level <= 0 && caustic) 
	{
		return; 
	}

	
    // Need to trace ray through scene
    // First need to find first primative
	
    trace(ray,object_list,hit); // this gives us closest hit for our photon.
	if(hit.flag == false)
	{
		return;
	}

	float material_eta = hit.what->material->index_refrac;
	
	float outside_eta = 1.0f;
	float NdotI = ray.direction.dot(hit.normal);

    // We have closest hit now we need to figure out probabilities.
	
    if(NdotI < 0) //makes sure normal is pointing correct direction.
    {
		NdotI = -NdotI;
	}
	else // if this happens we are inside the object and need to swap our index of refractions.
	{
        hit.normal.negate();
		std::swap(material_eta,outside_eta);
    }
	
   
   

    // update photon position for storage.
	photon.photon_ray.position = hit.position;
	
	//std::cout << hit.what->material->diffuse_photon.r << endl;
    
    Colour diffuse = hit.what->material->diffuse_photon;
	
	
    Colour specular = hit.what->material->specular_photon;

    float prob_reflect = max({diffuse.r + specular.r, diffuse.g + specular.g,diffuse.b + specular.b});
    float sum_diffuse = diffuse.r + diffuse.g + diffuse.b;
    float sum_specular = specular.r + specular.g + specular.b;
    float prob_diffuse = (sum_diffuse/(sum_diffuse + sum_specular))*prob_reflect;
    float prob_specular = prob_reflect - prob_diffuse;    
    float xi = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	
    if(xi>=0 && xi<=prob_diffuse)
    {
		if(caustic && (specular_level==0)) // need to have at least one specular bounce.
		{
			return;
		}
		diffuse_level = diffuse_level - 1;
        // diffuse reflection
        // store photon
		//std::cout << photon.power.r << photon.power.g << photon.power.b << endl;
		store_photon(photon,caustic);
        Colour reflect_power;
        Photon reflect_photon;
        // implementing lambertian diffusion. For this type of diffusion the power is constant.
        
		Vector L = ray.direction;
        L.negate();
		float LdotN = max(L.dot(hit.normal),0.0f);
        reflect_power.r = photon.power.r*LdotN*diffuse.r;
        reflect_power.g = photon.power.g*LdotN*diffuse.g;
        reflect_power.b = photon.power.b*LdotN*diffuse.b;
		

        Vector reflect_ray;
        hemisphere_diffuse(hit.normal, reflect_ray); // this function finds a new direction for the photon in a random direction which isn't into the surface.

        reflect_photon.photon_ray.direction = reflect_ray;
        reflect_photon.photon_ray.position.x = hit.position.x + 0.0001*reflect_photon.photon_ray.direction.x;
        reflect_photon.photon_ray.position.y = hit.position.y + 0.0001*reflect_photon.photon_ray.direction.y;
        reflect_photon.photon_ray.position.z = hit.position.z + 0.0001*reflect_photon.photon_ray.direction.z;
        reflect_photon.photon_ray.direction.normalise();
        reflect_photon.power = reflect_power; // update photon power.

		//std::cout << reflect_ray.x << reflect_ray.y << reflect_ray.z << endl;
        
        photon_trace(reflect_photon,level,diffuse_level,specular_level,caustic);


    }
    else if(xi>prob_diffuse && xi <= prob_specular + prob_diffuse)
    {
		specular_level = specular_level + 1;
        // specular (no storage)
		// this should generate caustics but would need to fire photons more directly at object
		// To do caustics we would also need to make sure we only look at certain interactions
		// We can do this by adding levels i.e. counters which change depending on type of interaction.
	

		// Implement reflection and refraction.
		Vector I = ray.direction;
		Vector N = hit.normal;
		float eta = outside_eta/material_eta;
		float costhetai = NdotI;
		float costhetat = sqrt(1.0f-(eta*eta)*(1.0f-pow(costhetai,2)));
		float kr = hit.what->material->reflect;
		float kt = hit.what->material->refract;
		float rand01 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float k = 1 - eta*eta*(1-costhetai*costhetai);
		Colour new_power;
		new_power.r = (photon.power.r*specular.r)/prob_specular;
		new_power.g = (photon.power.g*specular.g)/prob_specular;
		new_power.b = (photon.power.b*specular.b)/prob_specular;

		if(k<0.0f) // tir if material reflective
		{
			kt = 0.0f; // no refraction only reflection
			if(kr!=0.0f) // if material reflective make 100% reflective
			{
				kr = 1.0f; 
			}
		}
		else // do fresnel
		{
			if(kt!=0)
			{
				float r_par = (material_eta*costhetai-outside_eta*costhetat)/(material_eta*costhetai+outside_eta*costhetat); //eta2 = materialeta
				float r_per = (outside_eta*costhetai-material_eta*costhetat)/(outside_eta*costhetai+material_eta*costhetat); 
				kr = 0.5f*(r_par*r_par+r_per*r_per);
				kt = 1 - kr; 
			}
		}
		
		if(rand01 < kr) // reflect
		{
			Colour reflect_power;
			Photon reflect_photon;

			reflect_photon.photon_ray.direction = I - 2.0f*(I.dot(N))*N;
			reflect_photon.photon_ray.position.x = hit.position.x + 0.0001*reflect_photon.photon_ray.direction.x;
			reflect_photon.photon_ray.position.y = hit.position.y + 0.0001*reflect_photon.photon_ray.direction.y;
			reflect_photon.photon_ray.position.z = hit.position.z + 0.0001*reflect_photon.photon_ray.direction.z;
			reflect_photon.photon_ray.direction.normalise();
			reflect_photon.power = photon.power;
			photon_trace(reflect_photon,level,diffuse_level,specular_level,caustic);
	
		}
		else // refract
		{
			Colour refract_power;
			Photon refract_photon;
			
			refract_photon.power = photon.power;// photon.power;
			refract_photon.photon_ray.direction = eta*I-(eta*costhetai-costhetat)*N;
			refract_photon.photon_ray.position.x = hit.position.x + 0.0001*refract_photon.photon_ray.direction.x;
			refract_photon.photon_ray.position.y = hit.position.y + 0.0001*refract_photon.photon_ray.direction.y;
			refract_photon.photon_ray.position.z = hit.position.z + 0.0001*refract_photon.photon_ray.direction.z;
			photon_trace(refract_photon,level,diffuse_level,specular_level,caustic);

		}

    }
    else
    {
		
        // absorbtion (storage but no new photon)
		store_photon(photon,caustic);
		
    }
    // for shadow photon 
	// I decided against using shadow photons.
	
    Photon shadow_photon;
    Hit shadow_hit;
    //shadow_photon.photon_ray.position = hit.position;
	/*
    shadow_photon.photon_ray.direction = ray.direction;
    shadow_photon.power = Colour(0,0,0,0); // store 0 power (shadow).
	trace(shadow_photon.photon_ray,object_list,shadow_hit);
    if(shadow_hit.flag == true) // if our shadow ray intersects we store the photon at the hit position. No shadow rays for refractive
	{	
		shadow_photon.photon_ray.position = shadow_hit.position;
		store_photon(shadow_photon,caustic);
	}
	*/
	
}


void Scene::store_photon(Photon photon,bool caustic)
{

	
	PointCloud::Point point;

	
	point.x = photon.photon_ray.position.x;
	point.y = photon.photon_ray.position.y;
	point.z = photon.photon_ray.position.z;
	point.incident_direction = photon.photon_ray.direction;
	point.power = photon.power;
	if(caustic)
	{
		caustic_map.pts.push_back(point);	
	}
	else
	{
		map.pts.push_back(point);
	}

}


void Scene::photon_raytrace(Ray ray, Vertex eye, Colour &overall_colour, Vertex light_pos,int level,bool caustic,Colour &caustic_output, Colour &direct_output, Colour &global_map_output)
{
	level = level - 1;


	if(level<0)
	{
		return;
	}
	
	Hit hit;
	trace(ray,object_list,hit);
	
	if(hit.flag==false)
	{
		return;
	}

	// we have a hit.

	Ray tolight;
	
	tolight.direction.x = light_pos.x - hit.position.x;
	tolight.direction.y = light_pos.y - hit.position.y;
	tolight.direction.z = light_pos.z - hit.position.z;
	float distance;
	distance = tolight.direction.length();
	tolight.direction.normalise();
	tolight.position.x = hit.position.x + 0.0001*tolight.direction.x;
	tolight.position.y = hit.position.y + 0.0001*tolight.direction.y;
	tolight.position.z = hit.position.z + 0.0001*tolight.direction.z;
	
	Hit direct_hit;
	trace(tolight,object_list,direct_hit);
	if(direct_hit.t > distance)
	{
		//std::cout << "direct_hit.t: " << direct_hit.t << " distance to light: " << distance << endl; 
		direct_hit.flag = false;
	}
	Vector V,N,R,L,I;
	float red,green,blue,RdotV,LdotN;
	// reflection refraction stuff;
	float material_eta = hit.what->material->index_refrac;
	float outside_eta = 1.0f;
	float NdotI = hit.normal.dot(ray.direction);


	if(hit.normal.dot(ray.direction)<0)
	{
		NdotI = - NdotI;
	}
	else
	{
		hit.normal.negate();
		swap(material_eta,outside_eta);
	}
	float eta = outside_eta/material_eta; // so 1/eta = material/outside
	float costhetai = NdotI;
	float costhetat = sqrt(1.0f-(eta*eta)*(1.0f-pow(costhetai,2)));
	float k = 1 - eta*eta*(1-costhetai*costhetai);
	float kr,kt;
	kr = hit.what->material->reflect;
	kt = hit.what->material->refract;
	bool bool_refract = true;
	if(kt != 0) // If kt != 0 then to TIR and fresnel else keep kr the same to allow non-trasparent materials.
	{
		if(k<0.0f & bool_refract) // tir if material reflective
		{
			kt = 0.0f; // no refraction only reflection
			if(kr!=0.0f) // if material reflective make 100% reflective
			{
				kr = 1.0f; 
				
			}
		}
		else // do fresnel
		{
			if(kt!=0)
			{
				float r_par = (material_eta*costhetai-outside_eta*costhetat)/(material_eta*costhetai+outside_eta*costhetat); //eta2 = materialeta
				float r_per = (outside_eta*costhetai-material_eta*costhetat)/(outside_eta*costhetai+material_eta*costhetat); 
				kr = 0.5f*(r_par*r_par+r_per*r_per);
				kt = 1 - kr;
			}
			
		}
	}

	V.x = eye.x - hit.position.x;
	V.y = eye.y - hit.position.y;
	V.z = eye.z - hit.position.z;
	V.normalise();
	N = hit.normal;
	N.normalise();
	I = ray.direction;

	Colour diffuse = hit.what->material->diffuse_photon;
		
	Colour specular = hit.what->material->specular_photon;
	Colour direct = Colour(0,0,0,0);

	if(direct_hit.flag == false)//  if we have no hit we do direct lighting i.e. raytrace normally with no ambient
	{
		Vector eye_to_hit = V;
		eye_to_hit.negate();
		//std::cout<<"we have direct lighting" << endl;
		Vector ldir = tolight.direction;
		ldir.negate();
		hit.what->material->compute_light_colour(eye_to_hit,N,ldir,direct);
		// not sure if reflection and refraction should be outside this if statement but will keep inside at the moment.
		// Basically should you see reflection and refraction in shadow.
		direct.r = min(direct.r,1.0f);
		direct.g = min(direct.g,1.0f);
		direct.b = min(direct.b,1.0f);

	}
	Colour global = Colour(0,0,0,0);

	if(kr!=0) // if our material is reflective.
	{
		Colour r_colour;
		Ray reflect_ray;

		Vector R = I-2.0f*(N.dot(I))*N; 
		R.normalise();
		reflect_ray.position.x = hit.position.x + 0.0001f*N.x;
		reflect_ray.position.y = hit.position.y + 0.0001f*N.y;
		reflect_ray.position.z = hit.position.z + 0.0001f*N.z;

		reflect_ray.direction = R;

		photon_raytrace(reflect_ray,eye,r_colour,light_pos,level,caustic,caustic_output,direct_output,global_map_output);
		r_colour.r = kr*r_colour.r;
		r_colour.g = kr*r_colour.g;
		r_colour.b = kr*r_colour.b;

		global.add(r_colour); 
	}

	if(kt!=0) // if our material is refractive
	{
		Colour refract_colour;

		Vector T;
		T = eta*I-(eta*costhetai-costhetat)*N;
		T.normalise();
		Ray refract_ray;
		refract_ray.direction = T;

		refract_ray.position.x = hit.position.x + 0.0001*T.x;
		refract_ray.position.y = hit.position.y + 0.0001*T.y;
		refract_ray.position.z = hit.position.z + 0.0001*T.z; 

		photon_raytrace(refract_ray,eye,refract_colour,light_pos,level,caustic,direct_output,caustic_output,global_map_output);

		refract_colour.r = kt*refract_colour.r;
		refract_colour.g = kt*refract_colour.g;
		refract_colour.b = kt*refract_colour.b;

		global.add(refract_colour);
	}
	global.r = min(global.r,1.0f);
	global.g = min(global.g,1.0f);
	global.b = min(global.b,1.0f);
	
	const float query_point[3] = {hit.position.x,hit.position.y,hit.position.z};
	
	size_t num_closest = 200; // This should be larger I think
	std::vector<size_t> out_indices(num_closest);
	std::vector<float> out_distances_sq(num_closest);
	
	num_closest = tree->knnSearch(&query_point[0],num_closest,&out_indices[0],&out_distances_sq[0]);
	out_indices.resize(num_closest);
	out_distances_sq.resize(num_closest); // resizing means we don't search where we shouldn't in the for loop.
	float radius_sq = *max_element(out_distances_sq.begin(),out_distances_sq.end());
	float circle_area  = M_PI*radius_sq;

	
	Colour indirect = Colour(0.0f,0.0f,0.0f,0.0f);
	
	for(size_t i=0;i<num_closest;i++) 
	{	
		
		size_t index = out_indices[i];

		float x = map.pts[index].x;
		float y = map.pts[index].y;
		float z = map.pts[index].z;
	
		Vector incident_direction = map.pts[index].incident_direction;

		Colour photon_power = map.pts[index].power;
		int photon_map_size = map.pts.size();

		// BRDF
		
		L = incident_direction;
		L.negate();
		L.normalise();
		
		V.x = eye.x - x;
		V.y = eye.y - y;
		V.z = eye.z - z;
		V.normalise();
		R = incident_direction - 2.0f*(incident_direction.dot(N))*N;
		R.normalise();


		RdotV = max(R.dot(V),0.0f);
		LdotN = max(L.dot(N),0.0f);

		// implementing phong diffuse and specular.
		red 	= min(0.5f*(diffuse.r*(LdotN)*photon_power.r + specular.r*pow(RdotV,20.0f)*photon_power.r),1.0f);
		green 	= min(0.5f*(diffuse.g*(LdotN)*photon_power.g + specular.g*pow(RdotV,20.0f)*photon_power.g),1.0f);
		blue 	= min(0.5f*(diffuse.b*(LdotN)*photon_power.b + specular.b*pow(RdotV,20.0f)*photon_power.b),1.0f);


		indirect.r += red;
		indirect.g += green;
		indirect.b += blue;
		

	}
	indirect.r = min(indirect.r/circle_area,1.0f);
	indirect.g = min(indirect.g/circle_area,1.0f);
	indirect.b = min(indirect.b/circle_area,1.0f);
	
	Colour caustic_col = Colour(0,0,0,0);
	if(caustic) // if we have caustics compute caustic radiance.
	{
		float red_c,green_c,blue_c;
		red_c = 0;
		green_c = 0;
		blue_c = 0;
		size_t num_closest_c = 20; 
		std::vector<size_t> out_indices_c(num_closest_c);
		std::vector<float> out_distances_sq_c(num_closest_c);

		num_closest_c = caustic_tree->knnSearch(&query_point[0],num_closest_c,&out_indices_c[0],&out_distances_sq_c[0]);
		out_indices_c.resize(num_closest_c);
		out_distances_sq_c.resize(num_closest_c);
		float radius_sq_c = *max_element(out_distances_sq_c.begin(),out_distances_sq_c.end());
		float circle_area_c  = M_PI*radius_sq_c;
		float k = 1;
		for(size_t i = 0; i<num_closest_c;i++)
		{
			
			size_t current_index = out_indices_c[i];
			float x_c = caustic_map.pts[current_index].x;
			float y_c = caustic_map.pts[current_index].y;
			float z_c = caustic_map.pts[current_index].z;
			Vector dir_c = caustic_map.pts[current_index].incident_direction;
			Colour pow_c = caustic_map.pts[current_index].power;
			int caustic_map_size = caustic_map.pts.size();
			float dx = abs(x_c - hit.position.x);
			float dy = abs(y_c - hit.position.y);
			float dz = abs(z_c - hit.position.z);

			float dis = sqrt(dx*dx+dy*dy+dz*dz);
			float wpc = 1 - dis/(k*sqrt(radius_sq_c));
			L = dir_c;
			L.negate();
			L.normalise();

			V.x = eye.x - x_c;
			V.y = eye.y - y_c;
			V.z = eye.z - z_c;
			V.normalise();
			R = dir_c - 2.0f*(dir_c.dot(N))*N;
			R.normalise();

			RdotV = max(R.dot(V),0.0f);
			LdotN = max(L.dot(N),0.0f);
			// Phong reflection without ambient term.
			red_c 		= min(0.5f*(wpc*diffuse.r*(LdotN)*pow_c.r + specular.r*pow(RdotV,20.0f)*pow_c.r),1.0f); // caustics do not depend on the object they hit.
			green_c 	= min(0.5f*(wpc*diffuse.g*(LdotN)*pow_c.g + specular.g*pow(RdotV,20.0f)*pow_c.g),1.0f);
			blue_c 		= min(0.5f*(wpc*diffuse.b*(LdotN)*pow_c.b + specular.b*pow(RdotV,20.0f)*pow_c.b),1.0f);


		
			// maybe put a max on caustic contribution
			caustic_col.r += red_c;
			caustic_col.g += green_c;
			caustic_col.b += blue_c;
		}
		float denom = (1-2/(3*k))*circle_area_c;
		caustic_col.r = min(caustic_col.r/denom,1.0f);
		caustic_col.g = min(caustic_col.g/denom,1.0f);
		caustic_col.b = min(caustic_col.b/denom,1.0f);

	}

	overall_colour.add(direct);
	overall_colour.add(indirect);
	overall_colour.add(global);
	overall_colour.add(caustic_col);
	overall_colour.r = min(overall_colour.r,20.0f);
	overall_colour.g = min(overall_colour.g,20.0f); // 20 works well.
	overall_colour.b = min(overall_colour.b,20.0f); // These are unecessary 

	
	
	caustic_output.r = min(caustic_col.r,20.0f);
	caustic_output.g = min(caustic_col.g,20.0f);
	caustic_output.b = min(caustic_col.b,20.0f);
	
	global_map_output.r = min(indirect.r,20.0f);
	global_map_output.g = min(indirect.g,20.0f);
	global_map_output.b = min(indirect.b,20.0f);


	direct_output = direct;
	
}

