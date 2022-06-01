#ifndef RAYTRACING_PHOTONMAP_H
#define RAYTRACING_PHOTONMAP_H

#include <ANN/ANN.h>
#include "Vector.hpp"

const float M_PI = 3.141592653589793f;

typedef struct Photon {
    float pos[3];
    float power[3];
    //short plane
    unsigned char theta, phi;
} Photon;

typedef struct NearestPhotons {
	//int max; //I don't think we need this; it seems ANN does not use heap
	int found;
	//int got_heap;
	float pos[3];	
	float *dist2; // I think dist2[0] is to store the max_dist2 while following 
                  // elements are really storing distance of each nearest photon?
	const Photon **index;
} NearestPhotons;


class PhotonMap {
    public:
        PhotonMap(const int max_num_photons);
        ~PhotonMap() {
            free(photons);
            //free(kd_tree);
        }
        void store(const float power[3], const float pos[3], const float dir[3]);
        void photon_dir ( float *dir, const Photon *p) const;

        //my new function
        void buildKdtree(); //hopefully this can replace balance()
        void locate_photons(NearestPhotons* np) const;
        void irradiance_estimate(float irrad[3],
	                             const float pos[3],
	                             const float normal[3],
	                             const float max_dist
	                             /*,const int nphotons*/) const;

        int photon_stored;
        int photon_max;
        Photon* photons;

    private:
        //used for scaling
        int prev_scaled;

        //store sine and cosine values for convenience
        float costheta[256];
	    float sintheta[256];
	    float cosphi[256];
	    float sinphi[256];

        //I think this means the min and max coords of the scene
        float bbox_min[3];
        float bbox_max[3];

        //our kdTree
        ANNkd_tree* kd_tree;
};

#endif
