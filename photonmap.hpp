// referred to Jensen's code

#ifndef RAYTRACING_PHOTONMAP_H
#define RAYTRACING_PHOTONMAP_H

#include "Vector.hpp"
#include "Bounds3.hpp"

typedef struct Photon {
    float pos[3];
    float power[3];
    short plane; // splitting plane for kd-tree
    unsigned char theta, phi; // incoming direction
} Photon;

typedef struct NearestPhotons {
	int max;
	int found;
	int got_heap;
	float pos[3];	
	float *dist2; // max search distance d^2
	const Photon **index; // to keep track of the current subtree's root
} NearestPhotons;

class PhotonMap {
public:
    PhotonMap(const int max_phot);
    ~PhotonMap();
    void store(const float power[3], const float pos[3], const float dir[3]);
    
    // 1 / (# of emitted photons)
    void scale_photon_power(const float scale);
    
    // balance the kd-tree
    void balance(void);

    // compute the irradiance at a given position
    void irradiance_estimate(float irrad[3], // return irradiance
                                const float pos[3], // surface position
                                const float normal[3], // surface normal at pos
                                const float max_dist, // max distance to look for photons
                                const int nphotons ) const; // # of photons to use (maybe not necessary?)


    // locate_photons(1) initiate search at the root of the tree
    void locate_photons(NearestPhotons* const np, const int index) const;
    
    // return direction of photon
    void photon_dir(float *dir, const Photon *p) const;

    // int photon_stored;
    // int photon_max;
    // Photon* photons;

private:
    void balance_segment(Photon **pbal, Photon **porg, const int index, const int start, const int end);
    void median_split(Photon **p, const int start, const int end, const int median, const int axis);

    Photon *photons;

    int stored_photons;
    int half_stored_photons;
    int max_photons;
    int prev_scale; // used for scaling

    //store sine and cosine values for convenience
    float costheta[256];
    float sintheta[256];
    float cosphi[256];
    float sinphi[256];

    // not sure if these are necessary at the moment
    float bbox_min[3];
    float bbox_max[3];

};

#endif
