#include "photonmap.hpp"

//constructor
PhotonMap::PhotonMap (const int max_num_photons){
	photon_stored = 0;
	prev_scaled = 1;
	photon_max = max_num_photons;
	
	photons = (Photon*) malloc(sizeof(Photon)*(photon_max+1));

	if (photons == nullptr) {
		fprintf(stderr, "Out of memory initializing photon map\n");
		exit(-1);
	}
 
	bbox_min[0] = bbox_min[1] = bbox_min[2] = 1e8f;
	bbox_max[0] = bbox_max[1] = bbox_max[2] = -1e8f;

    //values for all angles
	for (int i = 0; i<256; i++) {
		double angle = double(i)*(1.0/256.0)*M_PI;
		costheta[i] = cos(angle);
		sintheta[i] = sin(angle);
		cosphi[i] = cos(2.0*angle);
		sinphi[i] = sin(2.0*angle);
	}
}

//store into Photon array
void PhotonMap::store(const float power[3], const float pos[3], const float dir[3]) {

    	if (photon_stored>photon_max)
		return;
	photon_stored++;
	Photon *const node = &photons[photon_stored];

	for (int i = 0; i<3 ;i++){
		node->pos[i] =  pos[i];
		
		if (node->pos[i]<bbox_min[i])
			bbox_min[i] = node->pos[i];
		if (node->pos[i]>bbox_max[i])
			bbox_max[i] = node->pos[i];

		node->power[i] = power[i];
	}

	int theta = int(acos(dir[2])*(256.0/M_PI) );
	if (theta>255)
		node->theta = 255;
	else
		node->theta = (unsigned char) theta;

	int phi = int(atan2(dir[1],dir[0])*(256.0/(2.0*M_PI)) );
	if (phi>255)
		node->phi = 255;
	else if (phi<0)
		node->phi = (unsigned char)(phi+256);
	else
		node->phi = (unsigned char)phi;
}

void PhotonMap::buildKdtree() {
    ANNpointArray dataPts = annAllocPts(photon_max, 3);

    //maybe unreadable points? dont know yet
    for (int i = 0; i < photon_max; i++) {
        for (int j = 0; j < 3; j++) {
            dataPts[i][j] = photons[i].pos[j];
        }
    }

    kd_tree = new ANNkd_tree(dataPts, photon_max, 3);
}

//I think this asks us to locate photons given a certain radius?
void PhotonMap::locate_photons(NearestPhotons* np) const {\
    ANNpoint q;
    for (int i = 0; i < 3; i++) {
        q[i] = np->pos[i];
    }
    double eps = 0.0;
    double sqRad = np->dist2[0];
    int k = 0; //to initiate the real radius-based search instead of k-nearest; see ANN doc
    ANNidxArray nn_idx = NULL;
    ANNdistArray dd = NULL;
    int ptsInRad = kd_tree->annkFRSearch(q, sqRad, k, nn_idx, dd, eps); //return result in nn_idx and dist_idx
    
    np->found = ptsInRad;
    for (int i = 0; i < ptsInRad; i++) {
        np->index[i] = &photons[nn_idx[i]]; //?
        np->dist2[i+1] = dd[i];
    }
}

void PhotonMap::photon_dir ( float *dir, const Photon *p) const {
	dir[0] = sintheta[p->theta]*cosphi[p->phi];
	dir[1] = sintheta[p->theta]*sinphi[p->phi];
	dir[2] = costheta[p->theta];
}

void PhotonMap::irradiance_estimate(
	float irrad[3],
	const float pos[3],
	const float normal[3],
	const float max_dist
	/*,const int nphotons*/) const 
{
	irrad[0] = irrad[1] = irrad[2] = 0.0;

	NearestPhotons np;
	np.dist2 = NULL;
	np.index = NULL;

	np.pos[0] = pos[0]; np.pos[1] = pos[1]; np.pos[2] = pos[2];
	//np.max = nphotons;
	np.found = 0;
	//np.got_heap = 0;
	np.dist2[0] = max_dist*max_dist;

	locate_photons(&np);
	
	if (np.found<8)
		return;
	float pdir[3];

	for (int i=1;i<=np.found; i++){
		const Photon *p = np.index[i];
		photon_dir(pdir,p);
		if ((pdir[0]*normal[0]+pdir[1]*normal[1]+pdir[2]*normal[2])<0.0f) {
			irrad[0] += p->power[0];
			irrad[1] += p->power[1];
			irrad[2] += p->power[2];
		}
	}
	
	const float tmp = (1.0f/M_PI)/(np.dist2[0]);
	
	irrad[0] *= tmp;
	irrad[1] *= tmp;
	irrad[2] *= tmp;
    
}