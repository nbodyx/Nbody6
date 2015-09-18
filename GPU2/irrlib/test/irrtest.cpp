#include <iostream>
#include <vector>
#include "vector3.h"
#include "plummer.h"
#include "gninit.h"
#include "gpuirr.h"

struct NBlist{
	enum{NBMAX = 128};
	int list[NBMAX];
	void print(const int i, const int nnb, FILE *fp = stdout){
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", list[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

int main(int argc, char **argv){
	if(argc < 2) return 1;

	const int nk = atoi(argv[1]);
	assert(1<=nk && nk<=1024);
	const int nbody = nk * 1024;

	const double rs = 0.15;

	Plummer pl(nbody); // has pos, mass

	std::vector<double> h2 (nbody, rs*rs);
	std::vector<int   > nnb(nbody);
	for(;;){
		gninit::count_neib(nbody, 0.0, pl.pos[0], &h2[0], &nnb[0]);
		bool overflow = false;
		for(int i=0; i<nbody; i++){
			if(nnb[i] > NBlist::NBMAX){
				overflow = true;
				const double ratio = double(NBlist::NBMAX) / double(nnb[i]);
				h2[i] *= pow(ratio, 1.5) * 0.95;
			}
		}
		if(!overflow) break;
	}

	std::vector<NBlist> nblist(nbody);
	std::vector<int * > nbptr (nbody);
	for(int i=0; i<nbody; i++){
		nbptr[i] = nblist[i].list;
	}
	gninit::get_neib(nbody, 0.0, pl.pos[0], &h2[0], &nnb[0], &nbptr[0]);

	std::vector<dvec3> acc(nbody);
	std::vector<dvec3> jrk(nbody);
	gninit::calc_force(nbody, 0.0, &pl.mass[0], pl.pos[0], pl.vel[0], 
	                   acc[0], jrk[0]);

	{
		int nmax = nbody;
		int lmax = NBlist::NBMAX + 1;
		gpuirr_open_(&nmax, &lmax);

		for(int i=0; i<nbody; i++){
			int addr = i+1;
			dvec3 pos = pl.pos[i];
			dvec3 vel = pl.vel[i];
			dvec3 acc2 = (1./2.) * acc[i];
			dvec3 jrk6 = (1./6.) * jrk[i];
			double mass = pl.mass[i];
			double time = 0.0;
			gpuirr_set_jp_(&addr, pos, vel, acc2, jrk6, &mass, &time);

			static int list[NBlist::NBMAX + 1];
			list[0] = nnb[i];
			const int len = nnb[i];
			for(int k=0; k<len; k++){
				list[k+1] = nblist[i].list[k] + 1;
			}
			gpuirr_set_list_(&addr, list);
		}
		int js = 1, je = nbody;
		double time = 1./64.;
		// double time = 0.;
		gpuirr_pred_all_(&js, &je, &time);

		int ni = 10;
		int addr[10] = {1, 2, 3, 4, 6, 5, 7, 8, 9, 10};
		dvec3 irracc[10], irrjrk[10];
		gpuirr_firr_vec_(&ni, addr, irracc[0], irrjrk[0]);

		for(int i=0; i<10; i++){
			std::cout << addr[i] << " : "
					  << nnb[addr[i]-1] << " : "
					  // << pl.pos[i] << ", "
			          << irracc[i] << ", "
			          << irrjrk[i] << std::endl;
		}

		gpuirr_close_();
	}

	return 0;
}
