#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include "vector3.h"

#define _out_
#define PROFILE

#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
#else
static double get_wtime(){
	return 0.0;
}
#endif

struct Particle{
	dvec3  pos;
	dvec3  vel;
	dvec3  acc2;
	dvec3  jrk6;
	double mass;
	double time;

	Particle(){}
	Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time) : 
		pos(_pos), vel(_vel), acc2(_acc2), jrk6(_jrk6), mass(_mass), time(_time)
	{
	}
};

struct Predictor{
	dvec3  pos;
	dvec3  vel;
	double mass;

	Predictor() {}
	Predictor(
			const Particle &p, 
			const double ti)
	{
		const double dt = ti - p.time;
		pos  = p.pos + dt*(p.vel + dt*(p.acc2 + dt*(p.jrk6)));
		vel  = p.vel + (2.0*dt)*(p.acc2 + (1.5*dt)*(p.jrk6));
		mass = p.mass;
	}
};

struct NBlist{
	enum{ NB_MAX = 511 };
	int nnb;
	int nb[NB_MAX];

	void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

static std::vector<Particle>  ptcl;
static std::vector<Predictor> pred;
static std::vector<NBlist   > list;
static std::vector<int>       counter;

static double time_pred, time_pact, time_grav, time_onep;
static unsigned long long num_inter;

static void gpuirr_open(
		const int nmax,
		const int lmax)
{
	assert(lmax <= 1 + NBlist::NB_MAX);

	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Opening GPUIRR lib. CPU ver. \n"); 
	fprintf(stderr, " nmax = %d, lmax = %d\n", nmax, lmax);
	fprintf(stderr, "**************************** \n"); 

	ptcl.resize(nmax);
	pred.resize(nmax);
	list.resize(nmax);
	counter.resize(NBlist::NB_MAX, 0);

	time_pred = time_pact =  time_grav = time_onep = 0.0;
	num_inter = 0;
}

static void gpuirr_close(){
	FILE *fp = fopen("NBSTAT", "w");
	fprintf(fp, "   #NB   #CALL\n"); 
	for(int i=0; i<NBlist::NB_MAX; i++){
		fprintf(fp, "%6d%8d\n", i, counter[i]);
	}
	fclose(fp);

	const double Gflops = 60.0 * double(num_inter) * 1.e-9 / time_grav;

	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Closing GPUIRR lib. CPU ver. \n"); 
	fprintf(stderr, "time pred  : %f sec\n", time_pred);
	fprintf(stderr, "time pact  : %f sec\n", time_pact);
	fprintf(stderr, "time grav  : %f sec\n", time_grav);
	fprintf(stderr, "time onep  : %f sec\n", time_onep);
	fprintf(stderr, "perf grav  : %f Gflops\n", Gflops);
	fprintf(stderr, "**************************** \n"); 
}

static void gpuirr_set_jp(
		const int addr,
		const double pos [3],
		const double vel [3],
		const double acc2[3],
		const double jrk6[3],
		const double mass,
		const double time)
{
	ptcl[addr] = Particle(pos, vel, acc2, jrk6, mass, time);
}

static void gpuirr_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
	assert(nnb <= NBlist::NB_MAX);
	list[addr].nnb = nnb;
	// int *dst = list[addr].nb;
	for(int k=0; k<nnb; k++){
		list[addr].nb[k] = nblist[k] - 1;
		// dst[k] = nblist[k] - 1;
	}
}

static void gpuirr_pred_all(
		const int    js,
		const int    je,
		const double ti)
{
	const int nmax = ptcl.size();
	assert(js >= 0);
	assert(je <= nmax);
	const double t0 = get_wtime();
#pragma omp parallel for
	for(int j=js; j<je; j++){
		pred[j] = Predictor(ptcl[j], ti);
	}
	const double t1 = get_wtime();
	::time_pred += t1-t0;
}

static void gpuirr_pred_act(
		const int ni,
		const int addr[],
		const double ti)
{
	static std::vector<int> slist(1024), ulist(1024);
	const double t0 = get_wtime();
	slist.clear();
	ulist.clear();
	for(int i=0; i<ni; i++){
		const int ii = addr[i]-1;
		slist.push_back(ii);
		for(int k=0; k<list[ii].nnb; k++){
			slist.push_back(list[ii].nb[k]);
		}
	}
	std::sort(slist.begin(), slist.end());
#if 0
	ulist = slist;
#else
	ulist.resize( slist.size() );
	std::vector<int>::iterator it
		= std::unique_copy(slist.begin(), slist.end(), ulist.begin());
	// cut the tail (resize)
	ulist.erase(it, ulist.end());
#endif

#if 0
	std::cerr << "slist : " << slist.size() << "    "
	          << "ulist : " << ulist.size() << std::endl;
#endif

	const int npred = ulist.size();
#pragma omp parallel for
	for(int i=0; i<npred; i++){
		const int j = ulist[i];
		const int j1 = ulist[i+1];
		__builtin_prefetch(&ptcl[j1]);
		__builtin_prefetch(&pred[j1]);
		pred[j] = Predictor(ptcl[j], ti);
	}
	const double t1 = get_wtime();
	::time_pact += t1-t0;
}

static void gpuirr_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3])
{
	// list[addr].print(addr);

	const int jmax = ptcl.size();
	const Predictor &ip = pred[addr];
	dvec3 acc(0.0), jrk(0.0);
	const int nnb = list[addr].nnb;
	counter[nnb]++;
	for(int k=0; k<nnb; k++){
		const int j = list[addr].nb[k];
		assert(j < jmax);
		Predictor &jp = pred[j];
		const dvec3  dr = jp.pos - ip.pos;
		const dvec3  dv = jp.vel - ip.vel;
		const double r2 = dr * dr;
		const double rv = dr * dv;
		const double rinv2 = 1.0 / r2;
		const double rinv  = sqrt(rinv2);
		const double alpha = -3.0 * rinv2 * rv;
		const double mrinv3 = jp.mass * rinv * rinv2;

		acc += mrinv3 * dr;
		jrk += mrinv3 * (dv + alpha * dr);
	}
	acc.store(accout);
	jrk.store(jrkout);
}

static void gpuirr_firr_vec(
		const int ni,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3])
{
	const double t0 = get_wtime();
	int nint = 0;
#pragma omp parallel for reduction(+: nint)
	for(int i=0; i<ni; i++){
		gpuirr_firr(addr[i]-1, accout[i], jrkout[i]);
		nint += list[addr[i]-1].nnb;
	}
	num_inter += nint;
	const double t1 = get_wtime();
	::time_grav += t1-t0;
}

extern "C"{
	void gpuirr_open_(int *nmax, int *lmax){
		gpuirr_open(*nmax, *lmax);
	}
	void gpuirr_close_(){
		gpuirr_close();
	}
	void gpuirr_set_jp_(
		int    *addr,
		double  pos [3],
		double  vel [3],
		double  acc2[3],
		double  jrk6[3],
		double *mass,
		double *time)
	{
		gpuirr_set_jp((*addr)-1, pos, vel, acc2, jrk6, *mass, *time);
	}
	void gpuirr_set_list_(
		int *addr,
		int *nblist)
	{
		gpuirr_set_list((*addr)-1, *nblist, nblist+1);
	}
	void gpuirr_pred_all_(
			int    *js,
			int    *je,
			double *ti)
	{
		gpuirr_pred_all((*js)-1, *(je),  *ti);
	}
	void gpuirr_pred_act_(
			int    *ni,
			int     addr[],
			double *ti)
	{
		gpuirr_pred_act(*ni, addr, *ti);
	}
	void gpuirr_firr_(
			int    *addr,
			double  acc[3],
			double  jrk[3])
	{
		const double t0 = get_wtime();
		gpuirr_firr((*addr)-1, acc, jrk);
		const double t1 = get_wtime();
		::time_onep += t1-t0;
	}
	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3])
	{
		gpuirr_firr_vec(*ni, addr, acc, jrk);
	}
}
