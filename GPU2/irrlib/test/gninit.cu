#include <cstdio>
#include <cassert>

#include <cutil.h>

#include "cuda_pointer.h"
#include "gnutil.h"
#include "gninit.h"

#define PROFILE
#ifdef PROFILE
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}
struct Timer{
	const char *name;
	FILE       *fp;
	const char *format;
	double tstart;

	Timer(
			const char *_name,
			FILE       *_fp     = stdout,
			const char *_format = " %-10s : %f sec\n")
	   : name(_name), fp(_fp), format(_format)
	{
		tstart = get_wtime();
	}
	~Timer(){
		double tend = get_wtime();
		fprintf(fp, format, name, tend - tstart);
		fflush(fp);
	}
};
#else
static double get_wtime(){
	return 0.0;
}
struct Timer{
	Timer(
			const char *_name,
			FILE       *_fp     = stdout,
			const char *_fotmat = " %s : %f sec\n")
	{}
	~Timer(){}
};
#endif

typedef Gvec3<float   > fvec3;
typedef Gvec3<DblFloat> dvec3;
typedef DblFloat        twofl;

struct PosH2{
	dvec3 pos;
	float h2;
	float pad;

	__host__ PosH2(const double _pos[3], const double _h2) :
		pos(_pos), h2(_h2) {}
};

__global__ void kernel_count_neib(
		const int    nbody,
		const float  eps2,
		const PosH2  ptcl[],
		_out_ int    nnb_out[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const PosH2 ip = ptcl[tid];
	int nnb = 0;

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const PosH2 jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const float r2 = eps2 + dr*dr;
		if((j != tid) && (r2 < ip.h2)) nnb++;
	}

	nnb_out[tid] = nnb;
}

__global__ void kernel_get_neib(
		const int    nbody,
		const float   eps2,
		const PosH2  ptcl   [],
		const int    nboff  [],
		_out_ int    nnb_out[],
		_out_ int    nblist [])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const PosH2 ip = ptcl[tid];
	const int ioff = nboff[tid];
	int *nbdst = nblist + ioff;
	int nnb = 0;

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const PosH2 jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const float r2 = eps2 + dr*dr;
		if((j != tid) && (r2 < ip.h2)){
			*nbdst++ = j;
			nnb++;
		}
	}

	nnb_out[tid] = nnb;
}

struct PosM{
	dvec3 pos;
	float mass;
	float pad;

	PosM(const double x[3], const double m) : pos(x), mass(m) {}
	PosM(int) : pos(dvec3(0.0)), mass(0.0f), pad(0.0f) {}
	__device__ PosM() {}
};

__global__ void kernel_calc_pot(
		const int    nbody,
		const float  eps2,
		const PosM   ptcl[],
		_out_ twofl  pot_out[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const PosM ip = ptcl[tid];
	twofl pot(0.f, 0.f);

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const PosM  jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const float r2 = eps2 + dr*dr;
		const float pij = jp.mass * rsqrtf(r2);
		if(j != tid) pot += pij;
	}

	pot_out[tid] = pot;
}

struct Particle{
	dvec3  pos;  // 6
	fvec3  vel;  // 9
	float  mass; // 10
	float2 pad;  // 12

	Particle(
			const double _pos[3], 
			const double _vel[3], 
			const double _mass) 
		: pos(_pos), vel(_vel), mass(_mass) {}
	__device__ Particle() {}
};

#if 0
struct Force{
	dvec3    acc; // 6
	fvec3    jrk; // 9
	float3   pad; // 12

	__device__ 
	Force() : acc(twofl(0.0f, 0.0f)), jrk(0.0f) {}

	__device__ 
	void accumulate(const fvec3 &a, const fvec3 &j)
	{
		acc += a;
		jrk += j;
	}
};
#else
struct Force{
	dvec3    acc; // 6
	dvec3    jrk; // 12

	__device__ 
	Force() : acc(twofl(0.0f, 0.0f)), jrk(twofl(0.0f, 0.0f)) {}

	__device__ 
	void accumulate(const fvec3 &a, const fvec3 &j)
	{
		acc += a;
		jrk += j;
	}
};
#endif

__global__ void kernel_calc_force(
		const int      nbody,
		const float    eps2,
		const Particle ptcl [],
		_out_ Force    force[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const Particle ip = ptcl[tid];
	Force fo; // initialized by the constructor

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const Particle jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const fvec3 dv = jp.vel - ip.vel;

		const float r2 = eps2 + dr * dr;
		const float rv = dr * dv;
		const float rinv1 = (tid != j) ? rsqrtf(r2) : 0.0f;

		const float rinv2 = rinv1 * rinv1;
		const float mrinv1 = jp.mass * rinv1;
		const float mrinv3 = mrinv1 * rinv2;
		const float alpha  = -3.f * rv * rinv2;

		const fvec3 acc = mrinv3 * dr;
		const fvec3 jrk = mrinv3 * (dv + alpha * dr);

		fo.accumulate(acc, jrk);
	}

	force[tid] = fo;
}

struct FatParticle{
	dvec3  pos;  // 6
	fvec3  vel;  // 9
	fvec3  acc;  // 12
	fvec3  jrk;  // 15
	float  h2;   // 16
	float  mass; // 17
	float3 pad;  // 20

	FatParticle(
			const double _pos[3], 
			const double _vel[3], 
			const double _acc[3], 
			const double _jrk[3], 
			const double _h2, 
			const double _mass) 
		: pos(_pos), vel(_vel), acc(_acc), jrk(_jrk), h2(_h2), mass(_mass) {}
	__device__ FatParticle() {}
};

struct FatForce{
	dvec3    acc; // 6
	fvec3    jrk; // 9
	fvec3    snp; // 12
	fvec3    crk; // 15
	float    pad; // 16

	__device__ 
	FatForce() : acc(twofl(0.0f, 0.0f)), jrk(0.0f), snp(0.0f), crk(0.0f) {}

	__device__ 
	void accumulate(
			const fvec3 &a, 
			const fvec3 &j, 
			const fvec3 &s, 
			const fvec3 &c)
	{
		acc += a;
		jrk += j;
		snp += s;
		crk += c;
	}
};

__global__ void kernel_calc_fpoly(
		const int         nbody,
		const float       eps2,
		const FatParticle ptcl[],
		_out_ FatForce    force_reg[],
		_out_ FatForce    force_irr[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const FatParticle ip = ptcl[tid];
	FatForce freg; // initialized by the constructor
	FatForce firr; // initialized by the constructor

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const FatParticle jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const fvec3 dv = jp.vel - ip.vel;
		const fvec3 da = jp.acc - ip.acc;
		const fvec3 dj = jp.jrk - ip.jrk;

		const float r2 = eps2 + dr * dr;
		const float rv = dr * dv;
		const float v2 = dv * dv;
		const float ra = dr * da;
		const float va = dv * da;
		const float rj = dr * dj;

		const float rinv1 = (tid != j) ? rsqrtf(r2) : 0.0f;

		const float rinv2 = rinv1 * rinv1;
		const float mrinv1 = jp.mass * rinv1;
		const float mrinv3 = mrinv1 * rinv2;

		const float alpha  = rv * rinv2;
		const float alphalpha = alpha * alpha;
		const float beta   = ((v2 + ra) * rinv2 + alphalpha);
		const float gamma  = ((3.f*va + rj) * rinv2 + alpha * (3.f*beta - 4.f*alphalpha));

		const fvec3 acc = mrinv3*dr;
		const fvec3 jrk = mrinv3*dv - (3.f*alpha)*acc;
		const fvec3 snp = mrinv3*da - (6.f*alpha)*jrk - (3.f*beta)*acc;
		const fvec3 crk = mrinv3*dj - (9.f*alpha)*snp - (9.f*beta)*jrk - (3.f*gamma)*acc;

		if(r2 < ip.h2){
			firr.accumulate(acc, jrk, snp, crk);
		}else{
			freg.accumulate(acc, jrk, snp, crk);
		}
	}

	force_reg[tid] = freg;
	force_irr[tid] = firr;
}

struct PosVelH2{
	dvec3 pos;
	fvec3 vel;
	float h2;
	float pad[2];

	PosVelH2(
			const double _pos[3], 
			const double _vel[3], 
			const double _h2) 
		: pos(_pos), vel(_vel), h2(_h2) {}
};

__global__ void kernel_calc_dtmin(
		const int         nbody,
		const float       eps2,
		const PosVelH2    ptcl[],
		_out_ float2      dtmin[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const PosVelH2 ip = ptcl[tid];
	float dtreg = 16777216.;
	float dtirr = 16777216.;

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const PosVelH2 jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const fvec3 dv = jp.vel - ip.vel;
		const float r2 = eps2 + dr * dr;
		const float v2 = dv * dv;
		const float dtij = r2 / v2;

		if(j != tid){
			if(r2 < ip.h2){
				dtirr = min(dtirr, dtij);
			}else{
				dtreg = min(dtreg, dtij);
			}
		}
	}
	dtmin[tid] = make_float2(sqrtf(dtreg), sqrtf(dtirr));
}

__global__ void kernel_count_friend(
		const int      nbody,
		const float    eps2,
		const PosVelH2 ptcl[],
		const float    dt_ov_eta,
		_out_ int      nfr_out[])
{
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid >= nbody) return;
	const PosVelH2 ip = ptcl[tid];
	const float crit = dt_ov_eta * dt_ov_eta;
	int nfr = 0;

#pragma unroll 4
	for(int j=0; j<nbody; j++){
		const PosVelH2 jp = ptcl[j];
		const fvec3 dr = jp.pos - ip.pos;
		const fvec3 dv = jp.vel - ip.vel;
		const float r2 = eps2 + dr * dr;
		const float v2 = dv * dv;
		const float dtij = r2 / v2;
		if((j != tid) && (dtij < crit)) nfr++;
	}

	nfr_out[tid] = nfr;
}

namespace gninit{
	void count_neib(
			const int    nbody,
			const double eps2,
			const double pos[][3],
			const double h2 [],
			_out_ int    nnb_out[])
	{
		Timer timer(__func__);

		cudaPointer<PosH2, true> posh2(nbody);
		cudaPointer<int  , true> nnb  (nbody);

		for(int i=0; i<nbody; i++){
			posh2[i] = PosH2(pos[i], h2[i]);
		}
		posh2.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_count_neib <<<nblock, nthread>>>
			(nbody, eps2, posh2, nnb);

		nnb.dtoh();
		for(int i=0; i<nbody; i++){
			nnb_out[i] = nnb[i];
		}
	}

	void get_neib(
			const int    nbody,
			const double eps2,
			const double pos[][3],
			const double h2 [],
			const int    nnb_in[],
			_out_ int * const list[]) // array of pointers
	{
		Timer timer(__func__);

		cudaPointer<PosH2, true> posh2 (nbody);
		cudaPointer<int  , true> nboff (nbody);
		cudaPointer<int  , true> nnb   (nbody);
		cudaPointer<int  , true> nblist;

		int nbsum = 0;
		for(int i=0; i<nbody; i++){
			posh2[i] = PosH2(pos[i], h2[i]);
			nboff[i] = nbsum;
			nbsum += nnb_in[i];
		}
		nblist.allocate(nbsum);
		fprintf(stderr, "nbsum = %d\n", nbsum);

		posh2.htod();
		nboff.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_get_neib <<<nblock, nthread>>>
			(nbody, eps2, posh2, nboff, nnb, nblist);

		nnb.dtoh();
		nblist.dtoh();

		for(int i=0; i<nbody; i++){
			assert(nnb_in[i] == nnb[i]);
			const int len = nnb[i];
			int *src = &nblist[nboff[i]];
			int *dst = list[i];
			for(int k=0; k<len; k++){
				dst[k] = src[k];
			}
		}
	}

	void calc_pot(
			const int    nbody,
			const double eps2,
			const double mass[],
			const double pos [][3],
			_out_ double pot_out[])
	{
		Timer timer(__func__);

		cudaPointer<PosM,  true> posm(nbody);
		cudaPointer<twofl, true> pot (nbody);

		for(int i=0; i<nbody; i++){
			posm[i] = PosM(pos[i], mass[i]);
		}
		posm.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_calc_pot <<<nblock, nthread>>>
			(nbody, eps2, posm, pot);

		pot.dtoh();
		for(int i=0; i<nbody; i++){
			pot_out[i] = -double(pot[i]);
		}
	}

	void calc_force(
			const int    nbody,
			const double eps2,
			const double mass[],
			const double pos [][3],
			const double vel [][3],
			_out_ double acc [][3],
			_out_ double jrk [][3])
	{
		Timer timer(__func__);

		cudaPointer<Particle, true> ptcl (nbody);
		cudaPointer<Force,    true> force(nbody);

		for(int i=0; i<nbody; i++){
			ptcl[i] = Particle(pos[i], vel[i], mass[i]);
		}
		ptcl.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_calc_force <<<nblock, nthread>>>
			(nbody, eps2, ptcl, force);

		force.dtoh();
		for(int i=0; i<nbody; i++){
			force[i].acc.write(acc[i]);
			force[i].jrk.write(jrk[i]);
		}
	}

	void calc_fpoly(
			const int    nbody,
			const double eps2,
			const double mass[],
			const double pos [][3],
			const double vel [][3],
			const double acc [][3],
			const double jrk [][3],
			const double h2 [],
			_out_ double fpoly_reg[][4][3],
			_out_ double fpoly_irr[][4][3])
	{
		Timer timer(__func__);

		cudaPointer<FatParticle, true> ptcl(nbody);
		cudaPointer<FatForce,    true> freg(nbody);
		cudaPointer<FatForce,    true> firr(nbody);

		for(int i=0; i<nbody; i++){
			ptcl[i] = FatParticle(pos[i], vel[i], 
					      acc[i], jrk[i], h2[i], mass[i]);
		}
		ptcl.htod();

		const int nthread = 64;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_calc_fpoly <<<nblock, nthread>>>
			(nbody, eps2, ptcl, freg, firr);

		freg.dtoh();
		firr.dtoh();
		for(int i=0; i<nbody; i++){
			freg[i].acc.write(fpoly_reg[i][0]);
			freg[i].jrk.write(fpoly_reg[i][1]);
			freg[i].snp.write(fpoly_reg[i][2]);
			freg[i].crk.write(fpoly_reg[i][3]);
			firr[i].acc.write(fpoly_irr[i][0]);
			firr[i].jrk.write(fpoly_irr[i][1]);
			firr[i].snp.write(fpoly_irr[i][2]);
			firr[i].crk.write(fpoly_irr[i][3]);
		}
	}
	void calc_dtmin(
			const int    nbody,
			const double eps2,
			const double pos [][3],
			const double vel [][3],
			const double h2 [],
			_out_ double dtreg[],
			_out_ double dtirr[])
	{
		Timer timer(__func__);

		cudaPointer<PosVelH2, true> ptcl (nbody);
		cudaPointer<float2,   true> dtmin(nbody);

		for(int i=0; i<nbody; i++){
			ptcl[i] = PosVelH2(pos[i], vel[i], h2[i]);
		}
		ptcl.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_calc_dtmin <<<nblock, nthread>>>
			(nbody, eps2, ptcl, dtmin);

		dtmin.dtoh();
		for(int i=0; i<nbody; i++){
			dtreg[i] = dtmin[i].x;
			dtirr[i] = dtmin[i].y;
		}
	}
	void count_friend(
			const int    nbody,
			const double eps2,
			const double pos [][3],
			const double vel [][3],
			const double dt_ov_eta,
			_out_ int    nfr_out[])
	{
		Timer timer(__func__);

		cudaPointer<PosVelH2, true> ptcl(nbody);
		cudaPointer<int   ,   true> nfr (nbody);

		for(int i=0; i<nbody; i++){
			ptcl[i] = PosVelH2(pos[i], vel[i], 0.0);
		}
		ptcl.htod();

		const int nthread = 128;
		const int nblock  = 1 + (nbody-1)/nthread;
		kernel_count_friend <<<nblock, nthread>>>
			(nbody, eps2, ptcl, dt_ov_eta, nfr);

		nfr.dtoh();
		for(int i=0; i<nbody; i++){
			nfr_out[i] = nfr[i];
		}
	}
}
