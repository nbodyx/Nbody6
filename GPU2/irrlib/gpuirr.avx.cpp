#define __USE_GNU
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <omp.h>
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

typedef float  v4sf __attribute__((vector_size(16)));
typedef double v2df __attribute__((vector_size(16)));
typedef float  v8sf __attribute__((vector_size(32)));
typedef double v4df __attribute__((vector_size(32)));

#define REP4(x) {x,x,x,x}
#define REP8(x) {x,x,x,x,x,x,x,x}

static inline v8sf rsqrt_NR(const v8sf x){
#if 1
	const v8sf y = __builtin_ia32_rsqrtps256(x);
	const v8sf c1 = REP8(-0.5f);
	const v8sf c2 = REP8(-3.0f);
	return (c1 * y) * (x*y*y + c2);
#else
	return __builtin_ia32_rsqrtps_nr256(x);
#endif
}

static inline v8sf pack_2_v4sf(const v4sf a, const v4sf b){
	// v8sf p;
	v8sf p = REP8(0.0f); // just avoid warning
	p = __builtin_ia32_vinsertf128_ps256(p, a, 0);
	p = __builtin_ia32_vinsertf128_ps256(p, b, 1);
	return p;
}

// more efficieent?
static inline v8sf pack_2_v4sf(const v4sf a, const v4sf b, v8sf &p){
	p = __builtin_ia32_vinsertf128_ps256(p, a, 0);
	p = __builtin_ia32_vinsertf128_ps256(p, b, 1);
	return p;
}

static inline v4df pack_2_v2df(const v2df a, const v2df b){
	v4df p = REP4(0.0f); // just avoid warning
	p = __builtin_ia32_vinsertf128_pd256(p, a, 0);
	p = __builtin_ia32_vinsertf128_pd256(p, b, 1);
	return p;
}

struct Particle{
	v4sf posH; // elemen w is for mass
	v4sf posL;
	v4sf vel;
	v4sf acc2;
	v4sf jrk6;
	v2df time; // 6 xmm words;

	void make_pos(const double *x, const double *m){
		const v4df xd = {x[0], x[1], x[2], m[0]};
		const v4sf xs = __builtin_ia32_cvtpd2ps256(xd);
		const v4df yd = xd - __builtin_ia32_cvtps2pd256(xs);
		const v4sf ys = __builtin_ia32_cvtpd2ps256(yd);
		posH = xs;
		posL = ys;
	}
	static v4sf make_v4sf(const double *x){
		const v4df xd = {x[0], x[1], x[2], 0.0};
		const v4sf xs = __builtin_ia32_cvtpd2ps256(xd);
		return xs;
	}
	Particle(
		const double _pos [3],
		const double _vel [3],
		const double _acc2[3],
		const double _jrk6[3],
		const double _mass,
		const double _time)
	{
		make_pos(_pos, &_mass);
		vel  = make_v4sf(_vel);
		acc2 = make_v4sf(_acc2);
		jrk6 = make_v4sf(_jrk6);
		time = (v2df){_time, _time};
	}

	Particle(int) // constructor for a dummy particle
	{
		posH = (v4sf){255.0f, 255.0f, 255.0f, 0.0f};
		posL = vel = acc2 = jrk6 = (v4sf)REP4(0.0f);
		time = (v2df){0.0, 0.0};
	}

	void prefetch(const int rw=0, const int locality=3) const {
		const char *addr = (char *)&posH;
		__builtin_prefetch( 0+addr, rw , locality);
		__builtin_prefetch(64+addr, rw , locality);
	}
	bool is_aligned() const {
		unsigned long addr = (unsigned long)(&posH);
		return (0 == addr%32);
	}
	void store(Particle *pdst) const {
		const v8sf ym0 = pack_2_v4sf(posH, posL);
		const v8sf ym1 = pack_2_v4sf(vel,  acc2);
		const v8sf ym2 = pack_2_v4sf(jrk6, (v4sf)time);
		v8sf *dst = (v8sf *)pdst;
		dst[0] = ym0;
		dst[1] = ym1;
		dst[2] = ym2;
	}
	void sstore(Particle *pdst) const {
		const v8sf ym0 = pack_2_v4sf(posH, posL);
		const v8sf ym1 = pack_2_v4sf(vel,  acc2);
		const v8sf ym2 = pack_2_v4sf(jrk6, (v4sf)time);
		v8sf *dst = (v8sf *)pdst;
		__builtin_ia32_movntps256((float *)(dst + 0), ym0);
		__builtin_ia32_movntps256((float *)(dst + 1), ym1);
		__builtin_ia32_movntps256((float *)(dst + 2), ym2);
	}
};

struct Pred2{ // couple of 2 predictors
	v8sf posH; // |m|z|y|x|m|z|y|x|
	v8sf posL; // |*|z|y|x|*|z|y|x|
	v8sf vel;  // |*|w|v|u|*|w|v|u|

	static inline v8sf gen_dt01(const v2df t0, const v2df t1, const v2df tnow){
		const v2df dt01d = tnow - __builtin_ia32_unpcklpd(t0, t1);
		const v4sf dt01s = __builtin_ia32_cvtpd2ps(dt01d);
		const v4sf dt0   = __builtin_ia32_shufps(dt01s, dt01s, 0x00);
		const v4sf dt1   = __builtin_ia32_shufps(dt01s, dt01s, 0x55);
		return pack_2_v4sf(dt0, dt1);
	}

	Pred2(const Particle &p0, const Particle &p1, const v2df tnow){
		const v8sf s0 = gen_dt01(p0.time, p1.time, tnow);
		const v8sf s1 = s0 + s0;
		const v8sf s2 = s0 * (v8sf)REP8(1.5f);

		const v8sf pos8 = pack_2_v4sf(p0.posL, p1.posL);
		const v8sf vel8 = pack_2_v4sf(p0.vel , p1.vel );
		const v8sf acc8 = pack_2_v4sf(p0.acc2, p1.acc2);
		const v8sf jrk8 = pack_2_v4sf(p0.jrk6, p1.jrk6);

		this->posH = pack_2_v4sf(p0.posH, p1.posH);
		// this->posH = pack_2_v4sf(p0.posH, p1.posH, this->posH);
		this->posL = pos8 + s0*(vel8 + s0*(acc8 + s0*(jrk8)));
		this->vel  = vel8 + s1*(acc8 + s2*(jrk8));
	}
};

struct Pred8{ // almost an array of structures
	v8sf xH, yH, zH, m;
	v8sf xL, yL, zL;
	v8sf vx, vy, vz;

	static void transpos(v8sf &v0, v8sf &v1, v8sf &v2, v8sf &v3){
		// transpose from
		// [ |w4|z4|y4|x4||w0|z0|y0|x0|, |w5|z5|y5|x5||w1|z1|y1|x1|, |w6|z6|y6|x6||w2|z2|y2|x2|, |w7|z7|y7|x7||w3|z3|y3|x3| ]
		// to
		// [ |x7|x6|x5|x4||x3|x2|x1|x0|, |y7|y6|y5|y4||y3|y2|y1|y0|, |z7|z6|z5|z4||z3|z2|z1|z0|, |w7|w6|w5|w4||w3|w2|w1|w0| ] 
		const v8sf y2y0x2x0 = __builtin_ia32_unpcklps256(v0, v2);
		const v8sf w2w0z2z0 = __builtin_ia32_unpckhps256(v0, v2);
		const v8sf y3y1x3x1 = __builtin_ia32_unpcklps256(v1, v3);
		const v8sf w3w1z3z1 = __builtin_ia32_unpckhps256(v1, v3);
		const v8sf xxxx = __builtin_ia32_unpcklps256(y2y0x2x0, y3y1x3x1);
		const v8sf yyyy = __builtin_ia32_unpckhps256(y2y0x2x0, y3y1x3x1);
		const v8sf zzzz = __builtin_ia32_unpcklps256(w2w0z2z0, w3w1z3z1);
		const v8sf wwww = __builtin_ia32_unpckhps256(w2w0z2z0, w3w1z3z1);
		v0 = xxxx;
		v1 = yyyy;
		v2 = zzzz;
		v3 = wwww;
	}

	Pred8(){}

	// for j-particle
	Pred8(const Pred2 &p0, const Pred2 &p1, const Pred2 &p2, const Pred2 &p3){
		{
			v8sf a0 = p0.posH;
			v8sf a1 = p1.posH;
			v8sf a2 = p2.posH;
			v8sf a3 = p3.posH;
			transpos(a0, a1, a2, a3);
			xH = a0;
			yH = a1;
			zH = a2;
			m  = a3;
		}
		{
			v8sf b0 = p0.posL;
			v8sf b1 = p1.posL;
			v8sf b2 = p2.posL;
			v8sf b3 = p3.posL;
			transpos(b0, b1, b2, b3);
			xL = b0;
			yL = b1;
			zL = b2;
		}
		{
			v8sf c0 = p0.vel;
			v8sf c1 = p1.vel;
			v8sf c2 = p2.vel;
			v8sf c3 = p3.vel;
			transpos(c0, c1, c2, c3);
			vx = c0;
			vy = c1;
			vz = c2;
		}
	}
	// for j-particle
	Pred8(const Pred2 &pr){
		xH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0x00);
		yH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0x55);
		zH = __builtin_ia32_shufps256(pr.posH, pr.posH, 0xaa);

		xL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0x00);
		yL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0x55);
		zL = __builtin_ia32_shufps256(pr.posL, pr.posL, 0xaa);

		vx = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0x00);
		vy = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0x55);
		vz = __builtin_ia32_shufps256(pr.vel,  pr.vel , 0xaa);
	}
};

// single precision version
struct Force{
	v8sf ax, ay, az;
	v8sf jx, jy, jz;

	void clear(){
		const v8sf zero = REP8(0.0f);
		ax = ay = az = zero;
		jx = jy = jz = zero;
	}
	static double reduce(const v8sf v){
		const v4df d0 = __builtin_ia32_cvtps2pd256(
		                  __builtin_ia32_vextractf128_ps256(v, 0));
		const v4df d1 = __builtin_ia32_cvtps2pd256(
		                  __builtin_ia32_vextractf128_ps256(v, 1));
		const v4df dd = d0 + d1;
		const v4df dh = __builtin_ia32_haddpd256(dd, dd);
		const v2df sum = __builtin_ia32_vextractf128_pd256(dh, 0)
		               + __builtin_ia32_vextractf128_pd256(dh, 1);
		return __builtin_ia32_vec_ext_v2df(sum, 0);
	}
	void write(double *acc, double *jrk) const{
#if 1
		const double a0 = reduce(ax);
		const double a1 = reduce(ay);
		const double a2 = reduce(az);
		const double j0 = reduce(jx);
		const double j1 = reduce(jy);
		const double j2 = reduce(jz);
		acc[0] = a0;
		acc[1] = a1;
		acc[2] = a2;
		jrk[0] = j0;
		jrk[1] = j1;
		jrk[2] = j2;
#else
		acc[0] = reduce(ax);
		acc[1] = reduce(ay);
		acc[2] = reduce(az);
		jrk[0] = reduce(jx);
		jrk[1] = reduce(jy);
		jrk[2] = reduce(jz);
#endif
	}
	void calc_and_accum(const Pred8 &pi, const Pred8 &pj){
		const v8sf dx = (pj.xH - pi.xH) + (pj.xL - pi.xL);
		const v8sf dy = (pj.yH - pi.yH) + (pj.yL - pi.yL);
		const v8sf dz = (pj.zH - pi.zH) + (pj.zL - pi.zL);

		const v8sf dvx = pj.vx - pi.vx;
		const v8sf dvy = pj.vy - pi.vy;
		const v8sf dvz = pj.vz - pi.vz;

		const v8sf r2 = dx*dx  + dy*dy  + dz*dz;
		const v8sf rv = dx*dvx + dy*dvy + dz*dvz;
		const v8sf rinv   = rsqrt_NR(r2);
		const v8sf rinv2  = rinv * rinv;
		const v8sf c1     = REP8(-3.0f);
		const v8sf alpha  = c1 * rinv2 * rv;
		const v8sf mrinv3 = pj.m * rinv * rinv2;

		ax += mrinv3 * dx;
		ay += mrinv3 * dy;
		az += mrinv3 * dz;

		jx += mrinv3 * (dvx + alpha * dx); 
		jy += mrinv3 * (dvy + alpha * dy); 
		jz += mrinv3 * (dvz + alpha * dz); 
	}
};

struct NBlist{
	// enum{ NB_MAX = 512 + 88 };
	enum{ NB_MAX = 768 + 32 };
	int pad[3];
	int nnb;
	int nb[NB_MAX];

	NBlist(){
		nnb = 0;
		for(int i=0; i<NB_MAX; i++){
			nb[i] = 0;
		}
	}

	void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

// The irregular force library
static bool       is_open = false;
static Particle  *ptcl = NULL;
static NBlist    *list = NULL;
static int        nmax;
static int        num_threads;
static double     time_grav;
static unsigned long long num_inter, num_fcall, num_steps;
static v2df vec_tnow;

static void gpuirr_open(
		const int nmax,
		const int lmax)
{
	if(is_open){
		fprintf(stderr, "gpuirr: it is already open\n");
		return;
	}

	assert(lmax <= 1 + NBlist::NB_MAX);

	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Opening GPUIRR lib. AVX ver. \n"); 
	fprintf(stderr, " nmax = %d, lmax = %d\n", nmax, lmax);
	// fprintf(stderr, " sizeof(NBlist) = %d\n", int(sizeof(NBlist))); 
	fprintf(stderr, "**************************** \n"); 
	assert(0 == sizeof(NBlist)%16);

	void *ptr;
	assert( 0 == posix_memalign(&ptr, 64, (1+nmax) * sizeof(Particle)) );
	memset(ptr, 0xff, (1+nmax) * sizeof(Particle));
	ptcl = (Particle *)ptr;

	// store dummy predictor to the last of the array
	ptcl[nmax] = Particle(0);

	// list.resize(nmax);
	assert( 0 == posix_memalign(&ptr, 64, nmax * sizeof(NBlist)) );
	memset(ptr, 0xff, nmax * sizeof(NBlist));
	list = (NBlist *)ptr;

	::nmax = nmax;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}
	time_grav = 0.0;
	num_inter = num_fcall = num_steps = 0;

	is_open = true;
}

static void gpuirr_close(){
	if(!is_open){
		fprintf(stderr, "gpuirr: it is already close\n");
		return;
	}

	free(ptcl); ptcl = NULL;
	free(list); list = NULL;

	const double Gflops = 60.0 * double(num_inter) * 1.e-9 / time_grav;
	const double usec_fcall    = 1.e6 * (time_grav / num_fcall);
	const double nnb_avr = double(num_inter) / double(num_steps);

	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Closing GPUIRR lib. AVX ver. \n"); 
	fprintf(stderr, "time grav  : %f sec\n", time_grav);
	fprintf(stderr, "\n");
	fprintf(stderr, "perf grav  : %f Gflops\n", Gflops);
	fprintf(stderr, "perf grav  : %f usec\n", usec_fcall);
	fprintf(stderr, "<#NB>      : %f \n",     nnb_avr);
	fprintf(stderr, "**************************** \n"); 

	is_open = false;
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
	Particle(pos, vel, acc2, jrk6, mass, time).store(&ptcl[addr]);
}

static void gpuirr_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
	assert(nnb <= NBlist::NB_MAX);
	list[addr].nnb = nnb;
#if 0
	for(int k=0; k<nnb; k++){
		list[addr].nb[k] = nblist[k] - 1;
	}
#else
	const int *src = nblist;
	      int *dst = list[addr].nb;
	for(int k=0; k<nnb; k+=8){
#if 0
		const int i0 = src[k+0] - 1;
		const int i1 = src[k+1] - 1;
		const int i2 = src[k+2] - 1;
		const int i3 = src[k+3] - 1;
		const int i4 = src[k+4] - 1;
		const int i5 = src[k+5] - 1;
		const int i6 = src[k+6] - 1;
		const int i7 = src[k+7] - 1;
		dst[k+0] = i0;
		dst[k+1] = i1;
		dst[k+2] = i2;
		dst[k+3] = i3;
		dst[k+4] = i4;
		dst[k+5] = i5;
		dst[k+6] = i6;
		dst[k+7] = i7;
#else
		// assert((unsigned long)dst %16 == 0);
		typedef int       v4si __attribute__((vector_size(16)));
		typedef long long v2di __attribute__ ((__vector_size__ (16)));
		const v4si one = REP4(1);
		const v4si idx0 = (v4si)__builtin_ia32_loaddqu((const char *)(src+k+0));
		const v4si idx1 = (v4si)__builtin_ia32_loaddqu((const char *)(src+k+4));
		__builtin_ia32_movntdq((v2di *)(dst+k+0), (v2di)(idx0-one));
		__builtin_ia32_movntdq((v2di *)(dst+k+4), (v2di)(idx1-one));
#endif
	}
	// fill dummy
	const int nmax = ::nmax;
	const int kmax = 8 * (1 + (nnb-1)/8);
	for(int k=nnb; k<kmax; k++){
		dst[k] = nmax;
	}
#endif
}

static void gpuirr_pred_all(
		const int    js,
		const int    je,
		const double ti)
{
	::vec_tnow = (v2df){ti, ti};
}
static void gpuirr_pred_act(
		const int ni,
		const int addr[],
		const double ti)
{
	::vec_tnow = (v2df){ti, ti};
}

static void gpuirr_firr(
		const int addr,
		_out_ double accout[3],
		_out_ double jrkout[3])
{
	const Particle *ptcl = ::ptcl;
	const v2df      tnow = ::vec_tnow;
	const int  nnb  = list[addr].nnb;

	Force force;
	force.clear();

	const Pred8 pri(Pred2(ptcl[addr], ptcl[addr], tnow));
#if 1
	for(int k=0; k<nnb; k++){
		const int jaddr = list[addr].nb[k];
		ptcl[jaddr].prefetch();
	}
#endif
#if 1
	for(int k=0; k<nnb; k+=8){
		const int *jptr = &(list[addr].nb[k]);
		const Pred2 p01(ptcl[jptr[0]], ptcl[jptr[1]], tnow);
		const Pred2 p23(ptcl[jptr[2]], ptcl[jptr[3]], tnow);
		const Pred2 p45(ptcl[jptr[4]], ptcl[jptr[5]], tnow);
		const Pred2 p67(ptcl[jptr[6]], ptcl[jptr[7]], tnow);

		const Pred8 prj(p01, p23, p45, p67);

		force.calc_and_accum(pri, prj);
	}
#else
	Pred8 pbuf[NBlist::NB_MAX/8];
	for(int k=0; k<nnb; k+=8){
		const int *jptr = &(list[addr].nb[k]);
		const Pred2 p01(ptcl[jptr[0]], ptcl[jptr[1]], tnow);
		const Pred2 p23(ptcl[jptr[2]], ptcl[jptr[3]], tnow);
		const Pred2 p45(ptcl[jptr[4]], ptcl[jptr[5]], tnow);
		const Pred2 p67(ptcl[jptr[6]], ptcl[jptr[7]], tnow);

		pbuf[k/8] =  Pred8(p01, p23, p45, p67);
	}
	for(int k=0; k<nnb; k+=8){
		const Pred8 prj = pbuf[k/8];
		force.calc_and_accum(pri, prj);
	}
#endif

	force.write(accout, jrkout);
}

static void gpuirr_firr_vec(
		const int ni,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3])
{
	const double t0 = get_wtime();
	int ninter = 0;
#if 1
#pragma omp parallel for reduction(+: ninter) schedule(guided)
	for(int i=0; i<ni; i++){
		gpuirr_firr(addr[i]-1, accout[i], jrkout[i]);
		ninter += list[addr[i]-1].nnb;
	}
#else // cache friendlly mode
#pragma omp parallel reduction(+: ninter)
	{
		const int nth = omp_get_num_threads();
		const int tid = omp_get_thread_num();
		for(int i=0; i<ni; i++){
			const int iaddr = addr[i]-1;
			if(tid == iaddr % nth){
				gpuirr_firr(iaddr, accout[i], jrkout[i]);
				ninter += list[iaddr].nnb;
			}
		}
	}
#endif
	::num_inter += ninter;
	const double t1 = get_wtime();
	::time_grav += t1-t0;
	::num_fcall++;
	::num_steps += ni;
}

// FORTRAN interface
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
	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3])
	{
		gpuirr_firr_vec(*ni, addr, acc, jrk);
	}
}

// test codes for the compiler
#if 0
extern void make_particle(const double *x, Particle *p){
	Particle(x, x+3, x+6, x+9, x[12], x[13]).sstore(p);
}


extern void make_pred2(const Particle ptcl[], Pred2 &pr, const v2df tnow){
	pr = Pred2(ptcl[0], ptcl[1], tnow);
}

extern void make_pred8(const Pred2 pr[], Pred8 &pr8){
	pr8 = Pred8(pr[0], pr[1], pr[2], pr[3]);
}

extern void make_pred8_from_list(
		const int      list[], 
		const Particle ptcl[], 
		const v2df     tnow,
		Pred8         &pr8)
{ 
	const Pred2 pr01 = Pred2(ptcl[list[0]], ptcl[list[1]], tnow);
	const Pred2 pr23 = Pred2(ptcl[list[2]], ptcl[list[3]], tnow);
	const Pred2 pr45 = Pred2(ptcl[list[4]], ptcl[list[5]], tnow);
	const Pred2 pr67 = Pred2(ptcl[list[6]], ptcl[list[7]], tnow);
	pr8 = Pred8(pr01, pr23, pr45, pr67);
}

extern void write_force(const Force &f, double *acc, double *jrk){
	f.write(acc, jrk);
}

extern void test_gravity(Force &f, const Pred8 pr[], const int N){
	f.clear();
	Pred8 pi = pr[N];
	for(int j=0; j<N; j++){
		f.calc_and_accum(pi, pr[j]);
	}
}
#endif
