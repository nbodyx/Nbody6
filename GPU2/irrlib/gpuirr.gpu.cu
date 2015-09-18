#include <algorithm>
#include <cstdio>
#include <cutil.h>
#include "cuda_pointer.h"
#include "gnutil.h"

#define _out_
#define PROFILE

typedef DblFloat twofl;
typedef Gvec3<twofl> dvec3;
typedef Gvec3<float> fvec3;

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

struct Jparticle{
	dvec3 pos;  // 6
	fvec3 vel;  // 9
	fvec3 acc2; // 12
	fvec3 jrk6; // 15
	float mass; // 16
	twofl time; // 18
	int   addr; // 19
	int   pad;  // 20

	__host__ Jparticle(
			const double _pos [3],
			const double _vel [3],
			const double _acc2[3],
			const double _jrk6[3],
			const double _mass,
			const double _time,
			const int    _addr) :
		pos(_pos), vel(_vel), acc2(_acc2), jrk6(_jrk6), mass(_mass), time(_time), addr(_addr) {}

	__device__ Jparticle(const void *vsrc){
		const float4 *src = (float4 *)vsrc;
		_out_ float4 *dst = (float4 *)this;
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		dst[4] = src[4];
	}
	__device__ void store(void *vdst) const{
		const float4 *src = (float4 *)this;
		_out_ float4 *dst = (float4 *)vdst;
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
		dst[3] = src[3];
		dst[4] = src[4];
	}
};

struct Predictor{
	dvec3 pos;   // 6
	fvec3 vel;   // 9
	float mass;  // 10
	twofl tlast; // 12, just for padding

	__device__ Predictor(const Jparticle &p, const twofl &time)
	{
		pos   = p.pos;
		vel   = p.vel;
		mass  = p.mass;
		// tlast = time;

		const float dt = time - p.time;
		pos += dt*(p.vel + dt*(p.acc2 + dt*(p.jrk6)));
		vel += (2.f*dt)*(p.acc2 + (1.5f*dt)*(p.jrk6));
	}

	__host__ void print(const int i, FILE *fp = stdout) const{
		double x = pos.x;
		double y = pos.y;
		double z = pos.z;
		double vx = vel.x;
		double vy = vel.y;
		double vz = vel.z;
		fprintf(fp, "%4d : %f %f %f, %f %f %f, %f\n",
				i, x, y, z, vx, vy, vz, mass);
		fflush(fp);
	}
};

struct Force{
	fvec3 acc; // 3
	fvec3 jrk; // 6

	__device__ Force() : acc(0.f), jrk(0.f) {}
	__device__ Force(const Predictor &ip, const Predictor &jp)
	{
		const fvec3 dr = jp.pos - ip.pos;
		const fvec3 dv = jp.vel - ip.vel;
		const float r2 = dr * dr;
		const float rv = dr * dv;

		const float rinv   = rsqrt(r2);
		const float rinv2  = rinv * rinv;
		const float alpha  = -3.f * rv * rinv2;
		const float mrinv3 = jp.mass * rinv * rinv2;

		acc = mrinv3 * dr;
		jrk = mrinv3 * (dv + alpha * dr);
	}

	__device__ void operator += (const Force &rhs){
		acc += rhs.acc;
		jrk += rhs.jrk;
	}
};

struct NBlist{
	enum{ NB_MAX = 511 };
	int nnb;
	int nb[NB_MAX]; // 2 kBi

	__host__ void print(const int i, FILE *fp = stdout) const{
		fprintf(fp, "%6d%6d :", i, nnb);
		for(int k=0; k<nnb; k++){
			fprintf(fp, " %d", nb[k]);
		}
		fprintf(fp, "\n");
		fflush(fp);
	}
};

__global__ void kernel_jp_flush(
		const int nj,
		const Jparticle jpsrc[],
		_out_ Jparticle jpdst[])
{
#if 0
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(tid < nj){
		const Jparticle jp = jpsrc[tid];
		jpdst[jp.addr] = jp;
	}
#else
	// 5 threads   / particle
	// 6 particles / warp
	__shared__ float4 f4share[64];
	const int tid = threadIdx.x;
	const int gid = threadIdx.x + blockDim.x * blockIdx.x;
	const int wbid = gid / 32;
	const int wtid = gid % 32;
	const int pid = wtid / 5;
	const int qid = wtid % 5;
	const int iaddr = 6 * wbid + pid;
	if((iaddr < nj) && (pid < 6)){
		float4 *src = (float4 *)(jpsrc + iaddr);
		f4share[tid] = src[qid];
		const int paddr = 32*(tid/32) + 5*pid;
		const Jparticle &jp = *(Jparticle *)(f4share + paddr);
		const int jaddr = jp.addr;
		float4 *dst = (float4 *)(jpdst + jaddr);
		dst[qid] = f4share[tid];
	}
#endif
}

__global__ void kernel_list_flush(
		const int4   task  [],
		const int    list  [],
		_out_ NBlist nblist[])
{
	const int  bid    = blockIdx.x;
	const int  tid    = threadIdx.x;
	const int4 mytask = task[bid];
	const int  addr   = mytask.x;
	const int  nnb    = mytask.y;
	const int  off    = mytask.z;
	const int *src    = list + off;
	_out_ int *dst    = nblist[addr].nb;

	if(tid == 0) nblist[addr].nnb = nnb; // store

	for(int k=0; k<nnb; k+=blockDim.x){
		if(k+tid < nnb) dst[k+tid] = src[k+tid] - 1;
	}
}

__global__ void kernel_predict(
		const int       js,
		const int       je,
		const twofl     ti,
		const Jparticle ptcl[],
		_out_ Predictor pred[])
{
#if 0
	const int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if(js<=tid && tid<je){
		pred[tid] = Predictor(ptcl[tid], ti);
	}
#else
	const int tid = threadIdx.x;
	const int off = blockDim.x * blockIdx.x;
	const int nth = blockDim.x;
	__shared__ float4 sbuf[128*5];
	Jparticle *sptcl = (Jparticle *)sbuf;
	Predictor *spred = (Predictor *)sbuf;
	{   // LOAD
		float4 *src = (float4 *)(ptcl + off);
		float4 *dst = (float4 *)(sptcl);
#pragma unroll
		for(int k=0; k<5; k++, src+=nth, dst+=nth){
			dst[tid] = src[tid];
		}
	}

	__syncthreads();
	Predictor pp(sptcl[tid], ti);
	__syncthreads();
	spred[tid] = pp;
	__syncthreads();

	{   // STORE
		float4 *src = (float4 *)(spred);
		float4 *dst = (float4 *)(pred + off);
#pragma unroll
		for(int k=0; k<3; k++, src+=nth, dst+=nth){
			dst[tid] = src[tid];
		}
	}
#endif
}

__global__ void kernel_gravity(
		const int       iaddr[],
		const NBlist    list [],
		const Predictor pred [],
		_out_ Force     force[])
{
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
	const int nth = blockDim.x;
	const int i   = iaddr[bid];
	const int nnb = list[i].nnb;
	const int *nblist = list[i].nb;

	__shared__ Force fshare[64];

	const Predictor ip = pred[i];
	Force fo; // flushed by constructor
	
	for(int k=0; k<nnb; k+=nth){
		if(k+tid < nnb){
			const int j = nblist[k+tid];
			const Predictor jp = pred[j];
			fo += Force(ip, jp);
		}
	}

	fshare[tid] = fo;

	if(nnb > 32){
		__syncthreads();
		if(tid < 32) fshare[tid] += fshare[tid + 32];
	}

	if(nnb > 16){
		if(tid < 16) fshare[tid] += fshare[tid + 16];
	}
	if(nnb > 8){
		if(tid <  8) fshare[tid] += fshare[tid +  8];
	}
	if(tid <  4) fshare[tid] += fshare[tid +  4];
	if(tid <  2) fshare[tid] += fshare[tid +  2];
	if(tid <  1) fshare[tid] += fshare[tid +  1];

	if(tid == 0){
		force[bid] = fshare[0];
	}
}

static double time_jp, time_list, time_pred, time_grav, time_onep;

struct Jparticle_que{
	enum{
		MAX_JP_QUE = (1<<14), // 1MBi
		NTHREAD    = 64,
	};
	typedef cudaPointer<Jparticle> pointer;

	int      count;
	pointer  buf;
	pointer *dst;

	Jparticle_que() : count(0), buf(), dst(NULL) {}
	void initialize(pointer * const _dst){
		count = 0;
		buf.allocate(MAX_JP_QUE);
		dst = _dst;
	}
	void finalize(){
		buf.free();
		count = 0;
		dst   = NULL;
	}
	void flush(){
		if(count > 0){
			const double t0 = get_wtime();

			buf.htod(count);
#if 0
			const int nblock = 1 + (count-1)/NTHREAD;
#else       // 5 threads   / particle
			// 6 particles / warp
			const int nwarp  = 1 + (count-1)/6;
			const int nblock = 1 + (nwarp-1)/(NTHREAD/32); 
#endif

			kernel_jp_flush <<<nblock, NTHREAD>>>
				(count, buf, (*dst));
			CUDA_SAFE_THREAD_SYNC();

			count = 0;

			const double t1 = get_wtime();
			::time_jp += t1-t0;
		}
	}
	void push(
		const int addr,
		const double pos [3],
		const double vel [3],
		const double acc2[3],
		const double jrk6[3],
		const double mass,
		const double time)
	{
		buf[count++] = Jparticle(pos, vel, acc2, jrk6, mass, time, addr);
		if(count >= MAX_JP_QUE) flush();
	}
};

struct List_que{
	enum{
		MAX_LIST_QUE = 4096,
		MAX_LIST_BUF = (1<<20), // 4 MBi
		NTHREAD      = 64,
	};
	int                  count;
	int                  list_len;
	cudaPointer<int>     list_buf;
	cudaPointer<int4>    task_buf; // addr, nnb, off, pad
	cudaPointer<NBlist> *dst;

	List_que() : count(0), list_len(0), list_buf(), task_buf(), dst(NULL) {}
	void initialize(cudaPointer<NBlist> * const _dst){
		dst      = _dst;
		count    = 0;
		list_len = 0;
		list_buf.allocate(MAX_LIST_BUF);
		task_buf.allocate(MAX_LIST_QUE);
	}
	void finalize(){
		dst      = NULL;
		count    = 0;
		list_len = 0;
		list_buf.free();
		task_buf.free();
	}
	void push(
			const int addr,
			const int nnb,
			const int nblist[])
	{
		if(list_len + nnb > MAX_LIST_BUF) flush();
		task_buf[count++] = make_int4(addr, nnb, list_len, 0);
		int *dst = &list_buf[list_len];
		for(int k=0; k<nnb; k++){
			dst[k] = nblist[k];
		}
		list_len += nnb;
		if(count >= MAX_LIST_QUE) flush();
	}
	void flush(){
		if(count > 0){
			const double t0 = get_wtime();

			list_buf.htod(list_len);
			task_buf.htod(count);

			const int nblock = count;
			kernel_list_flush <<<nblock, NTHREAD>>>
				(task_buf, list_buf, (*dst));
			CUDA_SAFE_THREAD_SYNC();

			count    = 0;
			list_len = 0;

			const double t1 = get_wtime();
			::time_list += t1-t0;
		}
	}
};

enum{
	NFORCE_MAX = 1024,
};

static cudaPointer<Jparticle> ptcl;
static cudaPointer<Predictor> pred;
static cudaPointer<Force    > force;
static cudaPointer<NBlist   > list;
static Jparticle_que          jpque;
static List_que               listque;
static cudaPointer<int      > iaddr;

static void gpuirr_open(
		const int nmax,
		const int lmax)
{
	assert(lmax <= 1 + NBlist::NB_MAX);

	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Opening GPUIRR lib. GPU ver. \n"); 
	fprintf(stderr, " nmax = %d, lmax = %d\n", nmax, lmax);
	fprintf(stderr, "**************************** \n"); 

	ptcl .allocate(nmax + 128);
	pred .allocate(nmax + 128);
	force.allocate(NFORCE_MAX);
	list .allocate(nmax);
	jpque  .initialize(&ptcl);
	listque.initialize(&list);
	iaddr.allocate(NFORCE_MAX);

	time_jp = time_list = time_pred = time_grav = time_onep = 0.0;
	fprintf(stderr, "Opened GPUIRR lib. GPU ver. \n"); 
}

static void gpuirr_close(){
	fprintf(stderr, "**************************** \n"); 
	fprintf(stderr, "Closing GPUIRR lib. GPU ver. \n"); 
	fprintf(stderr, "time jp    : %f sec\n", time_jp);
	fprintf(stderr, "time list  : %f sec\n", time_list);
	fprintf(stderr, "time pred  : %f sec\n", time_pred);
	fprintf(stderr, "time grav  : %f sec\n", time_grav);
	fprintf(stderr, "time onep  : %f sec\n", time_onep);
	fprintf(stderr, "**************************** \n"); 

	ptcl .free();
	pred .free();
	force.free();
	list .free();
	jpque  .finalize();
	listque.finalize();
	iaddr.free();
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
	jpque.push(addr, pos, vel, acc2, jrk6, mass, time);
}

static void gpuirr_set_list(
		const int addr,
		const int nnb,
		const int nblist[])
{
	assert(nnb <= NBlist::NB_MAX);
	listque.push(addr, nnb, nblist);
}

static void gpuirr_pred_all(
		const int    js,
		const int    je,
		const double ti){
	jpque.flush();

	const double t0 = get_wtime();

	const int nthread = 128;
	const int nblock  = 1 + (je-1)/nthread;
	kernel_predict <<<nblock, nthread>>>
		(js, je, twofl(ti), ptcl, pred);
	CUDA_SAFE_THREAD_SYNC();

	const double t1 = get_wtime();
	::time_pred += t1-t0;
}

static void gpuirr_firr_vec(
		const int nitot,
		const int addr[],
		_out_ double accout[][3],
		_out_ double jrkout[][3])
{
	listque.flush();

	const double t0 = get_wtime();

	for(int ii=0; ii<nitot; ii+=NFORCE_MAX){
		const int ni = nitot-ii < NFORCE_MAX ? nitot-ii : NFORCE_MAX;
		for(int i=0; i<ni; i++){
			iaddr[i] = addr[ii + i] - 1;
		}
		iaddr.htod(ni);

		kernel_gravity <<<ni, 64>>> (iaddr, list, pred, force);

		force.dtoh(ni);
		for(int i=0; i<ni; i++){
			force[i].acc.write(accout[ii + i]);
			force[i].jrk.write(jrkout[ii + i]);
		}
	}

	const double t1 = get_wtime();
	if(nitot == 1){
		::time_onep += t1-t0;
	}else{
		::time_grav += t1-t0;
	}
}


// FORTRAN interfface
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
#if 1
	void gpuirr_firr_(
			int    *addr,
			double  acc[3],
			double  jrk[3])
	{
		gpuirr_firr_vec(1, addr, (double (*)[3])acc, (double (*)[3])jrk);
	}
#endif
	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3])
	{
		gpuirr_firr_vec(*ni, addr, acc, jrk);
	}

	void DEBUG_list(const int count){
		listque.flush();
		list.dtoh(count);
		for(int i=0; i<count; i++){
			list[i].print(i);
		}
	}
	void DEBUG_pred(const int count){
		pred.dtoh(count);
		for(int i=0; i<count; i++){
			pred[i].print(i);
		}
	}
}
