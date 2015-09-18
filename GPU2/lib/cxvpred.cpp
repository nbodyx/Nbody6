#include <cstdio>
#include <cassert>

#define CHKALIGN(ptr, num) assert((unsigned long)(ptr) % num == 0);
#define _out_ 
#define inout 
#define RCPTR(x) (* __restrict const x)

void cxvpred(
		const int js,
		const int je,
		const double time,
		const double RCPTR(t0   ),
		const double RCPTR(x0   )[3],
		const double RCPTR(x0dot)[3],
		const double RCPTR(f    )[3],
		const double RCPTR(fdot )[3],
		_out_ double RCPTR(x    )[3],
		_out_ double RCPTR(xdot )[3],
		inout double RCPTR(tpred))
{
	// fprintf(stderr, "CXVPRED : js = %d, je = %d\n", js, je);
	CHKALIGN(js, 2);
	CHKALIGN(je, 2);
#pragma omp parallel for
	for(int j=js; j<je; j++){
		if(tpred[j] == time) continue;
		const double s = time - t0[j];
		const double s1 = 1.5 * s;
		const double s2 = 2.0 * s;
		for(int k=0; k<3; k++){
			x   [j][k] = x0   [j][k] + s *(x0dot[j][k] + s *(f   [j][k] + s*(fdot[j][k])));
			xdot[j][k] = x0dot[j][k] + s2*(f    [j][k] + s1*(fdot[j][k]));
		}
		tpred[j] = time;
	}
}


typedef double v2df __attribute__((vector_size(16)));
struct dbl6{
	v2df a, b, c;
	dbl6() {}
	dbl6(const v2df _a, const v2df _b, const v2df _c) 
		: a(_a), b(_b), c(_c) {}
	dbl6(const double *ptr){
		a = *(v2df*)(ptr+0); 
		b = *(v2df*)(ptr+2); 
		c = *(v2df*)(ptr+4); 
	}
	dbl6(const v2df dt){
		a = __builtin_ia32_unpcklpd(dt, dt);
		b = dt;
		c = __builtin_ia32_unpckhpd(dt, dt);
	}
	void sstore(double *ptr) const{
		__builtin_ia32_movntpd(ptr+0, a);
		__builtin_ia32_movntpd(ptr+2, b);
		__builtin_ia32_movntpd(ptr+4, c);
	}
	void store(double *ptr) const{
		*(v2df *)(ptr+0) = a;
		*(v2df *)(ptr+2) = b;
		*(v2df *)(ptr+4) = c;
	}
#if 0
	const dbl6 operator+(const dbl6 &rhs) const{
		return dbl6(a+rhs.a, b+rhs.b, c+rhs.c);
	}
	const dbl6 operator*(const dbl6 &rhs) const{
		return dbl6(a*rhs.a, b*rhs.b, c*rhs.c);
	}
	friend dbl6 operator*(const double s, const dbl6 &x){
		if(2.0 == s){
			return dbl6(x.a+x.a, x.b+x.b, x.c+x.c);
		}else{
			v2df v = {s, s};
			return dbl6(v*x.a, v*x.b, v*x.c);
		}
	}
#endif
	const dbl6 &operator+=(const dbl6 &rhs){
		a += rhs.a;
		b += rhs.b;
		c += rhs.c;
		return *this;
	}
	const dbl6 &operator*=(const dbl6 &rhs){
		a *= rhs.a;
		b *= rhs.b;
		c *= rhs.c;
		return *this;
	}
	const dbl6 &operator*=(const double s){
		v2df ss = {s, s};
		a *= ss;
		b *= ss;
		c *= ss;
		return *this;
	}
};
struct mask2{
	v2df mask_bits;
	mask2(const v2df mask) : mask_bits(mask) {}
	mask2(const v2df a, const v2df b) :
		mask_bits((v2df)__builtin_ia32_cmpneqpd(a,b)) {}
	v2df select(const v2df x, const v2df y) const{
		return __builtin_ia32_orpd(
				__builtin_ia32_andpd (mask_bits, x),
				__builtin_ia32_andnpd(mask_bits, y));
	}
	dbl6 select(const dbl6 &x, const dbl6 &y) const{
		dbl6 mask6(mask_bits);
		return dbl6(
				mask2(mask6.a).select(x.a, y.a),
				mask2(mask6.b).select(x.b, y.b),
				mask2(mask6.c).select(x.c, y.c));
	}
};

void cxvpred_sse(
		const int js,
		const int je,
		const double time,
		const double RCPTR(t0   ),
		const double RCPTR(x0   )[3],
		const double RCPTR(x0dot)[3],
		const double RCPTR(f    )[3],
		const double RCPTR(fdot )[3],
		inout double RCPTR(x    )[3],
		inout double RCPTR(xdot )[3],
		inout double RCPTR(tpred))
{
	CHKALIGN(js, 2);
	CHKALIGN(je, 2);
	// fprintf(stderr, "CXVPRED : %p %p %p %p %p %p %p\n", 
	// 		t0, x0, x0dot, f, fdot, x, xdot);
	CHKALIGN(t0,    16);
	CHKALIGN(x0,    16);
	CHKALIGN(x0dot, 16);
	CHKALIGN(f,     16);
	CHKALIGN(fdot,  16);
	CHKALIGN(x,     16);
	CHKALIGN(xdot,  16);
	CHKALIGN(tpred, 16);
#pragma omp parallel for
	for(int j=js; j<je; j+=2){
#if 0
		__builtin_prefetch(&t0   [j+2]);
		__builtin_prefetch(&x0   [j+2]);
		__builtin_prefetch(&x0dot[j+2]);
		__builtin_prefetch(&f    [j+2]);
		__builtin_prefetch(&fdot [j+2]);
		__builtin_prefetch(&x    [j+2]);
		__builtin_prefetch(&xdot [j+2]);
		__builtin_prefetch(&tpred[j+2]);
#endif
		const v2df v2_time = {time, time};
		const v2df v2_s = v2_time - *(v2df*)(t0+j);

		const dbl6 Hx0   (x0   [j]);
		const dbl6 Hx0dot(x0dot[j]);
		const dbl6 Hf    (f    [j]);
		const dbl6 Hfdot (fdot [j]);
		const dbl6 Hs (v2_s);
#if 0
		const dbl6 Hx    = Hx0    + Hs *(Hx0dot + Hs *(Hf + Hs*(Hfdot)));

		const dbl6 Hs1((v2df){1.5, 1.5} * v2_s);
		const dbl6 Hs2((v2df){2.0, 2.0} * v2_s);
		const dbl6 Hxdot = Hx0dot + Hs2*(Hf     + Hs1*(Hfdot));
#else
		dbl6 Hx, Hxdot;
		Hx     = Hfdot;
		Hx    *= Hs;
		Hxdot  = Hx;
		Hxdot *= 1.5;
		Hx    += Hf;
		Hxdot += Hf;
		Hx    *= Hs;
		Hxdot *= Hs;
		Hxdot += Hxdot; // multiply by 2
		Hx    += Hx0dot;
		Hxdot += Hx0dot;
		Hx    *= Hs;
		Hx    += Hx0;
#endif

#if 0
		Hx   .store(x   [j]);
		Hxdot.store(xdot[j]);
		__builtin_ia32_movntpd(&tpred[j], v2_time);
#else
		const v2df v2_tpred = *(v2df *)(tpred+j);
		const mask2 m2(v2_time, v2_tpred);
		Hx    = m2.select(Hx   , dbl6(x   [j]));
		Hxdot = m2.select(Hxdot, dbl6(xdot[j]));
		Hx   .store(x   [j]);
		Hxdot.store(xdot[j]);
		*(v2df *)(tpred+j) = v2_time;
#endif
	}
#if 0
#pragma omp parallel for
	for(int j=js; j<je; j+=2){
		const v2df v2_time = {time, time};
		__builtin_ia32_movntpd(&tpred[j], v2_time);
	}
#endif
}

extern "C" void cxvpred_(
		const int *js,
		const int *je,
		const double *time,
		const double t0   [],
		const double x0   [][3],
		const double x0dot[][3],
		const double f    [][3],
		const double fdot [][3],
		_out_ double x    [][3],
		_out_ double xdot [][3],
		_out_ double tpred[])
{
	int cjs =*js;
	cjs--;
	int cje = *je;
	cje += cje%2;
        cxvpred(cjs, cje, *time, t0, x0, x0dot, f, fdot, x, xdot, tpred);
        // cxvpred_sse(cjs, cje, *time, t0, x0, x0dot, f, fdot, x, xdot, tpred);
}
