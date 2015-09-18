#include <cstdio>
#include <cassert>

#define CHKALIGN(ptr, num) assert((unsigned long)(ptr) % num == 0);
#define _out_ 

void cxvpred(
		const int js,
		const int je,
		const double time,
		const double t0   [],
		const double x0   [][3],
		const double x0dot[][3],
		const double f    [][3],
		const double fdot [][3],
		_out_ double x    [][3],
		_out_ double xdot [][3],
		_out_ double tpred[])
{
	// fprintf(stderr, "CXVPRED : js = %d, je = %d\n", js, je);
	CHKALIGN(js, 2);
	CHKALIGN(je, 2);
#pragma omp parallel for
	for(int j=js; j<je; j++){
		const double s = time - t0[j];
		const double s1 = 1.5 * s;
		const double s2 = 2.0 * s;
		for(int k=0; k<3; k++){
			x   [j][k] = x0   [j][k] + s *(x0dot[j][k] + s *(f   [j][k] + s*(fdot[j][k])));
			xdot[j][k] = x0dot[j][k] + s2*(f    [j][k] + s1*(fdot[j][k]));
			tpred[j] = time;
		}
	}
}

void cxvpred_sse(
		const int js,
		const int je,
		const double time,
		const double t0   [],
		const double x0   [][3],
		const double x0dot[][3],
		const double f    [][3],
		const double fdot [][3],
		_out_ double x    [][3],
		_out_ double xdot [][3],
		_out_ double tpred[])
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

	typedef double v2df __attribute__((vector_size(16)));
	struct dbl6{
		v2df a, b, c;
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
		void store(double *ptr) const{
			__builtin_ia32_movntpd(ptr+0, a);
			__builtin_ia32_movntpd(ptr+2, b);
			__builtin_ia32_movntpd(ptr+4, c);
		}
		const dbl6 operator+(const dbl6 &rhs) const{
			return dbl6(a+rhs.a, b+rhs.b, c+rhs.c);
		}
		const dbl6 operator*(const dbl6 &rhs) const{
			return dbl6(a*rhs.a, b*rhs.b, c*rhs.c);
		}
	};

	const v2df v2_time = {time, time};
#pragma omp parallel for
	for(int j=js; j<je; j+=2){
		const v2df v2_s = v2_time - *(v2df*)(t0+j);

		const dbl6 Hs (v2_s);
		const dbl6 Hs1((v2df){1.5, 1.5} * v2_s);
		const dbl6 Hs2((v2df){2.0, 2.0} * v2_s);

		const dbl6 Hx0   (x0   [j]);
		const dbl6 Hx0dot(x0dot[j]);
		const dbl6 Hf    (f    [j]);
		const dbl6 Hfdot (fdot [j]);

		const dbl6 Hx    = Hx0    + Hs *(Hx0dot + Hs *(Hf + Hs*(Hfdot)));
		const dbl6 Hxdot = Hx0dot + Hs2*(Hf     + Hs1*(Hfdot));

		Hx   .store(x   [j]);
		Hxdot.store(xdot[j]);
		__builtin_ia32_movntpd(&tpred[j], v2_time);
	}
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
	cxvpred_sse(cjs, cje, *time, t0, x0, x0dot, f, fdot, x, xdot, tpred);
}
