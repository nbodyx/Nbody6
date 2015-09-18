#define _out_
namespace gninit{
	void count_neib(
			const int    nbody,
			const double eps2,
			const double pos[][3],
			const double h2 [],
			_out_ int    nnb[]);
	void get_neib(
			const int    nbody,
			const double eps2,
			const double pos[][3],
			const double h2 [],
			const int    nnb[],
			_out_ int * const list[]); // array of pointers
	void calc_pot(
			const int    nbody,
			const double eps2,
			const double mass[],
			const double pos [][3],
			_out_ double pot []);
	void calc_force(
			const int    nbody,
			const double eps2,
			const double mass[],
			const double pos [][3],
			const double vel [][3],
			_out_ double acc [][3],
			_out_ double jrk [][3]);
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
			_out_ double fpoly_irr[][4][3]);
}
