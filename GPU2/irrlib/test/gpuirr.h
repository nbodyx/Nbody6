extern "C"{
	void gpuirr_open_(int *nmax, int *lmax);
	void gpuirr_close_();
	void gpuirr_set_jp_(
		int    *addr,
		double  pos [3],
		double  vel [3],
		double  acc2[3],
		double  jrk6[3],
		double *mass,
		double *time);
	void gpuirr_set_list_(
		int *addr,
		int *nblist);
	void gpuirr_pred_all_(
			int    *js,
			int    *je,
			double *ti);
#if 0
	void gpuirr_firr_(
			int    *addr,
			double  acc[3],
			double  jrk[3]);
#endif
	void gpuirr_firr_vec_(
			int   *ni,
			int    addr[],
			double acc [][3],
			double jrk [][3]);

	void DEBUG_list(const int count);
	void DEBUG_pred(const int count);
}
