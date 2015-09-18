#include <iostream>
#include <sys/time.h>
static double get_wtime(){
#if 1
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
#else
	struct timespec tv;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
	return tv.tv_sec + 1.e-9 * tv.tv_nsec;
#endif
}

static double t0 = get_wtime();

extern "C" {
	void wtime_(int *ihour, int *imin, int *isec){
		double t = get_wtime() - t0;
		int ih = int(t / 3600.);
		t -= ih * 3600.;
		int im = int(t / 60.);
		t -= im * 60.;
		int is = int(t);
		*ihour = ih;
		*imin = im;
		*isec = is;
	}
	double dbltim_(){
		return get_wtime() - t0;
	}
}

#if 0
#include <unistd.h>
int main(){
	int ih, im, is;
	sleep(1);
	wtime_(&ih, &im, &is);
	std::cout << ih << ":" << im << ":" << is << std::endl;
	return 0;
}
#endif
