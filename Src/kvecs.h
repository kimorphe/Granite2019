#include <complex>
using namespace std;
class FSLICE{
	public:
		double freq;
		double Xa[2],Xb[2],Wd[2],dx[2];
		double **Kx,**Ky,**Phi,**Amp,**psi;
		double **Kxx,**Kyy,**Kxy;
		double k_mean,k_sig,a_mean,a_sig;
		int Nx,Ny,ndat,Nd[2];
		void init(int n, int m);
		void set_Xa(double *xa);
		void set_dx(double *dx);
		void set_Wd();
		void print_domain();
		void Grad();
		void Integrate();
		void Integrate2();
		void get_slice(complex<double> ***Z,int k);
		void export_phix(char fn[128]);
		void export_Grad(char fn[128]);
		void export_Hess(char fn[128]);
		void load_Grad(char fn[128]);
		void histogram(double kmin, double kmax, double nbin, double *prob_k, double *prob_a);
		int is_alloc;	// 0:No, 1:Yes
	private:
		void malloc_arrays();
		double** mem_alloc(int n, int m);
		int** mem_ialloc(int n, int m);
};
