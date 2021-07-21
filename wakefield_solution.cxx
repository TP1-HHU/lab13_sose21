#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

//---------------------------------------
using namespace std;
//-----------------------------------
//---------------------------------------
void init(double* const a0, double* const a1, double* const dn0, double* const dn1,
	        const double amp0, const double sigma, const double k0,
					const double om0, const double x0, const double dx, const double dt,
          const int N);
void stepA(double* const a2, const double* const a1, const double* const a0,
	         const double* dn1, const double dt, const double dx, const int N);
void stepdn(double* const dn2, const double* const dn1, const double* const dn0,
	          const double* const A1,  double* const h, const double dt,  const double dx,
            const int N);
void writeToFile(double* const psi, const string s, const double dx,
                 const int Nx, const double xmin);            
//---------------------------------------
int main(){

  const int N = 16000;
  const double L = 2000;
  const double Tend = 1900;
  const double dx = L/(N-1.0);
  const double dt = 0.9*dx;
  const int Na = 10;
  int Nk = int(Tend/Na/dt + 0.5);

  const double amp0 = 0.2;  // in m_e*c/e
  const double sigma = 17;   //in c/om_pe
  const double xmin = 0;
  const double x0 = 80; // in c/om_pe
  const double om0 = 2; //  in om_pe

  const double k0 = sqrt(om0*om0 -1);

  double* a0 = new double[N];
  double* a1 = new double[N];
  double* a2 = new double[N];
  
  double* dn0 = new double[N];
  double* dn1 = new double[N];
  double* dn2 = new double[N];
  
  double* h = new double[N];

  stringstream strm;

  init(a0,a1,dn0,dn1,amp0,sigma,k0,om0,x0,dx,dt,N);

  writeToFile(a1, "a_0", dx, N, xmin);
  writeToFile(dn1, "dn_0", dx, N, xmin);

  cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      stepdn(dn2,dn1,dn0,a1,h,dt,dx,N);
      stepA(a2,a1,a0,dn1,dt,dx,N);
      swap(a0,a1);
      swap(a1,a2);
      swap(dn0,dn1);
      swap(dn1,dn2);
   }
   strm.str("");
   strm << "a_" << i;
   writeToFile(a1, strm.str(), dx, N, xmin);
   strm.str("");
   strm << "dn_" << i;
   writeToFile(dn1, strm.str(), dx, N, xmin);

  }

  
  delete[] a0;
  delete[] a1;
  delete[] a2;

  delete[] dn0;
  delete[] dn1;
  delete[] dn2;

  delete[] h;

  return 0;
}
//-----------------------------------------------
void stepA(double* const a2, const double* const a1, const double* const a0,
	         const double* dn1, const double dt, const double dx, const int N)
{
  
  const double dx2 = dx*dx;
  const double dt2 = dt*dt;

  a2[0] = 0;
  for(int i=1; i<N-1; i++){
	  double rhs = -(1.0 + dn1[i] - 0.5*abs(a1[i])*abs(a1[i]) )* a1[i];
	  double Axx = (a1[i+1]-2*a1[i]+a1[i-1])/dx2;
	  a2[i] = (rhs + Axx)*dt2 +2*a1[i] -a0[i];
  }
  a2[N-1] = 0;

}
//-----------------------------------------------
void stepdn(double* const dn2, const double* const dn1, const double* const dn0,
	          const double* const a1,  double* const h, const double dt,  const double dx,
            const int N)
{
    const double dx2 = dx*dx;
    const double dt2 = dt*dt;

    for(int i=0; i<N; i++){
	      h[i] = abs(a1[i])*abs(a1[i]);
    }

    for(int i=1; i<N-1; i++){
       double rhs = 0.5 * (h[i+1] - 2*h[i] + h[i-1])/dx2;
       dn2[i] = (rhs - dn1[i])*dt2 - dn0[i] + 2*dn1[i];
    }

}
//-----------------------------------------------
void init(double* const a0, double* const a1, double* const dn0, double* const dn1,
	        const double amp0, const double sigma, const double k0,
					const double om0, const double x0, const double dx, const double dt,
          const int N)
{
	const double v0 = k0 / om0;
	for(int i=0; i<N; i++){
            double x = i*dx;
	          dn0[i] = 0;
            dn1[i] = 0;
	          a1[i] = amp0 * exp(-pow( (x-x0)/sigma,2)) * sin(-k0*x);
            a0[i] = amp0 * exp(-pow( (x-x0+ v0*dt)/sigma,2))*sin(-k0*x-om0*dt);
	}
}
//-----------------------------------
void writeToFile(double* const psi, const string s, const double dx,
                 const int Nx, const double xmin)
{
  ofstream out(s.c_str());

	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
    out << x << "\t" << psi[i] <<  endl;
	}
	out.close();
}