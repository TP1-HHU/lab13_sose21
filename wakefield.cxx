#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "vec.hxx"
//---------------------------------------
using namespace std;
//---------------------------------------
void init(Vec<double>& A0, Vec<double>& A1, Vec<double>& dn0, Vec<double>& dn1,
	        const double amp0, const double sigma, const double k0,
					const double om0, const double x0, const double dx, const double dt);
void stepA(Vec<double>& A2, const Vec<double>& A1, const Vec<double>& A0,
	         const Vec<double>& dn1, const double dt, const double dx);
void stepdn(Vec<double>& dn2, const Vec<double>& dn1, const Vec<double>& dn0,
	          const Vec<double>& A1, Vec<double>& h,
						const double dt,  const double dx);
//---------------------------------------
int main(){

  const int N = 16000;
  const double L = 2000;
  const double Tend = 2000;
  const double dx = L/(N-1.0);
  const double dt = 0.9*dx;
  const int Na = 10;
  int Nk = int(Tend/Na/dt + 0.5);

  const double amp0 = 0.2;
  const double sigma = 17;
  const double xmin = 0;
  const double x0 = 80;
  const double om0 = 2;
  const double k0 = sqrt(om0*om0 -1);

  Vec<double> A0(N);
  Vec<double> A1(N);
  Vec<double> A2(N);
  Vec<double> dn0(N);
  Vec<double> dn1(N);
  Vec<double> dn2(N);
  Vec<double> h(N);

  stringstream strm;

  init(A0,A1,dn0,dn1,amp0,sigma,k0,om0,x0,dx,dt);

  dn0.writeToFile( "dn_0", dx, xmin);
  A0.writeToFile( "A_0", dx, xmin);
  A1.writeToFile("a1_0", dx , xmin);

  cout << "Nk = " << Nk << endl;

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      stepdn(dn2,dn1,dn0,A1,h,dt,dx);
      stepA(A2,A1,A0,dn1,dt,dx);
      swap(A0,A1);
      swap(A1,A2);
      swap(dn0,dn1);
      swap(dn1,dn2);
   }
   strm.str("");
   strm << "A_" << i;
   A1.writeToFile(strm.str(), dx, xmin);
   strm.str("");
   strm << "dn_" << i;
   dn1.writeToFile(strm.str(), dx, xmin);
  }

  cout << dt << endl;
  return 0;
}
//-----------------------------------------------
void stepA(Vec<double>& A2, const Vec<double>& A1, const Vec<double>& A0,
	         const Vec<double>& dn1, const double dt, const double dx)
{
  int N = A0.size();
  const double dx2 = dx*dx;
  const double dt2 = dt*dt;

  A2(0) = 0;
  for(int i=1; i<N-1; i++){
	double rhs = -(1.0 + dn1(i) - 0.5*abs(A1(i))*abs(A1(i)) )* A1(i);
	double Axx = (A1(i+1)-2*A1(i)+A1(i-1))/dx2;
	A2(i) = (rhs + Axx)*dt2 +2*A1(i) -A0(i);
  }
  A2(N-1) = 0;

}
//-----------------------------------------------
void stepdn(Vec<double>& dn2, const Vec<double>& dn1, const Vec<double>& dn0,
	          const Vec<double>& A1,   Vec<double>& h,const double dt,  const double dx)
{
    int N = dn0.size();
    const double dx2 = dx*dx;
    const double dt2 = dt*dt;

    for(int i=0; i<N; i++){
	h(i) = abs(A1(i)) * abs(A1(i));
    }

    for(int i=1; i<N-1; i++){
       double rhs = 0.5 * (h(i+1) - 2*h(i) + h(i-1))/dx2;
       dn2(i) = (rhs - dn1(i))*dt2 - dn0(i) + 2*dn1(i);
    }

}
//-----------------------------------------------
void init(Vec<double>& A0, Vec<double>& A1, Vec<double>& dn0, Vec<double>& dn1,
	        const double amp0, const double sigma, const double k0,
					const double om0, const double x0, const double dx, const double dt)
{
	int N = dn0.size();
	const double v0 = k0 / om0;
	for(int i=0; i<N; i++){
            double x = i*dx;
	    dn0(i) = 0;
            dn1(i) = 0;

	    A0(i) = amp0 * exp(-pow( (x-x0)/sigma,2)) * sin(-k0*x);
            A1(i) = amp0 * exp(-pow( (x-x0+ v0*dt)/sigma,2))*sin(-k0*x-om0*dt);
	}
}
