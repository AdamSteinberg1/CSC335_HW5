//Solves problem 1 from section 5.9

#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

//4th order Runge-Kutta for a system of two differential equations
//returns a vector of triples (t,w1,w2) that represent a discrete approximation to u1(t) and u2(t)
std::vector<std::array<double,3>> rungeKutta( std::function<double(double,double,double)> f1, //derivative of u1
                                              std::function<double(double,double,double)> f2, //dericative of u2
                                              double a, double b, //interval over whcih t will vary
                                              double h, //step size
                                              double u1_0, double u2_0) //initial conditions
{
  int n = (b-a)/h; ////there will be n+1 points in the result
  std::vector<std::array<double,3>> result;
  result.reserve(n+1);
  double t = a;
  double w1 = u1_0;
  double w2 = u2_0;
  result.push_back({t,w1,w2});
  for(int i = 0; i < n; i++)
  {
    const double k1 = h*f1(t,w1,w2);
    const double j1 = h*f2(t,w1,w2);
    const double k2 = h*f1(t+h/2, w1+k1/2, w2+j1/2);
    const double j2 = h*f2(t+h/2, w1+k1/2, w2+j1/2);
    const double k3 = h*f1(t+h/2, w1+k2/2, w2+j2/2);
    const double j3 = h*f2(t+h/2, w1+k2/2, w2+j2/2);
    const double k4 = h*f1(t+h, w1+k3, w2+j3);
    const double j4 = h*f2(t+h, w1+k3, w2+j3);
    w1 += (k1 + 2*k2 + 2*k3 + k4)/6;
    w2 += (j1 + 2*j2 + 2*j3 + j4)/6;
    t += h;
    result.push_back({t,w1,w2});
  }
  return result;
}

int main()
{
  //derivative of u1 with respect to t
  auto f1 =  [](double t, double u1, double u2){
    return -4*u1 - 2*u2 + cos(t) + 4*sin(t);
  };
  double u1_0 = 0; //initial condition u1(0)=0

  //derivative of u2 with respect to t
  auto f2 =  [](double t, double u1, double u2){
    return 3*u1 + u2 - 3*sin(t);
  };
  double u2_0 = -1; //initial condition u2(0)=-1

  double h = 0.1;
  double a = 0;
  double b = 2;
  auto approximation = rungeKutta(f1, f2, a, b, h, u1_0, u2_0);
  std::cout.precision(5);
  std::cout << "t\tu1(t)\tu2(t)" << std::endl;
  for(auto& elem: approximation)
    std::cout << elem[0] << '\t' << elem[1] << '\t' << elem[2] << std::endl;
}
