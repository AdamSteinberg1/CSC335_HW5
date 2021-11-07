//solves problem 5 from section 5.9
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
  const double k1 = 3;
  const double k2 = 0.002;
  const double k3 = 0.0006;
  const double k4 = 0.5;

  //derivative of x1 with respect to t
  auto x1dot =  [&](double t, double x1, double x2){
    return k1*x1 - k2*x1*x2;
  };

  //derivative of x2 with respect to t
  auto x2dot =  [&](double t, double x1, double x2){
    return k3*x1*x2 - k4*x2;
  };

  double initialPrey = 1000;
  double initialPredator = 500;

  double h = 0.01;
  double a = 0;
  double b = 4;
  auto approximation = rungeKutta(x1dot, x2dot, a, b, h, initialPrey, initialPredator);
  std::cout.precision(10);
  std::cout << "t\tx1(t)\tx2(t)" << std::endl;
  for(auto& elem: approximation)
    std::cout << elem[0] << '\t' << elem[1] << '\t' << elem[2] << std::endl;
}
