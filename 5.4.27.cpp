//solves problem 27 from section 5.4
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

//returns a collection of points (t,w) that represent a discrete approximation to y(t)
//f(t,y) is the derivative of y with respect to t
//[a,b] is the interval over which t will vary
//there will be n+1 points in the result
//y0 is the initial condition y(0)
//uses 4th order Runge-Kutta
std::vector<std::pair<double,double>> rungeKutta(std::function<double(double,double)> f, double a, double b, int n, double y0)
{
  std::vector<std::pair<double,double>> points;
  points.reserve(n+1);
  double h = (b-a)/n;
  double t = a;
  double w = y0;
  points.push_back({t,w});
  for(int i = 0; i < n; i++)
  {
    double k1 = h*f(t,w);
    double k2 = h*f(t+h/2, w+k1/2);
    double k3 = h*f(t+h/2, w+k2/2);
    double k4 = h*f(t+h, w+k3);
    w += (k1 + 2*k2 + 2*k3 + k4)/6;
    t += h;
    points.push_back({t,w});
  }
  return points;
}



int main()
{
  const double n1 = 2e3;
  const double n2 = 2e3;
  const double n3 = 3e3;
  const double k = 6.22e-19;

  //dx/dt
  auto xdot = [&](double t, double x){
    return k * pow(n1-x/2, 2) * pow(n2-x/2, 2) * pow(n3-3*x/4, 3);
  };

  auto points = rungeKutta(xdot, 0, 0.2, 1000, 0);

  std::cout.precision(10);
  std::cout << "t\tx(t)" << std::endl;
  for(auto& point : points)
  {
    std::cout<< point.first << '\t' << point.second << std::endl;
  }
}
