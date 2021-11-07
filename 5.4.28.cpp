//solves problem 28 from section 5.4
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>

//returns a collection of points (t,w) that represent a discrete approximation to y(t)
//f(t,y) is the derivative of y with respect to t
//[a,b] is the interval over which t will vary
//h is the step size
//y0 is the initial condition y(0)
//uses 4th order Runge-Kutta
std::vector<std::pair<double,double>> rungeKutta(std::function<double(double,double)> f, double a, double b, double h, double y0)
{
  int n = (b-a)/h; ////there will be n+1 points in the result
  std::vector<std::pair<double,double>> points;
  points.reserve(n+1);
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
  const double r = 0.1;
  const double g = 32.1;
  const double x0 = 8;

  //dx/dt
  auto xdot = [&](double t, double x){
    return -0.6 * r*r * sqrt(2*g) * pow(x, -1.5);
  };

  //part a
  std::cout<< "Part a" << std::endl;

  auto points = rungeKutta(xdot, 0, 600, 20, 8);

  std::cout.precision(10);
  std::cout << "t\tx(t)" << std::endl;
  for(auto& point : points)
  {
    std::cout<< point.first << '\t' << point.second << std::endl;
  }

  //part b
  std::cout<< "\nPart b" << std::endl;
  points = rungeKutta(xdot, 0, 1500, 60, 8);
  std::cout << "t\tx(t)" << std::endl;
  for(auto& point : points)
  {
    std::cout<< point.first << '\t' << point.second << std::endl;
  }
}
