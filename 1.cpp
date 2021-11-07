//solves problem 1 from Dr. Pounds
#include <iostream>
#include <fstream>
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

void graph(std::vector<std::vector<std::pair<double,double>>> solutions)
{
  std::vector<double> h = {0.2, 0.1, 0.05, 0.01, 0.001};
  std::ofstream script("tmp.plt");
  script << "set terminal pngcairo\n"
         << "set output 'graph1.png'\n"
         << "set xlabel 'x'\n"
         << "set ylabel 'y'\n"
         << "set xrange[0:4]\n"
         << "set yrange[-1.1:1.1]\n"
         << "y(x) = exp(-x)\n"
         << "plot y(x) title 'exp(-x)', ";
  for(int i = 0; i < solutions.size()-1; i++)
    script << "'-' title 'h=" << h[i] << "' with lines, ";
  script << "'-' title 'h=" << h.back() << "' with lines\n";

  for(auto& solution : solutions)
  {
    for(auto& point : solution)
    {
      script << point.first << '\t' << point.second << '\n';
    }
    script << "e\n";
  }
  script.close();

  system("gnuplot tmp.plt");
  system("rm tmp.plt");
}

int main()
{
  const double y0 = 1;
  const double a = 0;
  const double b = 4;

  auto dydx = [&](double x, double y){
    return 5*y - 6*exp(-x);
  };

  auto points1 = rungeKutta(dydx, a, b, 0.2, y0);
  auto points2 = rungeKutta(dydx, a, b, 0.1, y0);
  auto points3 = rungeKutta(dydx, a, b, 0.05, y0);
  auto points4 = rungeKutta(dydx, a, b, 0.01, y0);
  auto points5 = rungeKutta(dydx, a, b, 0.001, y0);

  graph({points1, points2, points3, points4, points5});
}
