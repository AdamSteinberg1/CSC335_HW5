//solves problem 1 from Dr. Pounds
#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <vector>
#include <string>

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

//returns a collection of points (t,w) that represent a discrete approximation to y(t)
//f(t,y) is the derivative of y with respect to t
//[a,b] is the interval over which t will vary
//h is the step size
//y0 is the initial condition y(0)
//uses Euler's Method
std::vector<std::pair<double,double>> eulersMethod(std::function<double(double,double)> f, double a, double b, double h, double y0)
{
  int n = (b-a)/h; ////there will be n+1 points in the result
  std::vector<std::pair<double,double>> points;
  points.reserve(n+1);
  double t = a;
  double w = y0;
  points.push_back({t,w});
  for(int i = 0; i < n; i++)
  {
    w += h*f(t,w);
    t += h;
    points.push_back({t,w});
  }
  return points;
}

void graph(std::vector<std::pair<double,double>> points1, std::vector<std::pair<double,double>> points2, double h, std::string filename)
{
  std::ofstream script("tmp.plt");
  script << "set terminal pngcairo\n"
         << "set output '" << filename << ".png'\n"
         << "set xlabel 't'\n"
         << "set ylabel 'y'\n"
         << "set yrange [-10:10]\n"
         << "set title 'h = " << h << "'\n"
         << "plot '-' title 'Eulers Method' with lines, '-' title 'Runge-Kutta' with lines\n";


  for(auto& point : points1)
  {
    script << point.first << '\t' << point.second << '\n';
  }
  script << "e\n";
  for(auto& point : points2)
  {
    script << point.first << '\t' << point.second << '\n';
  }
  script.close();

  system("gnuplot tmp.plt");
  system("rm tmp.plt");
}

int main()
{
  const double y0 = 5;

  auto dydt = [&](double t, double y){
    return exp(t) * sin(y);
  };

  //part a
  double a = 0;
  double b = 3;
  double h = 0.1;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_a1");
  h = 0.01;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_a2");
  h = 0.001;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_a3");
  h = 0.0001;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_a4");


  //part b
  b = 10;
  h = 0.1;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_b1");
  h = 0.01;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_b2");
  h = 0.001;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_b3");
  h = 0.0001;
  graph(eulersMethod(dydt, a, b, h, y0), rungeKutta(dydt, a, b, h, y0), h, "graph_b4");
}
