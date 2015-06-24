#include <cmath>
#include <string>

#include "circle.hpp"

Circle::Circle(std::string identifier, double x, double y, double r)
{
  this->identifier = identifier;
  this->x = x;
  this->y = y;
  this->r = r;
}

double Circle::GetArea(void)
{
  return M_PI*r*r;
}

double Circle::GetPerimeter(void)
{
  return 2*M_PI*r;
}
