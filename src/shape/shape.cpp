#include <string>

#include "shape.hpp"

Shape::Shape(std::string identifier)
{
  this->identifier = identifier;
}

unsigned int Shape::GetColor(void)
{
  return color;
}
