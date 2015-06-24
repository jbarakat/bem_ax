#ifndef circle_hpp
#define circle_hpp

#include <string>
#include "shape.hpp"

class Circle: Shape
{
  private:
    std::string identifier;
    unsigned int color;
    unsigned int material;
    double x, y, r;
  public:
    Circle(std::string identifier);
    Circle(std::string identifier, double x, double y, double r);
    double GetArea(void);
    unsigned int GetColor(void);
    std::string GetIdentifier(void);
    unsigned int GetMaterial(void);
    double GetPerimeter(void);
    void SetCenter(double x, double y);
    void SetColor(unsigned int color);
    void SetIdentifier(std::string identifier);
    void SetMaterial(unsigned int material);
    void SetRadius(double r);
}

#endif /* circle_hpp */
