#ifndef shape_hpp
#define shape_hpp

#include <string>

class Shape
{
  private:
    std::string identifier;
    unsigned int color;
    unsigned int material;
  public:
    Shape(std::string identifier);
    unsigned int GetColor(void);
    std::string GetIdentifier(void);
    unsigned int GetMaterial(void);
    void SetColor(unsigned int color);
    void SetIdentifier(std::string identifier);
    void SetMaterial(unsigned int material);
};

#endif /* shape_hpp */
