
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_RECTANGLE_H
#define PBRT_SHAPES_RECTANGLE_H

// shapes/rectangle.h*
#include "shape.h"

namespace pbrt {

struct SphRectangleData {
  Vector3f x, y, z;
  Float ex1, ey1;
  Point3f o,s;
  Float z0, z0sq;
  Float x0, y0, y0sq;
  Float x1, y1, y1sq;
  Float b0, b1, b0sq, k;
  Float solidAngle;
};

void ComputeSphRectangleData(SphRectangleData *data, const Point3f &o);

// Rectangle Declarations
class Rectangle : public Shape {
  public:
    // Rectangle Public Methods
    Rectangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
         bool reverseOrientation, Float width, Float height, int samplingMode);
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    Interaction Sample(const Interaction &ref, const Point2f &u,
                       Float *pdf) const;
    Float Pdf(const Interaction &ref, const Vector3f &wi) const;

  private:
    // Rectangle Private Data
    const Float width, height;
    //samplingMode=1: muestreo uniforme respecto al área
    //samplingMode=2: muestreo uniforme respecto al ángulo sólido
    const int samplingMode;
};

std::shared_ptr<Rectangle> CreateRectangleShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_Rectangle_H
