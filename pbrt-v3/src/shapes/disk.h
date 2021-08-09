
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_DISK_H
#define PBRT_SHAPES_DISK_H

// shapes/disk.h*
#include "shape.h"

namespace pbrt {

struct SphEllipseData {
  Vector3f xd, yd, zd;
  Vector3f oc;
  Vector3f ye, ze;
  int signX, signY;
  Float a, b, alpha, beta, at, bt;
  double solidAngle;
};

// Disk Declarations
class Disk : public Shape {
  public:
    // Disk Public Methods
    Disk(const Transform *ObjectToWorld, const Transform *WorldToObject,
         bool reverseOrientation, Float height, Float radius, Float innerRadius,
         Float phiMax, int samplingMode)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          height(height),
          radius(radius),
          innerRadius(innerRadius),
          phiMax(Radians(Clamp(phiMax, 0, 360))),
          samplingMode(samplingMode) {}
    Bounds3f ObjectBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    Float Area() const;
    Interaction Sample(const Point2f &u, Float *pdf) const;
    Interaction Sample(const Interaction &ref, const Point2f &u,
                       Float *pdf) const;
    void ComputeSphEllipseData(SphEllipseData *data, Point3f o) const;
    Point3f ReprojectToDisk(const Vector3f &q, const SphEllipseData *data, const Point3f &o) const;

  private:
    // Disk Private Data
    const Float height, radius, innerRadius, phiMax;
    const int samplingMode;
};

std::shared_ptr<Disk> CreateDiskShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_DISK_H
