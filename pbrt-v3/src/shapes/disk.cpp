
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


// shapes/disk.cpp*
#include "shapes/disk.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"
#include "project/newtonRaphson.h"

namespace pbrt {

// Disk Method Definitions
Bounds3f Disk::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, height),
                    Point3f(radius, radius, height));
}

bool Disk::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                     bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute plane intersection for disk

    // Reject disk intersections for rays parallel to the disk's plane
    if (ray.d.z == 0) return false;
    Float tShapeHit = (height - ray.o.z) / ray.d.z;
    if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

    // See if hit point is inside disk radii and $\phimax$
    Point3f pHit = ray(tShapeHit);
    Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$
    Float phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;
    if (phi > phiMax) return false;

    // Find parametric representation of disk hit
    Float u = phi / phiMax;
    Float rHit = std::sqrt(dist2);
    Float v = (radius - rHit) / (radius - innerRadius);
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv =
        Vector3f(pHit.x, pHit.y, 0.) * (innerRadius - radius) / rHit;
    Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

    // Refine disk intersection point
    pHit.z = height;

    // Compute error bounds for disk intersection
    Vector3f pError(0, 0, 0);

    // Initialize _SurfaceInteraction_ from parametric information
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;
    return true;
}

bool Disk::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute plane intersection for disk

    // Reject disk intersections for rays parallel to the disk's plane
    if (ray.d.z == 0) return false;
    Float tShapeHit = (height - ray.o.z) / ray.d.z;
    if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

    // See if hit point is inside disk radii and $\phimax$
    Point3f pHit = ray(tShapeHit);
    Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$
    Float phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;
    if (phi > phiMax) return false;
    return true;
}

Float Disk::Area() const {
    return phiMax * 0.5 * (radius * radius - innerRadius * innerRadius);
}

Interaction Disk::Sample(const Point2f &u, Float *pdf) const {
    Point2f pd = ConcentricSampleDisk(u);
    Point3f pObj(pd.x * radius, pd.y * radius, height);
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
    if (reverseOrientation) it.n *= -1;
    it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
    *pdf = 1 / Area();
    return it;
}

void Disk::ComputeSphEllipseData(SphEllipseData *data, Point3f o) const{

    Vector3f up = data->oc + data->yd*radius;

    Vector3f down = data->oc - data->yd*radius;

    Float y02 = Dot(Normalize(down), data->yd), y12=Dot(Normalize(up), data->yd);
    Float z02 = Dot(Normalize(down), data->zd), z12=Dot(Normalize(up), data->zd);

    Float zey = (y12+y02)*0.5, zez = (z12+z02)*0.5;

    data->ze = Normalize(zey*data->yd + zez*data->zd);
    data->ye = Cross(data->ze, data->xd);

    //Calculate x1 and x1'
    Float zh = Dot(data->oc, data->zd);
    Float yh = zh * (zey/zez);

    Float cy = Dot(data->oc, data->yd);
    //std::cout << "\ncy" << cy << "\nyh" << yh << "\n";
    Float leg = yh - cy;
    Float x1 = std::sqrt(radius*radius - leg*leg);  

    Vector3f right = Normalize(Vector3f(x1, yh, zh)); 

    data->a = right.x;
    data->b = (y12-y02)*(y12-y02) + (z12-z02)*(z12-z02);
    data->b = 0.5 * std::sqrt(data->b);

    if(data->a > OneMinusEpsilon)
        data->a=1;
    if(data->b > OneMinusEpsilon)
        data->b=1;

    data->alpha = std::asin(data->a);
    data->beta = std::asin(data->b);

    data->at = std::tan(data->alpha);
    data->bt = std::tan(data->beta);

}

Point3f Disk::ReprojectToDisk(const Vector3f &q, const SphEllipseData *data, const Point3f &o) const{
    Vector3f s = Vector3f(Dot(q, data->xd), Dot(q, data->yd), Dot(q, data->zd));
    Float zh = Dot(data->oc, data->zd);
    s = (zh/s.z)*s;
    return o + s.x*data->xd + s.y*data->yd + zh*data->zd;
}

Interaction Disk::Sample(const Interaction &ref, const Point2f &u,
                       Float *pdf) const{
    if(samplingMode == 1){
        Interaction intr = Sample(u, pdf);
        Vector3f wi = intr.p - ref.p;
        if (wi.LengthSquared() == 0)
            *pdf = 0;
        else {
            wi = Normalize(wi);
            // Convert from area measure, as returned by the Sample() call
            // above, to solid angle measure.
            *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
            if (std::isinf(*pdf)) *pdf = 0.f;
        }
        return intr;
    }else{
        //Area preserving parametrization, radial map.
        Interaction it;
        it.p = (*ObjectToWorld)(Point3f(0, 0, height));
        it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
        if (reverseOrientation) it.n *= -1;

        Point3f aux = (*WorldToObject)(ref.p);
        Float rad = aux.x*aux.x + aux.y*aux.y;

        if (rad <= (radius*radius) && aux.z == height){
            *pdf = 1;
            it.p = ref.p;
            return it;
        }
        
        SphEllipseData *data;

        //Initialize constant spherical ellipse data
        data = new SphEllipseData;

        data->oc = (*ObjectToWorld)(Point3f(0, 0, height)) - ref.p;
        data->zd = - Normalize((*ObjectToWorld)(Vector3f(0, 0, 1)));
        if (Dot(data->zd, data->oc)<0) data->zd *= -1;

        data->xd = Normalize(Cross(data->oc, data->zd));
        data->yd = Cross(data->zd, data->xd);  

        ComputeSphEllipseData(data, ref.p);

        if(data->a == 1 || data->b == 1){
            *pdf = 0;      
            return it;
        }

        Float a2 = data->a * data->a;
        Float b2 = data->b * data->b;
        Float oneb2 = 1 - b2;

        double n = (a2 - b2) / (a2 * oneb2);
        double m = (a2 - b2) / oneb2;
        double p = data->b * (1-a2) / (data->a * std::sqrt(oneb2));
        double at_bt = data->at / data->bt;

        double u0;
        if(u[0]<0.25){
            data->signX = 1;
            data->signY = 1;
            u0 = std::min(4*u[0], OneMinusEpsilon);
        }else if(u[0]<0.5){
            data->signX = -1;
            data->signY = 1;
            u0 = std::min(1-4*(u[0]-0.25f), OneMinusEpsilon);
        }else if(u[0]<0.75){
            data->signX = -1;
            data->signY = -1;
            u0 = std::min((u[0]-0.5f)*4, OneMinusEpsilon);
        }else{
            data->signX = 1;
            data->signY = -1;
            u0 = std::min(1-(u[0]-0.75f)*4, OneMinusEpsilon);
        }
   
        data->solidAngle = omega_r2(M_PI*0.5, n, m, p, at_bt);
        Float phi = NR_diskRad(1e-10, n, m, p, at_bt, u0*data->solidAngle);
        

        Float sinphi = std::sin(phi);
        Float cosphi = std::cos(phi);

        Float r = data->a * data->b / std::sqrt(a2*sinphi*sinphi + b2*cosphi*cosphi);
        Float h = std::sqrt(1-r*r);

        h = (1-u[1])*h + u[1];

        r = std::sqrt(1-h*h);

        Vector3f q = data->signX*r*cosphi*data->xd + data->signY*r*sinphi*data->ye + h*data->ze;

        it.p = ReprojectToDisk(q, data, ref.p);

        if ((it.p-ref.p).LengthSquared() == 0)
            *pdf = 0;
        else {
            *pdf = 1/(4*data->solidAngle);
            if (std::isinf(*pdf)) *pdf = 0;
        }
        return it;

    }
}

std::shared_ptr<Disk> CreateDiskShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params) {
    Float height = params.FindOneFloat("height", 0.);
    Float radius = params.FindOneFloat("radius", 1);
    Float inner_radius = params.FindOneFloat("innerradius", 0);
    Float phimax = params.FindOneFloat("phimax", 360);
    int samplingMode = params.FindOneInt("samplingMode", 1);
    return std::make_shared<Disk>(o2w, w2o, reverseOrientation, height, radius,
                                  inner_radius, phimax, samplingMode);
}

}  // namespace pbrt
