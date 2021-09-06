// shapes/rectangle.cpp*
#include "project/rectangle.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

Rectangle::Rectangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
                     bool reverseOrientation, Float width, Float height, int samplingMode)
    : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        width(width), height(height), samplingMode(samplingMode){}

// Rectangle Method Definitions
Bounds3f Rectangle::ObjectBound() const {
    return Bounds3f(Point3f(-width*0.5, -height*0.5, 0),
                    Point3f(width*0.5, height*0.5, 0));
}

bool Rectangle::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                     bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);
    
    Bounds3f rectangle= Bounds3f(Point3f(-width*0.5, -height*0.5, 0),
                    							Point3f(width*0.5, height*0.5, 0));

    // Reject rectangle intersections for rays parallel to the rectangle's plane
    if (ray.d.z == 0) return false;
    Float tShapeHit, t2;
    if(!rectangle.IntersectP(ray,&tShapeHit,&t2)){
        return false;
    }

    // See if hit point is inside rectangle radii and $\phimax$
    Point3f pHit = ray(tShapeHit);

    // Find parametric representation of rectangle hit
    Float u=pHit.x/width + 0.5;
    Float v=pHit.y/height + 0.5;
    Vector3f dpdu(width, 0., 0.);
    Vector3f dpdv(0., height, 0.);
    Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

    // Refine rectangle intersection point
    pHit.z = 0;

    // Compute error bounds for rectangle intersection
    Vector3f pError(0, 0, 0);

    // Initialize _SurfaceInteraction_ from parametric information
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;
    return true;
}

bool Rectangle::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);
    
    Bounds3f rectangle= Bounds3f(Point3f(-width*0.5, -height*0.5, 0),
                    			Point3f(width*0.5, height*0.5, 0));

    // Reject rectangle intersections for rays parallel to the rectangle's plane
    if (ray.d.z == 0) return false;
    Float tShapeHit, t2;
		return rectangle.IntersectP(ray,&tShapeHit,&t2);
}

Float Rectangle::Area() const {
    return width*height;
}

Interaction Rectangle::Sample(const Point2f &u, Float *pdf) const {
    Point3f pObj(width*(u[0]-0.5), height*(u[1]-0.5), 0);   
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
    if (reverseOrientation) it.n *= -1;
    it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
    *pdf = 1 / Area();
    return it;
}

//Función que calcula los parámetros usados posteriormente
void ComputeSphRectangleData(SphRectangleData *data, const Point3f &o){
    
    data->o= o;
    Vector3f d = data->s - data->o;
    data->z0 = Dot(d, data->z);
    if(data->z0 > 0){
        data->z *= -1;
        data->z0 *= -1;
    }

    data->z0sq = data->z0 * data->z0;
    data->x0 = Dot(d, data->x);
    data->y0 = Dot(d, data->y);
    data->x1 = data->x0 + data->ex1;
    data->y1 = data->y0 + data->ey1;
    data->y0sq = data->y0 * data->y0;
    data->y1sq = data->y1 * data->y1;

    Vector3f v00 = Vector3f(data->x0, data->y0, data->z0);
    Vector3f v01 = Vector3f(data->x0, data->y1, data->z0);
    Vector3f v10 = Vector3f(data->x1, data->y0, data->z0);
    Vector3f v11 = Vector3f(data->x1, data->y1, data->z0);
    //Calcular las normales a los planos laterales de la pirámide
    Vector3f n0 = Normalize(Cross(v00, v10));
    Vector3f n1 = Normalize(Cross(v10, v11));
    Vector3f n2 = Normalize(Cross(v11, v01));
    Vector3f n3 = Normalize(Cross(v01, v00));
    //Calcular los ángulos internos de la pirámide (gamma_i)
    Float g0 = std::acos(-Dot(n0,n1));
    Float g1 = std::acos(-Dot(n1,n2));
    Float g2 = std::acos(-Dot(n2,n3));
    Float g3 = std::acos(-Dot(n3,n0));

    data->b0 = n0.z;
    data->b1 = n2.z;
    data->b0sq = data->b0 * data->b0;
    data->k = 2*Pi - g2 - g3;
    data->solidAngle = g0 + g1 -data->k;
}

Interaction Rectangle::Sample(const Interaction &ref, const Point2f &u, Float *pdf) const{
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
        //Area-preserving parametrization

        SphRectangleData *sphRectData;

        //1. Initializamos el sistema de referencia y los parámetros usados 
        sphRectData = new SphRectangleData;
        sphRectData->x=(*ObjectToWorld)(Vector3f(1, 0, 0));
        sphRectData->y=(*ObjectToWorld)(Vector3f(0, 1, 0));
        sphRectData->ex1=width;
        sphRectData->ey1=height;
        sphRectData->z=Cross(sphRectData->x, sphRectData->y);
        sphRectData->s=(*ObjectToWorld)(Point3f(-width*0.5, -height*0.5, 0));

        ComputeSphRectangleData(sphRectData, ref.p);
        
        Float au = u[0] * sphRectData->solidAngle + sphRectData->k;
        Float fu = (std::cos(au) * sphRectData->b0 - sphRectData->b1) / std::sin(au);
        Float cu = 1/std::sqrt(fu*fu + sphRectData->b0sq) * (fu>0 ? +1 : -1);
        cu = Clamp(cu, -1, 1); // avoid NaNs
        // 2. Calculamos ’xu’
        Float xu = -(cu * sphRectData->z0) / std::sqrt(1 - cu*cu);
        xu = Clamp(xu, sphRectData->x0, sphRectData->x1); // avoid Infs
        // 3. Calculamos ’yv’
        Float d = std::sqrt(xu*xu + sphRectData->z0sq);
        Float h0 = sphRectData->y0 / std::sqrt(d*d + sphRectData->y0sq);
        Float h1 = sphRectData->y1 / std::sqrt(d*d + sphRectData->y1sq);
        Float hv = h0 + u[1] * (h1-h0), hv2 = hv*hv;
        Float yv = (hv2 < OneMinusEpsilon) ? (hv*d)/std::sqrt(1-hv2) : sphRectData->y1;
        // 4. Transformamos (xu,yv,z0) a coordenadas del mundo
        Interaction it;
        it.p = sphRectData->o + sphRectData->x*xu + yv*sphRectData->y + sphRectData->z0*sphRectData->z;

        it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
        if (reverseOrientation) it.n *= -1;

        //Si el punto pertenece a la fuente de luz no tenemos en cuenta la muestra
        if ((it.p-ref.p).LengthSquared() == 0)
            *pdf = 0;
        else {
            *pdf = 1/sphRectData->solidAngle;
            if (std::isinf(*pdf)) *pdf = 0.f;
        }

        return it;
    }
}

Float Rectangle::Pdf(const Interaction &ref, const Vector3f &wi) const {
    if(samplingMode==1){
        // Intersect sample ray with area light geometry
        Ray ray = ref.SpawnRay(wi);
        Float tHit;
        SurfaceInteraction isectLight;
        // Ignore any alpha textures used for trimming the shape when performing
        // this intersection. Hack for the "San Miguel" scene, where this is used
        // to make an invisible area light.
        if (!Intersect(ray, &tHit, &isectLight, false)) return 0;

        // Convert light sample weight to solid angle measure
        Float pdf = DistanceSquared(ref.p, isectLight.p) /
                    (AbsDot(isectLight.n, -wi) * Area());
        if (std::isinf(pdf)) pdf = 0.f;
        return pdf;
    }else{
        Ray ray = ref.SpawnRay(wi);
        Float tHit;
        SurfaceInteraction isectLight;
        if (!Intersect(ray, &tHit, &isectLight, false)) return 0;

        //Calculamos la función de densidad respecto al ángulo sólido
        // en la dirección wi
        SphRectangleData *sphRectData;

        //1. Initializamos el sistema de referencia y los parámetros usados 
        sphRectData = new SphRectangleData;
        sphRectData->x=(*ObjectToWorld)(Vector3f(1, 0, 0));
        sphRectData->y=(*ObjectToWorld)(Vector3f(0, 1, 0));
        sphRectData->ex1=width;
        sphRectData->ey1=height;
        sphRectData->z=Cross(sphRectData->x, sphRectData->y);
        sphRectData->s=(*ObjectToWorld)(Point3f(-width*0.5, -height*0.5, 0));

        ComputeSphRectangleData(sphRectData, ref.p);

        Float pdf = 1 / sphRectData->solidAngle;
        if (std::isinf(pdf)) pdf = 0.f;
        return pdf;
    }
}

std::shared_ptr<Rectangle> CreateRectangleShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params) {
	Float width = params.FindOneFloat("width", 1);
    Float height = params.FindOneFloat("height", 1);
    int samplingMode = params.FindOneInt("samplingMode", 1);
    return std::make_shared<Rectangle>(o2w, w2o, reverseOrientation, width, height, samplingMode);
}

}  // namespace pbrt
