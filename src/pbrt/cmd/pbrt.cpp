// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/pbrt.h>
#include <iomanip>
#include <pbrt/cpu/render.h>

#ifdef PBRT_BUILD_GPU_RENDERER
#include <pbrt/gpu/memory.h>
#endif  // PBRT_BUILD_GPU_RENDERER

#include <pbrt/options.h>
#include <pbrt/parser.h>
#include <pbrt/scene.h>
#include <pbrt/util/args.h>
#include <pbrt/util/check.h>
#include <pbrt/util/error.h>
#include <pbrt/util/log.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/print.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/string.h>
#include <pbrt/wavefront/wavefront.h>
#include <pbrt/pbrt.h>

#include <string>
#include <vector>

#include <iostream>
#include <array>
#include "pbrt/util/splines.h"

#include "pbrt/shapes.h"


using namespace pbrt;
using namespace std;


// main program
void test_bezier_shape() {
    cout << endl << "================================================================================" << endl;
    cout << "    Bezier Shape " << endl;
    std::array<Point3f, 4> cp;
    Bounds3f b = BoundCubicBezier(pstd::MakeConstSpan(cp), 0.f, 1.f);
    b = Expand(b, 1e-3 * Length(b.Diagonal()));
    for (Float u = 0; u <= 1.f; u += 1.f / 1024.f) {
        Point3f p = EvaluateCubicBezier(pstd::MakeConstSpan(cp), u);
        bool inside = Inside(p, b);
        if (inside) {
//            cout << "inside   " << p << " @ u = " << u << " not in " << b << endl;
        } else {
//            cout << "outside   " << p << " @ u = " << u << " not in " << b << endl;
        }
    }

//
//    RNG rng;
//    for (int i = 0; i < 1000; ++i) {
//        for (int j = 0; j < 4; ++j)
//            for (int c = 0; c < 3; ++c)
//                cp[j][c] = -5.f + 10.f * rng.Uniform<Float>();
//
//        Bounds3f b = BoundCubicBezier(pstd::MakeConstSpan(cp), 0.f, 1.f);
//        b = Expand(b, 1e-3 * Length(b.Diagonal()));
//        for (Float u = 0; u <= 1.f; u += 1.f / 1024.f) {
//            Point3f p = EvaluateCubicBezier(pstd::MakeConstSpan(cp), u);
//            cout << (Inside(p, b)) << p << " @ u = " << u << " not in " << b << endl;
//        }
//    }
//
    cout << endl << "================================================================================" << endl;

}

void test_clamp() {
    cout << endl << "================================================================================" << endl;
    cout << "clamp     test" << endl << endl;

    Float low = 0.0;
    Float high = 2.0;

    Float values[10] = {-0.000001, 0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0, 3.0};

    for (int i = 0; i < 10; i++) {
        Float res = Clamp(values[i], low, high);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "low " << low << "   high:  " << high
             << "   value " << values[i] << "  clamped to  ==>      " << res << endl;
    }
    cout << endl << "================================================================================" << endl;

}

void test_lerp_f64() {
    cout << endl << "================================================================================" << endl;
    cout << "lerp f64    test" << endl << endl;

    Float min = 0.0;
    Float max = 2.0;

    Float arr[8] = {0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0};

    for (int i = 0; i < 8; i++) {
        Float res = Lerp(arr[i], min, max);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "min " << min << "   max:  " << max
             << "  res     " << res << endl;
    }
    cout << endl << "================================================================================" << endl;
}


void test_lerp_vec() {
    cout << endl << "================================================================================" << endl;
    cout << "lerp point3d    test" << endl << endl;

    Point3f a = Point3f(1.0f, 2.0f, 3.0f);
    Point3f b = Point3f(14.0f, -15.0f, 16.0f);

    Float arr[8] = {0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0};

    for (int i = 0; i < 8; i++) {
        Point3f res = Lerp(arr[i], a, b);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "a " << a << "   b:  " << b
             << "  res     " << res << endl;
    }
    cout << endl << "================================================================================" << endl;

}


void test_blossom_bezier() {
    cout << endl << "================================================================================" << endl;
    cout << "bloosom bezier test" << endl << endl;

    pstd::span<const Point3f> cp = {Point3f(0.0f, 1.0f, 1.0f),
                                    Point3f(1.0f, 1.0f, 2.0f),
                                    Point3f(2.0f, 2.0f, 2.0f),
                                    Point3f(3.0f, 3.0f, 3.0f)
    };

    Float u0 = 0.0;
    Float u1 = 0.25;
    Float u2 = 1.0;

    Point3f f1 = BlossomCubicBezier(cp, u0, u0, u0);
    Point3f f2 = BlossomCubicBezier(cp, u0, u0, u1);
    Point3f f3 = BlossomCubicBezier(cp, u1, u1, u2);
    Point3f f4 = BlossomCubicBezier(cp, u2, u2, u2);

    int i = 1;
    for (Point3f p: cp) {
        cout << "point  " << i << "   x " << p.x << "  y   " << p.y << "   p.z  " << p.z << endl;
        i++;
    }
    cout << "u0  " << u0 << endl;
    cout << "u1  " << u1 << endl;
    cout << "u2  " << u2 << endl;

    cout << "u0 " << u0 << "   u0  " << u0 << "  u0   " << u0 << "   f1       " << f1 << endl;
    cout << "u0 " << u0 << "   u0  " << u0 << "  u1   " << u1 << "   f2       " << f2 << endl;
    cout << "u1 " << u1 << "   u1  " << u1 << "  u2   " << u2 << "   f3       " << f3 << endl;
    cout << "u2 " << u2 << "   u2  " << u2 << "  u2   " << u2 << "   f4       " << f4 << endl;
    cout << endl << "================================================================================" << endl;

}


void test_coordinate_system() {
    cout << endl << "================================================================================" << endl;
    cout << "coordinate_system   " << endl << endl;

    Vector3f v1[4] = {Vector3f(1.0, 0.0, 0.0), Vector3f(1.0, 1.0, 0.0), Vector3f(1.0, 0.0, 1.0),
                      Vector3f(1.0, 1.0, 1.0)};
    Vector3f v2[4] = {Vector3f(0.0, 0.0, 0.0), Vector3f(0.0, 0.0, 0.0), Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0)};
    Vector3f v3[4] = {Vector3f(0.0, 0.0, 0.0), Vector3f(0.0, 0.0, 0.0), Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0)};

    for (int i = 0; i < 4; i++) {
        CoordinateSystem(v1[i], &v2[i], &v3[i]);
        cout << "v1 " << v1[i].x << "     " << v1[i].x << "    " << v1[i].z << endl;
        cout << "v2 " << v2[i].x << "     " << v2[i].x << "    " << v2[i].z << endl;
        cout << "v3 " << v3[i].x << "     " << v3[i].x << "    " << v3[i].z << endl;


        cout << endl;
    }


    cout << endl << "================================================================================" << endl;

}


void test_look_at() {
    cout << endl << "================================================================================" << endl;
    cout << "   LookAt   " << endl << endl;


    Point3f pos = Point3f(1.0, 2.0, 3.0);
    Point3f look = Point3f(5.0, 6.0, 7.0);
    Vector3f up[3] = {Vector3f(1.0, 0.0, 0.0), Vector3f(0.0, 1.0, 0.0), Vector3f(0.0, 0.0, 1.0)};

    cout << "pos:   " << pos.x << "   " << pos.y << "    " << pos.z << endl;
    cout << "look:   " << look.x << "   " << look.y << "    " << look.z << endl;

    for (int i = 0; i < 3; i++) {
        Transform t = LookAt(pos, look, up[i]);

        SquareMatrix<4> m = t.GetMatrix();
        SquareMatrix<4> m_inv = t.GetInverseMatrix();

        cout << "up:   " << up[i].x << "   " << up[i].y << "    " << up[i].z << endl;

        cout << "m[..][0]:   " << m[0][0] << "   " << m[1][0] << "    " << m[2][0] << "  " << m[3][0] << endl;
        cout << "m[..][1]:   " << m[1][1] << "   " << m[1][1] << "    " << m[2][1] << "  " << m[3][1] << endl;
        cout << "m[..][2]:   " << m[1][2] << "   " << m[1][2] << "    " << m[2][2] << "  " << m[3][2] << endl;
        cout << "m[..][3]:   " << m[1][3] << "   " << m[1][3] << "    " << m[2][3] << "  " << m[3][3] << endl;


        cout << "m_inv[..][0]:   " << m_inv[0][0] << "   " << m_inv[1][0] << "    " << m_inv[2][0] << "  "
             << m_inv[3][0] << endl;
        cout << "m_inv[..][1]:   " << m_inv[1][1] << "   " << m_inv[1][1] << "    " << m_inv[2][1] << "  "
             << m_inv[3][1] << endl;
        cout << "m_inv[..][2]:   " << m_inv[1][2] << "   " << m_inv[1][2] << "    " << m_inv[2][2] << "  "
             << m_inv[3][2] << endl;
        cout << "m_inv[..][3]:   " << m_inv[1][3] << "   " << m_inv[1][3] << "    " << m_inv[2][3] << "  "
             << m_inv[3][3] << endl;
        cout << endl;

    }


    cout << endl << "================================================================================" << endl;

}


void test_union() {
    cout << endl << "================================================================================" << endl;
    cout << "   Union   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Point3f p2 = Point3f(0.5, -0.5, -1.5);
    Point3f p3 = Point3f(1.5, 13.0, 3.0);

    Bounds3f bb1 = Bounds3f(p0, p1);
    Bounds3f bb2 = Bounds3f(p2, p3);

    Bounds3f u = Union(bb1, bb2);


    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl;

    cout << "bb2.min:   " << bb2.pMin.x << "   " << bb2.pMin.y << "    " << bb2.pMin.z << endl;
    cout << "bb2.max:   " << bb2.pMax.x << "   " << bb2.pMax.y << "    " << bb2.pMax.z << endl;

    cout << "union.min:   " << u.pMin.x << "   " << u.pMin.y << "    " << u.pMin.z << endl;
    cout << "union.max:   " << u.pMax.x << "   " << u.pMax.y << "    " << u.pMax.z << endl;


    cout << endl << "================================================================================" << endl;

}


void test_expand() {
    cout << endl << "================================================================================" << endl;
    cout << "   Expand   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Bounds3f bb1 = Bounds3f(p0, p1);

    Bounds3f e = Expand(bb1, 0.75);

    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl;

    cout << "expand.min:   " << e.pMin.x << "   " << e.pMin.y << "    " << e.pMin.z << endl;
    cout << "expand.max:   " << e.pMax.x << "   " << e.pMax.y << "    " << e.pMax.z << endl;

    cout << endl << "================================================================================" << endl;
}


void test_ray_bounds() {
    cout << endl << "================================================================================" << endl;
    cout << "   rayBounds   " << endl << endl;

    Ray r = Ray(Point3f(1.0, 2.0, 3.0), Vector3f(-2.0, -2.0, -1.5));
    Float tMax = 1.0;
    Bounds3f rayBounds(Point3f(0, 0, 0), Point3f(0, 0, Length(r.d) * tMax));


    cout << "ray.p:   " << r.o.x << "   " << r.o.y << "    " << r.o.z << endl;
    cout << "ray.d:   " << r.d.x << "   " << r.d.y << "    " << r.d.z << "    length(r.d)  " << Length(r.d) << endl;


    cout << "rayBounds.min:   " << rayBounds.pMin.x << "   " << rayBounds.pMin.y << "    " << rayBounds.pMin.z << endl;
    cout << "rayBounds.max:   " << rayBounds.pMax.x << "   " << rayBounds.pMax.y << "    " << rayBounds.pMax.z << endl;

    cout << endl << "================================================================================" << endl;
}


void test_overlaps() {
    cout << endl << "================================================================================" << endl;
    cout << "   Overlaps   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Point3f p2 = Point3f(0.5, 0.5, 0.5);
    Point3f p3 = Point3f(1.5, 13.0, 3.0);

    Point3f p4 = Point3f(10.5, 10.5, 10.5);
    Point3f p5 = Point3f(1.5, 13.0, 13.0);

    Bounds3f bb1 = Bounds3f(p0, p1);
    Bounds3f bb2 = Bounds3f(p2, p3);
    Bounds3f bb3 = Bounds3f(p4, p5);

    bool overlaps = Overlaps(bb1, bb2);
    bool overlaps2 = Overlaps(bb1, bb3);

    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl << endl;

    cout << "bb2.min:   " << bb2.pMin.x << "   " << bb2.pMin.y << "    " << bb2.pMin.z << endl;
    cout << "bb2.max:   " << bb2.pMax.x << "   " << bb2.pMax.y << "    " << bb2.pMax.z << endl << endl;

    cout << "bb3.min:   " << bb3.pMin.x << "   " << bb3.pMin.y << "    " << bb3.pMin.z << endl;
    cout << "bb3.max:   " << bb3.pMax.x << "   " << bb3.pMax.y << "    " << bb3.pMax.z << endl << endl;

    cout << "overlaps:   " << std::boolalpha << overlaps << endl;
    cout << "overlaps2:   " << std::boolalpha << overlaps2 << endl;

    cout << endl << "================================================================================" << endl;
}


void test_evaluate_cubic_bezier() {
    cout << endl << "================================================================================" << endl;
    cout << "evaluate_cubic_bezier " << endl << endl;

    Vector3f deriv[3] = {Vector3f(1.0, 0.0, 0.0), Vector3f(1.0, 0.0, 1.0), Vector3f(1.0, 1.0, 0.0)};

    pstd::span<const Point3f> cp = {Point3f(0.0f, 1.0f, 1.0f),
                                    Point3f(1.0f, 1.0f, 2.0f),
                                    Point3f(2.0f, 2.0f, 2.0f),
                                    Point3f(3.0f, 3.0f, 3.0f)
    };

    int i = 1;
    for (Point3f p: cp) {
        cout << "point  " << i << "   x " << p.x << "  y   " << p.y << "   p.z  " << p.z << endl;
        i++;
    }

    Float u = 0.5;

    for (int j = 0; j < 3; j++) {
        Point3f eval = EvaluateCubicBezier(cp, u, &deriv[j]);
        cout << "evaluated point  " << j << "   x " << eval.x << "  y   " << eval.y << "   p.z  " << eval.z << endl;
    }

    Point3f eval = EvaluateCubicBezier(cp, u, nullptr);
    cout << "evaluated point dreiv = nullptr   " << "   x " << eval.x << "  y   " << eval.y << "   p.z  " << eval.z
         << endl;

    cout << endl << "================================================================================" << endl;
}


void test_intersect_curve() {
    cout << endl << "================================================================================" << endl;
    cout << "intersect curve" << endl << endl;

// Shape::Create()
////     CurveCommon cc =   CurveCommon;
//   auto curve =pbrt::Curve::Create();
//    const Transform renderFromObject = Transform();
//    const Transform objectFromRender = Transform();
//    bool reverseOrientation = false;
//    // const ParameterDictionary parameters = ParameterDictionary();
//    pstd::span<const Point3f> cp = {Point3f(0.0f, 1.0f, 1.0f),
//                                    Point3f(1.0f, 1.0f, 2.0f),
//                                    Point3f(2.0f, 2.0f, 2.0f),
//                                    Point3f(3.0f, 3.0f, 3.0f)
//    };
//    ParsedParameterVector parameterVector;
//    FileLoc f;
//    ParsedParameter *param = new ParsedParameter( f);
//param.
//    const ParameterDictionary parameters = ParameterDictionary({"P", cp}, RGBColorSpace::sRGB);
//
//    const Allocator alloc = Allocator{};
//
//   auto shapes = Curve::Create(&renderFromObject, &objectFromRender, reverseOrientation,
//                                parameters, nullptr, alloc);
//




//
//    Vector3f deriv[3] = {Vector3f(1.0, 0.0, 0.0), Vector3f(1.0, 0.0, 1.0), Vector3f(1.0, 1.0, 0.0)};
//
//
//    int i = 1;
//    for (Point3f p: cp) {
//        cout << "point  " << i << "   x " << p.x << "  y   " << p.y << "   p.z  " << p.z << endl;
//        i++;
//    }
//
//    Float u = 0.5;
//
//    for (int j = 0; j < 3; j++) {
//        Point3f eval = EvaluateCubicBezier(cp, u, &deriv[j]);
//        cout << "evaluated point  " << j << "   x " << eval.x << "  y   " << eval.y << "   p.z  " << eval.z << endl;
//    }
//
//    Point3f eval = EvaluateCubicBezier(cp, u, nullptr);
//    cout << "evaluated point dreiv = nullptr   " << "   x " << eval.x << "  y   " << eval.y << "   p.z  " << eval.z << endl;
//

    cout << endl << "================================================================================" << endl;

}

// main program
int main(int argc, char *argv[]) {
    test_lerp_f64();
    test_lerp_vec();
    test_blossom_bezier();
    test_coordinate_system();
    test_look_at();
    test_union();
    test_bezier_shape();
    test_expand();
    test_ray_bounds();
    test_overlaps();
    test_clamp();
    test_evaluate_cubic_bezier();
    test_intersect_curve();

    return 0;
}

