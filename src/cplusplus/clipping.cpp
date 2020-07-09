# include "clipping.h"


struct Point { double x, y; };
struct DVector { double x, y; };
struct IntersectionResult { bool succes; Point p; };


double cross_product(DVector U, DVector V) {
    return (U.x * V.y - U.y * V.x);
}


double dot_product(DVector U, DVector V) {
    return (U.x * V.x + U.y * V.y);
}


bool inside(Point p, Point r, DVector U) {
    return (U.x * (p.y - r.y) > U.y* (p.x - r.x));
}


IntersectionResult intersection(Point a, DVector V, Point r, DVector N) {
    DVector W = { r.x - a.x, r.y - a.y };
    double nw = dot_product(N, W);
    double nv = dot_product(N, V);
    if (nv != 0.0) {
        double t = nw / nv;
        Point p = { a.x + t * V.x, a.y + t * V.y };
        return { true, p };
    }
    else {
        Point p = { 0.0, 0.0 };
        return { false, p };
    }
}


double polygon_area(std::vector<Point> polygon, int length) {
    double area = 0.0;
    Point a = polygon[0];
    Point b = polygon[1];
    DVector U = { b.x - a.x, b.y - a.y };
    for (int i = 2; i < length; i++) {
        Point c = polygon[i];
        DVector V = { a.x - c.x, a.y - c.y };
        area += abs(cross_product(U, V));
        b = c;
        U = V;
    }
    return 0.5 * area;
}


void copyto(std::vector<Point> src, std::vector<Point> dst, int n) {
    for (int i = 0; i < n; i++) {
        dst[i] = src[i];
    };
};


double clip_polygons(std::vector<Point> polygon, std::vector<Point> clipper, int n_max) {
    std::vector<Point> subject(n_max), output(n_max);
    Point r, s, a, b;
    DVector U, V, N;
    bool a_inside, b_inside;
    int n_poly, n_clip, n_output, length;
    IntersectionResult int_res;

    n_poly = polygon.size();
    n_clip = clipper.size();
    n_output = n_poly;
    for (int ii = 0; ii < n_output; ii++) {
        output[ii] = polygon[ii];
    };

    r = clipper[n_clip - 1];
    for (int i = 0; i < n_clip; i++) {
        s = clipper[i];
        U = { s.x - r.x, s.y - r.y };
        N = { -U.y, U.x };

        if (U.x == 0.0 && U.y == 0.0) {
            continue;
        }

        length = n_output;
        for (int ii = 0; ii < length; ii++) {
            subject[ii] = output[ii];
        };
        n_output = 0;

        a = subject[length - 1];
        a_inside = inside(a, r, U);
        for (int j = 0; j < length; j++) {
            b = subject[j];
            V = { b.x - a.x, b.y - a.y };

            if (V.x == 0.0 && V.y == 0.0) {
                continue;
            }

            b_inside = inside(b, r, U);
            if (b_inside) {
                if (!a_inside) {
                    int_res = intersection(a, V, r, N);
                    if (int_res.succes) { 
                        output[n_output++] = int_res.p;
                    };
                }
                output[n_output++] = b;
            }
            else if (a_inside) {
                int_res = intersection(a, V, r, N);
                if (int_res.succes) {
                    output[n_output++] = int_res.p;
                } else {
                    b_inside = true; 
                    output[n_output++] = b; 
                }
            }
            a = b;
            a_inside = b_inside;
        }
        if (n_output < 3) return 0.0;
        r = s;
    }
    double area = polygon_area(output, n_output);
    return area;
}


void cpp_area_of_intersection(
    long ntriangles,
    long nvertex,
    long ndim,
    double* p,//olygons
    double* c,//lippers
    double* areas
) {
    long size = ndim * nvertex * ntriangles;
    std::vector<Point> polygon(nvertex), clipper(nvertex);
    for (auto i = 0; i < ntriangles; i++) {
        long start = i * 6;
        for (auto j = 0; j < nvertex; j++) {
            long first = start + j * 2;
            long second = start + j * 2 + 1;
            polygon[j] = { p[first], p[second] };
            clipper[j] = { c[first], c[second] };
        }
        areas[i] = clip_polygons(polygon, clipper, 2 * nvertex);
    }
}
