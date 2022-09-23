use crate::common::{Point, Vector};
use numpy::ndarray::{Array, ArrayView2, ArrayViewMut2};

fn push(mut array: ArrayViewMut2<f64>, n: &mut usize, value: Point) {
    array[[*n, 0]] = value.x;
    array[[*n, 1]] = value.y;
    *n += 1;
}

fn copy_to(src: ArrayView2<f64>, mut dst: ArrayViewMut2<f64>, n: usize) {
    for i in 0..n {
        dst[[i, 0]] = src[[i, 0]];
        dst[[i, 1]] = src[[i, 1]];
    }
}

fn inside(p: Point, r: Point, U: Vector) -> bool {
    U.x * (p.y - r.y) > U.y * (p.x - r.x)
}

fn intersection(a: Point, V: Vector, r: Point, N: Vector) -> Option<Point> {
    let W = Vector::new(r.x - a.x, r.y - a.y);
    let nw = N.dot(W);
    let nv = N.dot(V);
    if nv != 0.0 {
        let t = nw / nv;
        Some(Point::new(a.x + t * V.x, a.y + t * V.y))
    } else {
        None
    }
}

fn polygon_area(polygon: ArrayView2<f64>, length: usize) -> f64 {
    let mut area = 0.0;
    let a = Point::new(polygon[[0, 0]], polygon[[0, 1]]);
    let b = Point::new(polygon[[1, 0]], polygon[[1, 1]]);
    let mut U = Vector::new(b.x - a.x, b.y - a.y);
    for i in 2..length {
        let c = Point::new(polygon[[i, 0]], polygon[[i, 1]]);
        let V = Vector::new(a.x - c.x, a.y - c.y);
        area += U.cross(V).abs();
        U = V;
    }

    0.5 * area
}

pub fn clip_polygons(polygon: ArrayView2<'_, f64>, clipper: ArrayView2<'_, f64>) -> f64 {
    let n_max = polygon.shape()[0] * 2;
    let mut output = Array::zeros((n_max, 2));
    let mut subject = Array::zeros((n_max, 2));

    let mut n_output = polygon.shape()[0];
    copy_to(polygon, output.view_mut(), n_output);

    let n_clip = clipper.shape()[0];
    let mut r = Point::new(clipper[[n_clip - 1, 0]], clipper[[n_clip - 1, 1]]);
    for i in 0..n_clip {
        let s = Point::new(clipper[[i, 0]], clipper[[i, 1]]);
        let U = Vector::new(s.x - r.x, s.y - r.y);
        let N = Vector::new(-U.y, U.x);

        if U.x == 0. && U.y == 0. {
            continue;
        }

        let length = n_output;
        copy_to(output.view(), subject.view_mut(), length);
        n_output = 0;
        let mut a = Point::new(subject[[length - 1, 0]], subject[[length - 1, 1]]);
        let mut a_inside = inside(a, r, U);

        for j in 0..length {
            let b = Point::new(subject[[j, 0]], subject[[j, 1]]);
            let V = Vector::new(b.x - a.x, b.y - a.y);

            if V.x == 0. && V.y == 0. {
                continue;
            }

            let mut b_inside = inside(b, r, U);

            if b_inside {
                if !a_inside {
                    if let Some(p) = intersection(a, V, r, N) {
                        push(output.view_mut(), &mut n_output, p);
                    }
                }
                push(output.view_mut(), &mut n_output, b);
            } else if a_inside {
                if let Some(p) = intersection(a, V, r, N) {
                    push(output.view_mut(), &mut n_output, p);
                } else {
                    b_inside = true;
                    push(output.view_mut(), &mut n_output, b);
                }
            }
            a = b;
            a_inside = b_inside;
        }
        if n_output < 3 {
            return 0.0;
        }
        r = s;
    }

    return polygon_area(output.view(), n_output);
}
