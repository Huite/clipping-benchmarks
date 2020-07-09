import numpy as np
cimport numpy as cnp


cdef extern from "clipping.h":
    void area_of_intersection(long, long, long, double*, double*, double*)


cdef _cpp_area_of_intersection(int ntriangles, int nvertex, int ndim, a, b):
    cdef cnp.ndarray[double, ndim=3, mode="c"] aa
    cdef cnp.ndarray[double, ndim=3, mode="c"] bb 
    cdef cnp.ndarray[double, ndim=1, mode="c"] out
    a = np.ascontiguousarray(a, dtype=np.float64)
    b = np.ascontiguousarray(b, dtype=np.float64)
    aa = a
    bb = b
    out = np.zeros((ntriangles,), dtype=np.float64)
    area_of_intersection(ntriangles, nvertex, ndim, &aa[0, 0, 0], &bb[0, 0, 0], &out[0])
    return out


def cpp_area_of_intersection(a, b):
    ntriangles, nvertex, ndim = a.shape
    return _cpp_area_of_intersection(ntriangles, nvertex, ndim, a, b)
