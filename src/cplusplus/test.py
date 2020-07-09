import numpy as np
from cpp_clipping import area_of_intersection


a = np.random.rand(100, 3, 2)
b = np.random.rand(100, 3, 2)
c = area_of_intersection(a, b)
print(c)