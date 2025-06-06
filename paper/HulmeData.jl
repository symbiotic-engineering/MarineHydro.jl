"""
Hulme A. The wave forces acting on a floating hemisphere undergoing forced periodic oscillations. Journal of Fluid Mechanics. 1982;121:443-463. doi:10.1017/S0022112082001980
"""
#  A and B of surging spehre table2
K_surge = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.4, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 10.0]

A_surge = [0.5223, 0.5515, 0.5848, 0.6175, 0.6439, 0.6586, 0.6682, 0.6421, 0.6127, 0.5740, 0.4860, 0.4038, 0.3371, 0.2866, 0.2493, 0.1961, 0.1720, 0.1634, 0.1620, 0.1641, 0.1679, 0.1772, 0.1865, 0.1949, 0.2022, 0.2085, 0.2732]

B_surge = [0.001, 0.0082, 0.0255, 0.0557, 0.0987, 0.1516, 0.2092, 0.2653, 0.3145, 0.3535, 0.3978, 0.406, 0.3929, 0.3695, 0.3424, 0.2769, 0.2237, 0.1826, 0.151, 0.1266, 0.1073, 0.0794, 0.0608, 0.0479, 0.0386, 0.0317, 0.0]

# Added mass and dampingg coeff of heaving hemisphere : Table1

K_heave = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

A_heave = [0.8764, 0.8627, 0.7938, 0.7157, 0.6452, 0.5861, 0.5381, 0.4999, 0.4698, 0.4464, 0.4284, 0.4047, 0.3924, 0.3871, 0.3864, 0.3884, 0.3988, 0.4111, 0.4322, 0.4471, 0.4574, 0.4647, 0.4700, 0.4740, 0.4771]

B_heave = [0.1036, 0.1816, 0.2793, 0.3254, 0.3410, 0.3391, 0.3271, 0.3098, 0.2899, 0.2691, 0.2484, 0.2096, 0.1756, 0.1469, 0.1229, 0.1031, 0.0674, 0.0452, 0.0219, 0.0116, 0.0066, 0.0040, 0.0026, 0.0017, 0.0012]


### diffraction forces analyytical by https://www.sciencedirect.com/science/article/pii/S002980182202813X
