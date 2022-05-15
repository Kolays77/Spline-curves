from geomdl import NURBS

from geomdl import utilities
from geomdl import exchange
from geomdl import operations
from geomdl import multi
from geomdl.visualization import VisMPL

"""
    Examples for the NURBS-Python Package
    Released under MIT License
    Developed by Onur Rauf Bingol (c) 2018
"""


# Fix file path
#os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Create a B-Spline curve instance
crv = NURBS.Curve()

# Set degree
crv.degree = 3

# Set unweights control points with weights
crv.ctrlpts = [[-70.0, -76.0, 0.0], [-70.0, 75.0, 0.0], [74, 75, 0.0], [74, -77, 0.0], [-40, -76, 0.0]]
crv.weights = [1.0, 0.5, 4.0, 5.0, 1.0]
# Set knot vector
crv.knotvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]

# Split the curve
t_cup = 0.7

curve_list = operations.split_curve(crv, t_cup)
curves = multi.CurveContainer(curve_list)

# Move the 1st curve a little bit far away from the 2nd curve
#c2tan = curves[1].tangent(0.0, normalize=True)
#c2tanvec = [-3 * p for p in c2tan[1]]
#operations.translate(curves[0], c2tanvec, inplace=True)

# Plot the curves using the curve container
curves.sample_size = 100

# Plot the curve
vis_comp = VisMPL.VisCurve2D()
curves.vis = vis_comp
curves.render()

