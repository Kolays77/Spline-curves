import numpy as np
import matplotlib.pyplot as plt
 
 
class BaseSpline:
    def __init__(self, p: int, points, knots, weights=None, closed=False):
        """The Base class Spline stores spline data 
        and  reports an error with incorrect data.

        Args:
            p (int): degree of Spline
            points (list): control points
            knots (list): knot vector

        Raises:
            ValueError: Checking for degree positivity
            ValueError: Checking for the number of points and degrees
            ValueError: Checking the length of the knot  vector
        """
        if p >= (len(points)):
            raise ValueError(
                "Degree must be less than the number of control points.")

        if not (p > 0 and isinstance(p, int)):
            raise ValueError(
                "Spline degree must be a positive integer, not %s." % p)
        
        if knots is not None:
            if len(knots) != p + len(points) + 1:
                raise ValueError(
                "The length Spline's knots must be equal p + n + 1, not %s." % p)

        self.p = p
  
        self.points = np.array(points, dtype=np.float128)
        self.n = len(points)
      
        self.create_weights(weights)
        self.create_knots(knots)
        self.points = np.array([p*w for p, w in zip(self.points, self.weights)])
            
        
        # if closed:
        #     self.points = np.concatenate((self.points,self.points[:p]))
        #     self.knots = list(range(len(self.points) + self.p + 1))    
        #     self.weights = np.ones(len(self.points))
        #     self.n = len(self.points)
        # else:
        #     self.n = len(self.points)
        #     self.create_knots(knots)
            
        self.dim = len(self.points.shape)        
        self.t_start, self.t_end = self.knots[self.p], self.knots[-self.p - 1]

    def create_weights(self, weights):
        if weights is not None:
            if len(weights) != self.n:
                raise ValueError(
                f"The length {self.n}  of points must be equal length of weights {len(weights)}.")
            else:
                self.weights = weights
        else:
            self.weights = np.ones(self.n)  
    
    def create_knots(self, knots):
        if knots is None:
            interiors = [(i)/(self.n - self.p) for i in range(1, self.n - self.p)]
            knots = [0]*(self.p+1) + interiors + [1]*(self.p+1)
        self.knots = np.array(knots, dtype=np.float128)
        
    def construct_data_to_derivative(self):
        n = self.n - 1
        new_points = []
        new_knots = self.knots[1:-1]
        p = self.p - 1
        for i in range(n):
            new_points.append(self.p * (self.points[i+1] - self.points[i]) /
                              (self.knots[i+self.p + 1] - self.knots[i+1]))

        return (p, new_points, new_knots)