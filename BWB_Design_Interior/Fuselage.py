from kbeutils.geom import PositionedFittedCurve
from parapy.core import *
from parapy.geom import *
import numpy as np


class Fuselage(GeomBase):
    length = Input()
    diameter = Input()
    n_points = Input(100)

    @Attribute
    def eta_ordinate(self):
        return np.linspace(0, 1, self.n_points)

    @Attribute
    def x_ordinate(self):
        return [ordinate*self.length for ordinate in self.eta_ordinate]

    @Attribute
    def sears_haack_radius(self):
        return [self.diameter/2*(4*eta*(1-eta))**(3/4) for eta in self.eta_ordinate]

    @Attribute
    def crown_coordinates(self):
        x = self.x_ordinate
        # y = np.zeros_like(x)
        y = [0]*len(x)
        z = self.sears_haack_radius
        return list(zip(x, y, z))

    @Part
    def crown_curve(self):
        return PositionedFittedCurve(coordinates=self.crown_coordinates)

    @Part
    def surface(self):
        return RevolvedSurface(basis_curve=self.crown_curve,
                               center=self.position.point,
                               direction=self.position.Vx,
                               mesh_deflection=0.0001)

    # @Part
    # def surface1(self):
    #     return RevolvedSurface(basis_curve=self.crown_curve,
    #                            center=self.crown_curve.position.point,
    #                            direction=self.crown_curve.position.Vx,
    #                            mesh_deflection=0.0001)
    #
    # @Part
    # def surface2(self):
    #     return RevolvedSurface(basis_curve=self.crown_curve,
    #                            center=self.crown_curve.start,
    #                            direction=self.crown_curve.start - self.crown_curve.end,
    #                            mesh_deflection=0.0001)


if __name__ == '__main__':
    from parapy.gui import display
    obj = Fuselage(diameter=5,
                   length=30)
    display(obj)