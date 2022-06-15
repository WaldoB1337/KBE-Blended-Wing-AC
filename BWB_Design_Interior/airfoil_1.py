from parapy.core import *
from parapy.geom import *
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve



import numpy as np
from math import *

class Airfoil(GeomBase):
    airfoil_name = Input()
    pt_le = Input()
    pt_te = Input()
    dihedral_angle = Input()
    height_factor = Input()
    nb_eq_pts_airfoils = Input()

    @Attribute
    def airfoil(self):
        """

        Returns the curve of the airfoil 'airfoil_name' (NACA 4,5 or .dat file) placed at a location so that 'pt_le' and 'pt_te' are the
        leading and trailing edge points of the airfoil.
        It considers the dihedral angle and height factor provided.

        How it works?
        - NACA cases:
        An initial curve is created using 'Naca4/5AirfoilCurve' function. The height factor is applied by translating
        any points equispaced on the curve along the z-axis (no scaling). The obtained curve is then translated to the
        leading edge point and scaled with the given chord.
        - .dat file case:
        A curve is obtained by reading and fitting the points from the files. The surface is split between upper and
        lower parts to better fit the points.
        The airfoil is then rotated around the z-axis (cant angle), y-axis (twist angle), airfoil axis (dihedral angle)

        """
        # airfoil_name = "whitcomb"
        # pt_le = Point(0,0,0)
        # pt_te = Point(1,0,2)
        chord = -self.pt_le.distance(self.pt_te)

        if len(self.airfoil_name) == 4 and (type(int(char) == int) for char in self.airfoil_name):
            #print('NACA4')
            crv_airfoil = DynamicType(type = Naca4AirfoilCurve,
                                      designation = self.airfoil_name,
                                      mesh_deflection = 0.00001,
                                      hidden = True)
            if self.height_factor != 1:
                pts_airfoil = crv_airfoil.equispaced_points(self.nb_eq_pts_airfoils)
                pts_airfoil_wthout_fact = pts_airfoil
                for i in range(self.nb_eq_pts_airfoils):
                    pts_airfoil[i] = pts_airfoil_wthout_fact[i].translate(Vector(0, 0, 1), self.height_factor *
                                                                          pts_airfoil_wthout_fact[i][2])
                crv_airfoil = FittedCurve(pts_airfoil)
            crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = OXY, to_position = Position(self.pt_le))
            crv_airfoil = ScaledCurve(crv_airfoil,
                                      self.pt_le,
                                      chord,
                                      mesh_deflection = 0.00001)

        if len(self.airfoil_name) == 5 and (type(int(char) == int) for char in self.airfoil_name):
            #print('NACA5')
            crv_airfoil = DynamicType(type = Naca5AirfoilCurve,
                                      designation = self.airfoil_name,
                                      mesh_deflection = 0.00001,
                                      hidden = True)
            if self.height_factor != 1:
                pts_airfoil = crv_airfoil.equispaced_points(self.nb_eq_pts_airfoils)
                pts_airfoil_wthout_fact = pts_airfoil
                for i in range(self.nb_eq_pts_airfoils):
                    pts_airfoil[i] = pts_airfoil_wthout_fact[i].translate(Vector(0, 0, 1), self.height_factor *
                                                                          pts_airfoil_wthout_fact[i][2])
                crv_airfoil = FittedCurve(pts_airfoil)
            crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = OXY, to_position = Position(self.pt_le))
            crv_airfoil = ScaledCurve(crv_airfoil,
                                      self.pt_le,
                                      chord,
                                      mesh_deflection = 0.00001
                                      )

        if not ((len(self.airfoil_name) == 4 or len(self.airfoil_name) == 5) and (type(int(char) == int) for char in
                                                                        self.airfoil_name)):
            #print('Airfoil specified')
            with open(self.airfoil_name + ".dat", 'r') as f:
                pts_airfoil_lst = []
                for line in f:
                    x, z = line.split(' ',
                                      1)  # the cartesian coordinates are directly interpreted as X and Z coordinates
                    pts_airfoil_lst.append(self.pt_le.translate(
                        "x", float(x) * chord,  # the x points are scaled according to the airfoil chord length
                        "z", float(z) * chord * self.height_factor))  # the y points are scaled according to the
            # Creation of a curve te->le(up)+le->te(lo)

            idx = pts_airfoil_lst.index(self.pt_le)
            upper_pts_lst = pts_airfoil_lst[:idx + 1]
            lower_pts_lst = pts_airfoil_lst[idx:]

            crv_upper = FittedCurve(points = upper_pts_lst, max_degree = 8)
            crv_lower = FittedCurve(points = lower_pts_lst, max_degree = 8)
            crv_airfoil = ComposedCurve(built_from = [crv_upper, crv_lower], allow_multiple = 'True')


        rot_angle_z = np.arcsin((self.pt_te[1] - self.pt_le[1]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = self.pt_le,
                                   vector = Vector(0, 0, 1), angle = rot_angle_z)
        rot_angle_y = np.arcsin((self.pt_te[2] - self.pt_le[2]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = self.pt_le,
                                   vector = Vector(0, 1, 0), angle = -rot_angle_y)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = self.pt_le,
                                   vector = Vector(self.pt_te[0] - self.pt_le[0], self.pt_te[1] - self.pt_le[1], self.pt_te[2] - self.pt_le[2]),
                                   angle = -np.deg2rad(self.dihedral_angle))


        # display(crv_airfoil)
        return crv_airfoil