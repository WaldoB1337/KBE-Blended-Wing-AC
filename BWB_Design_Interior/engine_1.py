from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from copy import deepcopy
from parapy.gui import display

import numpy as np
from math import *

from surface_1 import Surface_BWB


class Engine(GeomBase):
    fact_engine = Input()
    AoA_engine = Input()
    x_vb = Input()
    y_vb = Input()
    chords_vb = Input()
    span_vb = Input()
    crvs_surface_to_points_vb = Input()
    position_surf1_on_surf2_vb = Input()
    nb_eq_pts_airfoils = Input()
    is_there_engine = Input()


    @Attribute
    def calculate_x_y_engine(self):
        """
        Returns the position of the engine in the (OXY) plane.
        The engine is placed on the already defined vertical base. The good positioning of the engine with respect to
        the vertical base has been determined experimentally (trials and errors) given the geometry of the engine casing
        (read from a .dat file). The positioning also depends on a scale factor 'fact_engine', the original length of
        the engine being of 120.847 cm.

        """
        print("Engine placement started")
        x_engine = self.crvs_surface_to_points_vb[-1][int(self.nb_eq_pts_airfoils / 2)][0] - 120.847 / (
                    4 * self.fact_engine)  # change here with inputs
        y_engine = self.crvs_surface_to_points_vb[-1][int(self.nb_eq_pts_airfoils / 2)][1]
        return x_engine, y_engine

    @Attribute
    def create_engine_l_casing(self):
        """
        Returns the curves of the engine casing.
        The curves are circles of radius (z_scaled-original_radius/scaling_factor) placed at the location (center)
        (x_scaled+x_engine,y_scaled+y_engine,inf). The engine is then oriented with a certain angle of attack.
        The z-coordinate is chosen to be infinite to be sure to place the engine above any type of generated BWB.
        Then, a downward vertical projection of the vertical tail on the main surface will be operated to place it on
        the BWB. The attachment point on the surface of the BWB is the midpoint of the 4th curve defining the engine.

        """
        [x_engine, y_engine] = self.calculate_x_y_engine
        with open("coord_engine.dat", 'r') as f:
            crvs_engine_rot = []
            center_engine = Point(120.847 / (2 * self.fact_engine) + x_engine, y_engine / self.fact_engine, 0)
            for line in f:
                x, y, z = line.split(' ', 2)
                x = float(x) / self.fact_engine + x_engine
                y = float(y) / self.fact_engine + y_engine
                z = float(z.strip()) / self.fact_engine
                crv = Circle(z - 40 / self.fact_engine, Point(x, 0, 0))
                crv = RotatedCurve(crv, Point(x, 0, 0), Vector(0, 1, 0), -np.pi / 2)
                crv = RotatedCurve(crv, center_engine, Vector(0, 1, 0), np.deg2rad(self.AoA_engine))
                crv = TranslatedCurve(curve_in = crv, displacement = Vector(0, 0, 1e10))
                crv = TranslatedCurve(curve_in = crv, displacement = Vector(0, y, 0))
                crvs_engine_rot.append(crv)

            pt_attachment_engine = crvs_engine_rot[4].midpoint
            crvs_engine_l = crvs_engine_rot
            for i in range(len(crvs_engine_l)):
                crvs_engine_l[i] = RotatedCurve(crvs_engine_rot[i], pt_attachment_engine, Vector(0, 0, 1), -np.pi)

        return crvs_engine_l

    @Attribute
    def create_base_projection_engine(self):
        """
        Returns the rectangular projection base that is placed at an adequate location in the vertical base to
        serve as a reference surface to weld the engine (especially its attachment point). The appropriate location of
        the center of the rectangle is :
        - x-coordinate : The leading edge x location of tip airfoil of the vertical base
        - y-coordinate : The y location of the vertical base, i.e. the y-location of the leading-edge of root airfoil
        of the vertical base
        - z-coordinate : The leading edge z location of tip airfoil of the vertical base which is lowered of around 6%
        of the vertical base height to avoid leakages between the engine and the vertical base tip.

        """
        pt_center_rect_placement_engine = Point(self.crvs_surface_to_points_vb[-1][int(self.nb_eq_pts_airfoils / 2)][0],
                                                self.y_vb,
                                                self.crvs_surface_to_points_vb[-1][int(self.nb_eq_pts_airfoils / 2)][
                                                    2] - self.span_vb / 17)
        surface_base_rect_placement_engine = RectangularSurface(width = self.chords_vb[0] * 2,
                                                                length = self.chords_vb[0],
                                                                position = Position(pt_center_rect_placement_engine))
        surface_base_rect_placement_engine = RotatedSurface(surface_base_rect_placement_engine,
                                                            pt_center_rect_placement_engine, Vector(0, 1, 0),
                                                            -np.deg2rad(self.position_surf1_on_surf2_vb[1]))
        return surface_base_rect_placement_engine

    @Attribute
    def define_attachment_pt_object(self):
        return self.create_engine_l_casing[4].midpoint

    @Attribute
    def position_surf1_on_surf2(self):
        """
        Returns the projected curves of the engine on the surface of the wing.

        """
        pt_attachment_engine = self.define_attachment_pt_object
        pts_surface_vb_lst = self.crvs_surface_to_points_vb
        crvs_engine_l = self.create_engine_l_casing[:]
        surface_base_rect_placement_engine = self.create_base_projection_engine
        lst_inter_te = LineSegment(start = pt_attachment_engine,
                                   end = Point(pt_attachment_engine[0], pt_attachment_engine[1],
                                               -5)).surface_intersections(
            surface_base_rect_placement_engine)
        pt_attachment_engine_projected = lst_inter_te[0]['point']
        vect_displ_engine = Vector(pt_attachment_engine_projected[0] - pt_attachment_engine[0],
                                   pt_attachment_engine_projected[1] - pt_attachment_engine[1],
                                   pt_attachment_engine_projected[2] - pt_attachment_engine[2])

        for i in range(len(crvs_engine_l)):
            crvs_engine_l[i] = TranslatedCurve(curve_in = crvs_engine_l[i], displacement = vect_displ_engine)

        return crvs_engine_l

    @Attribute
    def mirror_crvs(self):
        crvs_engine_l = self.position_surf1_on_surf2
        crvs_engine_r = []
        for i in range(len(crvs_engine_l)):
            crvs_engine_r.append(
                MirroredCurve(curve_in = crvs_engine_l[i], reference_point = Point(0, 0, 0), vector1 = Vector(1, 0, 0),
                              vector2 = Vector(0, 0, 1)))
        print("Engine placed")
        return crvs_engine_r

    @Part
    def surface(self):
        """
        Returns the fusion of the lofted curves defining the cowl of the casing and of the other lofted curves defining
        the outlet of the engine.

        """
        return Fused(LoftedSurface(profiles = self.position_surf1_on_surf2[0:3], min_degree = 25),
                     LoftedSurface(profiles = self.position_surf1_on_surf2[2:], min_degree = 25))

    @Part
    def surface_mirrored(self):
        return Fused(LoftedSurface(profiles = self.mirror_crvs[0:3], min_degree = 25),
                     LoftedSurface(profiles = self.mirror_crvs[2:], min_degree = 25))


if __name__ == '__main__':
    from parapy.gui import display

    obj = Engine()
    display(obj)
