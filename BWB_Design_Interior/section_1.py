from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl

import numpy as np
from math import *
from airfoil_1 import Airfoil


class Section(GeomBase):
    airfoil_name = Input()
    y_airfoil = Input()
    cant_angle = Input()
    dihedral_angle = Input()
    eq_pts_le_ss_lst = Input()
    eq_pts_te_ss_lst = Input()
    height_factor = Input()
    nb_eq_pts_rails = Input()
    nb_eq_pts_airfoils = Input()

    @Attribute
    def curve(self):
        """
        Returns the curve of the given section using 'Airfoil' class.
        The leading edge point of the airfoil curve is the point of the line defined by the two points of the
        leading edge rails located on both sides of the y-location of the section at the y-location.
        A same mechanism is used to evaluate the position of the trailing edge point of the airfoil curve. The cant
        angle is then added to reevaluate the position of the trailing edge of the airfoil on the trailing edge rail.

        """
        pln_le = Plane(reference = Point(0, self.y_airfoil, 0), normal = Vector(0, 1, 0))

        p = 0
        while (p < self.nb_eq_pts_rails - 1) and not (
                self.eq_pts_le_ss_lst[p][1] <= self.y_airfoil <= self.eq_pts_le_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_le = LineSegment(start = self.eq_pts_le_ss_lst[p],
                                   end = self.eq_pts_le_ss_lst[p + 1]).surface_intersections(
            pln_le)
        pt_le = lst_inter_le[0]['point']
        p = 0
        while (p < self.nb_eq_pts_rails - 1) and not (
                self.eq_pts_te_ss_lst[p][1] <= self.y_airfoil <= self.eq_pts_te_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_te_no_cant = LineSegment(start = self.eq_pts_te_ss_lst[p],
                                           end = self.eq_pts_te_ss_lst[p + 1]).surface_intersections(pln_le)
        pt_te_no_cant = lst_inter_te_no_cant[0]['point']
        pt_te_projected = pt_te_no_cant.translate(Vector(0, 1, 0),
                                                  np.sin(np.deg2rad(self.cant_angle)) * pt_le.distance(
                                                      pt_te_no_cant))

        pln_te = Plane(reference = pt_te_projected,
                       normal = Vector(pt_le[2] - pt_te_projected[2] - (pt_le[1] - pt_te_projected[1]),
                                       pt_le[0] - pt_te_projected[0], 0))

        p = 0
        while (p < self.nb_eq_pts_rails - 1) and not (
                self.eq_pts_te_ss_lst[p][1] <= pt_te_projected[1] <= self.eq_pts_te_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_te = LineSegment(start = self.eq_pts_te_ss_lst[p],
                                   end = self.eq_pts_te_ss_lst[p + 1]).surface_intersections(
            pln_te)
        pt_te = lst_inter_te[0]['point']

        return Airfoil(airfoil_name = self.airfoil_name, pt_le = pt_le, pt_te = pt_te,
                       dihedral_angle = self.dihedral_angle, height_factor = self.height_factor,
                       nb_eq_pts_airfoils = self.nb_eq_pts_airfoils).airfoil, pt_le.distance(pt_te)

    @Part  # the camber of the airfoil is accounted. Any curve is allowed. It includes
    # avl_section_by_points
    def avl_section(self):
        """
        Returns the AVL section corresponding to the section.
        """
        return avl.SectionFromCurve(curve_in = self.curve[0])
