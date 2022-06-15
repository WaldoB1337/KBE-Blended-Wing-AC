from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl


import numpy as np
from math import *
from airfoil_1 import Airfoil
from section_1 import Section


class SubSection(GeomBase):
    y_airfoil = Input()
    nb_eq_pts_rails = Input()
    nb_eq_pts_airfoils = Input()
    eq_pts_le_ss_lst = Input()
    eq_pts_te_ss_lst = Input()
    crv_defined_airfoils_ss_1 = Input()
    crv_defined_airfoils_ss_2 = Input()
    y_defined_airfoils = Input()
    cant_angles_airfoils = Input()
    dihedral_angles_airfoils = Input()
    chords_defined_airfoils = Input()
    j = Input()

    @Attribute
    def put_pts_upper_lower_lst(self):
        """
        Returns the upper and lower curves corresponding to the upper and lower surface of the intermediate airfoil or
        subsection.
        This enables to have a clear idea of the locations of the leading and trailing edge points of the airfoil, as
        well as a same number of points on both side of the airfoil when a '.equispaced_points' method is applied
        (further information thereafter).

        """
        pts_airfoil_1_lst = self.crv_defined_airfoils_ss_1.equispaced_points(self.nb_eq_pts_airfoils)
        pts_airfoil_2_lst = self.crv_defined_airfoils_ss_2.equispaced_points(self.nb_eq_pts_airfoils)
        x_max = -1e3
        index_max = 0
        for i in range(self.nb_eq_pts_airfoils):
            if pts_airfoil_1_lst[i][0] >= x_max:
                x_max = pts_airfoil_1_lst[i][0]
                index_max = i
        pts_airfoil_1_lst_lower = pts_airfoil_1_lst[:index_max + 1]
        crv_airfoil_1_lst_lower = FittedCurve(pts_airfoil_1_lst_lower)
        pts_airfoil_1_lst_upper = pts_airfoil_1_lst[index_max:]
        crv_airfoil_1_lst_upper = FittedCurve(pts_airfoil_1_lst_upper)

        x_max = -1e3
        index_max = 0
        for i in range(self.nb_eq_pts_airfoils):
            if pts_airfoil_2_lst[i][0] >= x_max:
                x_max = pts_airfoil_2_lst[i][0]
                index_max = i
        pts_airfoil_2_lst_lower = pts_airfoil_2_lst[:index_max + 1]
        crv_airfoil_2_lst_lower = FittedCurve(pts_airfoil_2_lst_lower)

        pts_airfoil_2_lst_upper = pts_airfoil_2_lst[index_max:]
        crv_airfoil_2_lst_upper = FittedCurve(pts_airfoil_2_lst_upper)

        return crv_airfoil_1_lst_lower, crv_airfoil_1_lst_upper, crv_airfoil_2_lst_lower, crv_airfoil_2_lst_upper

    @Attribute
    def curve(self):
        """
        Returns the curve of the given subsection airfoil by interpolating the coordinates of the sections on both sides
        of the airfoil with the y-distance separating them.
        The advantage of having a curve for the upper and lower surface of the airfoil is the ability to distribute
        an equal number of points on both sides, which enables to be sure to consider the right points at section j and
        j+1 for the interpolation. Indeed, if a global '.equispaced_points' method along all the airfoil curve is
        directly applied, the difference of arc lengths between the upper and lower surface will create a shifting of
        the points that will distort the linear interpolation between 3 points (of the j section, the airfoil and j+1
        section) that should be on the same spanwise line.

        As for the 'Section.curve' function (see for more information), the leading and trailing edge point of the
        subsection airfoil are found on the rails. Then, the coordinates of the airfoils are interpolated from the
        sections j and j+1 on both sides. The x-coordinates are finally scaled with the chord of the airfoil dived by
        the chord of the directly linearly interpolated airfoil.
        Finally and as for 'Airfoil.airfoil' function (see for more information), the airfoil is rotated along the z,y
        and proper axis to respect the twist, cant and dihedral conditions (which are for the two last linearly
        interpolated from the conditions in section j and j+1).

        """
        crv_airfoil_1_lst_lower = self.put_pts_upper_lower_lst[0]
        pts_airfoil_1_lst_lower = crv_airfoil_1_lst_lower.equispaced_points(int(self.nb_eq_pts_airfoils / 2))
        crv_airfoil_1_lst_upper = self.put_pts_upper_lower_lst[1]
        pts_airfoil_1_lst_upper = crv_airfoil_1_lst_upper.equispaced_points(int(self.nb_eq_pts_airfoils / 2))
        crv_airfoil_2_lst_lower = self.put_pts_upper_lower_lst[2]
        pts_airfoil_2_lst_lower = crv_airfoil_2_lst_lower.equispaced_points(int(self.nb_eq_pts_airfoils / 2))
        crv_airfoil_2_lst_upper = self.put_pts_upper_lower_lst[3]
        pts_airfoil_2_lst_upper = crv_airfoil_2_lst_upper.equispaced_points(int(self.nb_eq_pts_airfoils / 2))

        pts_airfoil_1_lst = self.crv_defined_airfoils_ss_1.equispaced_points(self.nb_eq_pts_airfoils)
        pts_airfoil_2_lst = self.crv_defined_airfoils_ss_2.equispaced_points(self.nb_eq_pts_airfoils)

        crv_airfoil_lst = []
        pts_airfoil_lst = []
        cant_angle = self.cant_angles_airfoils[self.j] + ((self.y_airfoil - self.y_defined_airfoils[self.j]) / (
                self.y_defined_airfoils[self.j + 1] - self.y_defined_airfoils[self.j])) * (
                             self.cant_angles_airfoils[self.j + 1] - self.cant_angles_airfoils[self.j])
        dihedral_angle = self.dihedral_angles_airfoils[self.j] + ((self.y_airfoil - self.y_defined_airfoils[self.j]) / (
                self.y_defined_airfoils[self.j + 1] - self.y_defined_airfoils[self.j])) * (
                                 self.dihedral_angles_airfoils[self.j + 1] - self.dihedral_angles_airfoils[self.j])

        pln_le = Plane(reference = Point(0, self.y_airfoil, 0), normal = Vector(0, 1, 0))
        p = 0
        # print("eq_pts_le_ss_lst=",eq_pts_le_ss_lst)
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
                                                  np.sin(np.deg2rad(cant_angle)) * pt_le.distance(
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
        # display([crv_le_ss,crv_te_ss,pts])

        chord = -pt_le.distance(pt_te)

        chord_1 = self.chords_defined_airfoils[self.j]
        chord_2 = self.chords_defined_airfoils[self.j+1]

        for k1 in range(int(self.nb_eq_pts_airfoils / 2)):
            z_airfoil_interp = (pts_airfoil_1_lst_lower[k1][2]) + (
                    (self.y_airfoil - pts_airfoil_1_lst_lower[k1][1]) / (
                    pts_airfoil_2_lst_lower[k1][1] - pts_airfoil_1_lst_lower[k1][1])) * (
                                       (pts_airfoil_2_lst_lower[k1][2]) -
                                       (pts_airfoil_1_lst_lower[k1][2]))

            x_airfoil_interp = (pts_airfoil_1_lst_lower[k1][0]) + (
                    (self.y_airfoil - pts_airfoil_1_lst_lower[k1][1]) / (
                    pts_airfoil_2_lst_lower[k1][1] - pts_airfoil_1_lst_lower[k1][1])) * (
                                       (pts_airfoil_2_lst_lower[k1][0]) -
                                       (pts_airfoil_1_lst_lower[k1][0]))

            chord_airfoil_interp = chord_1 + (
                    (self.y_airfoil - pts_airfoil_1_lst[k1][1]) / (
                    pts_airfoil_2_lst[k1][1] - pts_airfoil_1_lst[k1][1])) * (chord_2 - chord_1)

            z_airfoil_linear = z_airfoil_interp
            x_airfoil_linear = x_airfoil_interp / np.abs(chord_airfoil_interp) * np.abs(chord)

            pt = Point(x_airfoil_linear, self.y_airfoil, z_airfoil_linear)

            pts_airfoil_lst.append(pt)

        for k2 in range(int(self.nb_eq_pts_airfoils / 2)):
            z_airfoil_interp = (pts_airfoil_1_lst_upper[k2][2]) + (
                    (self.y_airfoil - pts_airfoil_1_lst_upper[k2][1]) / (
                    pts_airfoil_2_lst_upper[k2][1] - pts_airfoil_1_lst_upper[k2][1])) * (
                                       (pts_airfoil_2_lst_upper[k2][2]) -
                                       (pts_airfoil_1_lst_upper[k2][2]))

            x_airfoil_interp = (pts_airfoil_1_lst_upper[k2][0]) + (
                    (self.y_airfoil - pts_airfoil_1_lst_upper[k2][1]) / (
                    pts_airfoil_2_lst_upper[k2][1] - pts_airfoil_1_lst_upper[k2][1])) * (
                                       (pts_airfoil_2_lst_upper[k2][0]) -
                                       (pts_airfoil_1_lst_upper[k2][0]))

            chord_airfoil_interp = chord_1 + (
                    (self.y_airfoil - pts_airfoil_1_lst[k2][1]) / (
                    pts_airfoil_2_lst[k2][1] - pts_airfoil_1_lst[k2][1])) * (chord_2 - chord_1)

            z_airfoil_linear = z_airfoil_interp
            x_airfoil_linear = x_airfoil_interp / chord_airfoil_interp * np.abs(chord)

            pt = Point(x_airfoil_linear, self.y_airfoil, z_airfoil_linear)

            if k2 == 0:
                pt_le_interp = pt
            pts_airfoil_lst.append(pt)

        # display([pts_airfoil_1_lst_lower,pts_airfoil_1_lst_upper,pts_airfoil_lst,pts_airfoil_2_lst_lower,pts_airfoil_2_lst_upper])
        crv_airfoil = FittedCurve(pts_airfoil_lst)
        crv_airfoil_lst.append(crv_airfoil)
        crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = Position(pt_le_interp),
                                       to_position = Position(pt_le))
        crv_airfoil_lst.append(crv_airfoil)
        rot_angle_y = np.arcsin((crv_airfoil.start[2] - pt_te[2]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(0, 1, 0), angle = rot_angle_y)
        # crv_airfoil_lst.append(crv_airfoil)
        rot_angle_z = np.arcsin((crv_airfoil.start[1] - pt_le[1]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(0, 0, 1), angle = rot_angle_z)
        # crv_airfoil_lst.append(crv_airfoil)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(pt_te[0] - pt_le[0], pt_te[1] - pt_le[1],
                                                   pt_te[2] - pt_le[2]),
                                   angle = -np.deg2rad(dihedral_angle))
        crv_airfoil_lst.append(crv_airfoil)
        # display([pts_airfoil_1_lst_lower,pts_airfoil_1_lst_upper,crv_airfoil_lst,pts_airfoil_2_lst_lower,pts_airfoil_2_lst_upper])

        return crv_airfoil, crv_airfoil_lst, pts_airfoil_lst, pt_le.distance(pt_te)

    @Part  # the camber of the airfoil is accounted. Any curve is allowed. It includes
    # avl_section_by_points
    def avl_section(self):
        """
        Returns the AVL section corresponding to the subsection.
        """
        return avl.SectionFromCurve(curve_in = self.curve[0])
