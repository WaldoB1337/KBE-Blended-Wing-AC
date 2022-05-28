from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from parapy.gui import display

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
    j = Input()
    crv_defined_airfoils_ss_1 = Input()
    crv_defined_airfoils_ss_2 = Input()
    y_defined_airfoils = Input()
    cant_angles_airfoils = Input()
    dihedral_angles_airfoils = Input()

    @Attribute
    def put_pts_upper_lower_lst(self):
        pts_airfoil_1_lst = self.crv_defined_airfoils_ss_1.equispaced_points(self.nb_eq_pts_airfoils)
        pts_airfoil_2_lst = self.crv_defined_airfoils_ss_2.equispaced_points(self.nb_eq_pts_airfoils)
        x_max = -1e3
        index_max = 0
        for i in range(self.nb_eq_pts_airfoils):
            if pts_airfoil_1_lst[i][0] >= x_max:
                x_max = pts_airfoil_1_lst[i][0]
                index_max = i
        pts_airfoil_1_lst_lower = pts_airfoil_1_lst[:index_max+1]
        crv_airfoil_1_lst_lower = FittedCurve(pts_airfoil_1_lst_lower)
        pts_airfoil_1_lst_upper = pts_airfoil_1_lst[index_max:]
        crv_airfoil_1_lst_upper = FittedCurve(pts_airfoil_1_lst_upper)

        x_max = -1e3
        index_max = 0
        for i in range(self.nb_eq_pts_airfoils):
            if pts_airfoil_2_lst[i][0] >= x_max:
                x_max = pts_airfoil_2_lst[i][0]
                index_max = i
        pts_airfoil_2_lst_lower = pts_airfoil_2_lst[:index_max]
        crv_airfoil_2_lst_lower = FittedCurve(pts_airfoil_2_lst_lower)

        pts_airfoil_2_lst_upper = pts_airfoil_2_lst[index_max:]
        crv_airfoil_2_lst_upper = FittedCurve(pts_airfoil_2_lst_upper)

        return crv_airfoil_1_lst_lower, crv_airfoil_1_lst_upper, crv_airfoil_2_lst_lower, crv_airfoil_2_lst_upper

    @Attribute
    def curve(self):
        crv_airfoil_1_lst_lower = self.put_pts_upper_lower_lst[0]
        pts_airfoil_1_lst_lower = crv_airfoil_1_lst_lower.equispaced_points(int(self.nb_eq_pts_airfoils/2))
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
        if not lst_inter_te == []:
            pt_te = lst_inter_te[0]['point']
        # display([crv_le_ss,crv_te_ss,pts])

        chord = -pt_le.distance(pt_te)

        chord_interp = -pt_le.distance(pt_te_no_cant)

        # te->le(up)+le->te(lo) for x
        x_upper = np.linspace(pt_te[0], pt_le[0], int(self.nb_eq_pts_airfoils / 2))
        x_lower = np.linspace(pt_le[0], pt_te[0], int(self.nb_eq_pts_airfoils / 2))
        x_airfoil = np.concatenate((x_upper, x_lower))

        for k1 in range(int(self.nb_eq_pts_airfoils / 2)):
            z_airfoil_interp = pts_airfoil_1_lst_lower[k1][2] + (
                    (self.y_airfoil - pts_airfoil_1_lst_lower[k1][1]) / (
                    pts_airfoil_2_lst_lower[k1][1] - pts_airfoil_1_lst_lower[k1][1])) * (
                                       pts_airfoil_2_lst_lower[k1][2] -
                                       pts_airfoil_1_lst_lower[k1][2])
            z_airfoil_norm = z_airfoil_interp / np.abs(chord_interp)
            z_airfoil_linear = z_airfoil_norm * np.abs(chord)

            pt = Point(x_airfoil[k1], self.y_airfoil, z_airfoil_linear)

            pts_airfoil_lst.append(pt)

        for k2 in range(int(self.nb_eq_pts_airfoils / 2)):
            z_airfoil_interp = pts_airfoil_1_lst_upper[k2][2] + (
                    (self.y_airfoil - pts_airfoil_1_lst_upper[k2][1]) / (
                    pts_airfoil_2_lst_upper[k2][1] - pts_airfoil_1_lst_upper[k2][1])) * (
                                       pts_airfoil_2_lst_upper[k2][2] -
                                       pts_airfoil_1_lst_upper[k2][2])
            z_airfoil_norm = z_airfoil_interp / np.abs(chord_interp)
            z_airfoil_linear = z_airfoil_norm * np.abs(chord)

            pt = Point(x_airfoil[k2+k1+1], self.y_airfoil, z_airfoil_linear)
            if k2 == 0:
                    pt_le_interp = pt
            pts_airfoil_lst.append(pt)

        crv_airfoil = FittedCurve(pts_airfoil_lst)
        #crv_airfoil_lst.append(crv_airfoil)
        crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = Position(pt_le_interp),
                                       to_position = Position(pt_le))
        #crv_airfoil_lst.append(crv_airfoil)
        rot_angle_y = np.arcsin((crv_airfoil.start[2] - pt_te[2]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(0, 1, 0), angle = rot_angle_y)
        #crv_airfoil_lst.append(crv_airfoil)
        rot_angle_z = np.arcsin((crv_airfoil.start[1] - pt_le[1]) / chord)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(0, 0, 1), angle = rot_angle_z)
        #crv_airfoil_lst.append(crv_airfoil)
        crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                   vector = Vector(pt_te[0] - pt_le[0], pt_te[1] - pt_le[1],
                                                   pt_te[2] - pt_le[2]),
                                   angle = -np.deg2rad(dihedral_angle))
        #crv_airfoil_lst.append(crv_airfoil)

        return crv_airfoil
