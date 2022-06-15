from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl

from parapy.gui import display

import numpy as np
from math import *

from section_1 import Section
from sub_section_3 import SubSection
from airfoil_1 import Airfoil


class Surface_BWB(GeomBase):
    nb_eq_pts_rails = Input()
    nb_eq_pts_airfoils = Input()

    nb_fixed_control_points = Input()
    sweeps_le = Input()
    sweeps_te = Input()
    sspan = Input()
    x_fcp_except_last = Input()
    y_fcp_except_last = Input()
    chords = Input()
    rho_le = Input()
    rho_te = Input()

    nb_pts_twist = Input()
    angles_twist = Input()  # degrees of twist
    y_twist_except_last = Input()  # location of the definition of twist
    nb_pts_twist_axis = Input()
    x_twist_axis = Input()
    y_twist_axis_except_last = Input()  # 1st and last points must be at the root and tip

    nb_pts_dihedral = Input()
    angles_dihedral = Input()  # angle defined on the right of the section
    spans_dihedral_except_last = Input()

    weights_dihedral = Input()

    nb_defined_airfoils = Input()
    airfoils_names = Input()
    y_defined_airfoils_except_last = Input()  # 1st and last correspond to root and tip
    cant_angles_airfoils = Input()  # degree
    dihedral_angles_airfoils = Input()
    heights_factors = Input()

    nb_airfoils_ss = Input()

    surface_type = Input()
    is_twisted = Input()
    is_there_dihedral = Input()

    name = Input()
    nb_chordwise_vortices = Input()
    nb_spanwise_vortices = Input()
    is_mirrored = Input()
    nb_sections_AVL = Input()

    # tol = 1e-5
    # round_nb = 7

    @Attribute
    def spans_sections_dihedral(self):
        """
        Returns the lengths between each defined airfoil along the semi-span of the surface including the last section

        """
        print("Surface '"+self.surface_type+"' treatment started")
        sum = 0
        for i in range(len(self.spans_dihedral_except_last)):
            sum += self.spans_dihedral_except_last[i]
        return self.spans_dihedral_except_last + [self.sspan - sum]

    @Attribute
    def dihedral_angles_to_y_z(self):
        """
        Returns the z(y) dihedral distribution (along the span) given the dihedral angles filled out by the user

        """
        y_dihedral = [
            0, ]  # location of the definition of dihedral. 1st and last points must be defined at root and tip
        z_dihedral = [0, ]
        #print(y_dihedral)
        for i in range(self.nb_pts_dihedral - 1):
            #print([y_dihedral[-1] + cos(np.deg2rad(self.angles_dihedral[i + 1])) * self.spans_sections_dihedral[i]])
            y_dihedral += [
                y_dihedral[-1] + cos(np.deg2rad(self.angles_dihedral[i+1])) * self.spans_sections_dihedral[i]]
            z_dihedral += [
                z_dihedral[-1] + sin(np.deg2rad(self.angles_dihedral[i+1])) * self.spans_sections_dihedral[i]]
        return y_dihedral, z_dihedral

    @Attribute
    def y_defined_airfoils(self):
        """
        Returns the y-location (along the span) of the defined airfoils (including the last one, being the tip of
        the surface, and thus depending on the tip dihedral). The points coordinates are rounder for further processing
        purposes.

        """
        y_defined_airfoils = self.y_defined_airfoils_except_last + [self.sspan]
        return [round(y_defined_airfoils[i], 7) for i in range(self.nb_defined_airfoils)]

    @Attribute
    def x_fcp(self):
        """
        Returns the x-coordinates of the fixed control points defining the sweep distribution of the surface. If only
        one fixed control point is defined, the sweep of the surface is a line between this point an another defined
        with the condition in sweep given by the user.

        """
        if len(self.x_fcp_except_last) == 1:
            return self.x_fcp_except_last + [-tan(np.deg2rad(self.sweeps_le[0])) * self.sspan]
        return self.x_fcp_except_last

    @Attribute
    def y_fcp(self):
        """
        Returns the y-coordinates of the fixed control points (including the last one, being the tip of the surface, and
        thus depending on the tip dihedral)

        """
        return self.y_fcp_except_last + [self.sspan]

    @Attribute
    def y_twist(self):
        """
        Returns the y-coordinates of the locations of the twist defined by the user (including the last one, being the
        tip of the surface, and thus depending on the tip dihedral)

        """
        return self.y_twist_except_last + [self.sspan]

    @Attribute
    def y_twist_axis(self):
        """
        Returns the y-coordinates of the locations of the twist axis defined by the user (including the last one, being
        the tip of the surface, and thus depending on the tip dihedral)

        """
        return self.y_twist_axis_except_last + [self.sspan]

    @Attribute
    def planform_area(self):
        """
        Returns the area of the surface by approximating the planform surface between each defined airfoil as a
        trapezoidal form.

        """
        area = 0
        for i in range(self.nb_airfoils_ss - 1):
            area += (self.chords_all_airfoils[i] + self.chords_all_airfoils[i + 1]) * (
                        self.determine_position_all_airfoils[i + 1] - self.determine_position_all_airfoils[i]) / 2
        #print('area=', area)
        return area * 2

    @Attribute
    def mac(self):
        """
        Returns the MAC (Mean Average Chord) of the surface by using a formula for the multi-cranked flying surface"

        """
        mac_num = 0
        mac_denum = 0
        for i in range(self.nb_airfoils_ss - 1):
            lambda_i = self.chords_all_airfoils[i + 1] / self.chords_all_airfoils[i]
            zeta_i = (self.determine_position_all_airfoils[i + 1] - self.determine_position_all_airfoils[
                i]) / self.sspan
            mac_num += self.chords_all_airfoils[i] ** 2 * (1 + lambda_i + lambda_i ** 2) * zeta_i
            mac_denum += self.chords_all_airfoils[i] * (1 + lambda_i) * zeta_i
        return 2 / 3 * mac_num / mac_denum

    @Attribute
    def twist_axis(self):
        if self.nb_pts_twist == 2:
            x_twist_axis = [-self.chords[0]/2, -np.tan(np.deg2rad(self.sweeps_le[0]))*self.sspan-self.chords[1]/2]
            y_twist_axis = [0,self.sspan]
            return x_twist_axis,y_twist_axis

    @Attribute
    def sweep_crvs(self):
        """
        Returns 'nb_eq_pts_rails' equally distributed points along the leading and trailing edges curves. These curves
        are created using the conditions in sweep given by the user at the defined fixed control points along the
        surface span.
        By adding weights ('rho_le' and 'rho_te'), the conditions are more or less saturated in the close
        proximity of the variables control points thus defined. Bezier curves are then built with the whole points.
        The final curves are decomposed into points which are equally spaced along the y-axis (the points are thereafter
        processed along the y-axis for twist and dihedral reasons).
        If the surface considered in a vertical base (only 2 fixed control points), the leading and trailing edge curves
        are lines defined between these 2 points.

        """
        # Creation of the sweep curves
        lambda_le = np.deg2rad(self.sweeps_le)
        if self.surface_type == 'vb' or self.surface_type == 'vt':
            crv_le = LineSegment(start = Point(self.x_fcp[0], self.y_fcp[0], 0),
                                 end = Point(self.x_fcp[1], self.y_fcp[1], 0))
            crv_te = LineSegment(start = Point(self.x_fcp[0] - self.chords[0], self.y_fcp[0], 0),
                                 end = Point(self.x_fcp[1] - self.chords[1], self.y_fcp[1], 0))
        else:
            lambda_te = np.deg2rad(self.sweeps_te)
            crvs_le_lst = []
            crvs_te_lst = []

            for i in range(self.nb_fixed_control_points - 1):
                xc_1_le = self.x_fcp[i] + self.rho_le[i] * -sin(lambda_le[i])
                yc_1_le = self.y_fcp[i] + self.rho_le[i] * cos(lambda_le[i])
                xc_1_te = self.x_fcp[i] - self.chords[i] + self.rho_te[i] * -sin(lambda_te[i])
                yc_1_te = self.y_fcp[i] + self.rho_te[i] * cos(lambda_te[i])

                xc_2_le = self.x_fcp[i + 1] + self.rho_le[i + 1] * -sin(lambda_le[i + 1] + pi)
                yc_2_le = self.y_fcp[i + 1] + self.rho_le[i + 1] * cos(lambda_le[i + 1] + pi)
                xc_2_te = self.x_fcp[i + 1] - self.chords[i + 1] + self.rho_te[i + 1] * -sin(lambda_te[i + 1] + pi)
                yc_2_te = self.y_fcp[i + 1] + self.rho_te[i + 1] * cos(lambda_te[i + 1] + pi)

                FCP0 = Point(self.x_fcp[i], self.y_fcp[i], 0)
                VCP0 = Point(xc_1_le, yc_1_le, 0)
                VCP1 = Point(xc_2_le, yc_2_le, 0)
                FCP1 = Point(self.x_fcp[i + 1], self.y_fcp[i + 1], 0)
                pts = [FCP0, VCP0, VCP1, FCP1]
                crv = BezierCurve(pts)

                crvs_le_lst.append(crv)

                FCP0 = Point(self.x_fcp[i] - self.chords[i], self.y_fcp[i], 0)
                VCP0 = Point(xc_1_te, yc_1_te, 0)
                VCP1 = Point(xc_2_te, yc_2_te, 0)
                FCP1 = Point(self.x_fcp[i + 1] - self.chords[i + 1], self.y_fcp[i + 1], 0)
                pts = [FCP0, VCP0, VCP1, FCP1]
                crv = BezierCurve(pts)

                crvs_te_lst.append(crv)

            crv_le = ComposedCurve(built_from = crvs_le_lst, allow_multiple = 'True')
            # display(crv_le)
            crv_te = ComposedCurve(built_from = crvs_te_lst, allow_multiple = 'True')
            # display(crv_te)

        eq_ss_lst = np.linspace(0, self.sspan, self.nb_eq_pts_rails)
        eq_pts_rails_lst = [Point(0, eq_ss_lst[i], 0) for i in range(self.nb_eq_pts_rails)]
        #display([eq_pts_rails_lst,crv_le,crv_te])
        eq_pts_le_ss_lst = []
        eq_pts_te_ss_lst = []
        for i in range(self.nb_eq_pts_rails):
            pt_le = crv_le.intersection_point(Plane(reference = eq_pts_rails_lst[i], normal = Vector(0, 1, 0)))
            eq_pts_le_ss_lst.append(pt_le)
            pt_te = crv_te.intersection_point(Plane(reference = eq_pts_rails_lst[i], normal = Vector(0, 1, 0)))
            eq_pts_te_ss_lst.append(pt_te)
        #display([eq_pts_le_ss_lst, eq_pts_te_ss_lst])
        return eq_pts_le_ss_lst, eq_pts_te_ss_lst

    @Attribute
    def calculate_twist_distrib(self):
        """"
        Returns 'nb_eq_pts_rails' equally distributed points along the y-axis (and not along the span) respecting the
        twist conditions provided by the user. A twist distribution is generated by
        fitting the twist conditions along the span of the surface. The twist axis is defined by joined segments between
        the points given by the user.

        """
        if not self.is_twisted:
            return []
        # Twist distribution
        pts_twist_lst = []
        for i in range(self.nb_pts_twist):
            pt = Point(self.y_twist[i], self.angles_twist[i], 0)
            pts_twist_lst.append(pt)
        crv_twist_distrib = FittedCurve(points = pts_twist_lst)

        # Twist axis definition
        crvs_twist_axis = []
        if self.nb_pts_twist == 2:
            x_twist_axis = self.twist_axis[0]
            y_twist_axis = self.twist_axis[1]
            print(x_twist_axis,y_twist_axis)
        else:
            x_twist_axis = self.x_twist_axis
            y_twist_axis = self.y_twist_axis
        for i in range(self.nb_pts_twist_axis - 1):
            crv = LineSegment(Point(x_twist_axis[i], y_twist_axis[i], 0),
                              Point(x_twist_axis[i + 1], y_twist_axis[i + 1], 0))
            crvs_twist_axis.append(crv)

        crv_twist_axis = ComposedCurve(built_from = crvs_twist_axis)

        eq_ss_lst = np.linspace(0, self.y_twist[-1], self.nb_eq_pts_rails)
        eq_pts_twist_distrib = [Point(eq_ss_lst[i], 0, 0) for i in range(self.nb_eq_pts_rails)]
        eq_pts_twist_axis = [Point(0, eq_ss_lst[i], 0) for i in range(self.nb_eq_pts_rails)]

        eq_pts_twist_distrib_lst = []
        eq_pts_twist_axis_lst = []

        for i in range(self.nb_eq_pts_rails):
            pt_twist_distrib = crv_twist_distrib.intersection_point(
                Plane(reference = eq_pts_twist_distrib[i], normal = Vector(1, 0, 0)))
            eq_pts_twist_distrib_lst.append(pt_twist_distrib)
            pt_twist_axis = crv_twist_axis.intersection_point(
                Plane(reference = eq_pts_twist_axis[i], normal = Vector(0, 1, 0)))
            eq_pts_twist_axis_lst.append(pt_twist_axis)

        return eq_pts_twist_distrib_lst, eq_pts_twist_axis_lst

    @Attribute
    def calculate_dihedral_distrib(self):
        """
        Returns 'nb_eq_pts_rails' equally distributed points along the y-axis (and not along the span) respecting the
        dihedral conditions provided by the user and treated with the 'dihedral_angles_to_y_z' function. The dihedral
        distribution is generated by using Bezier curve which enables to control the slope of the curve more
        exactly and hence to create more realistic dihedral distribution.
        In the region of high dihedral of a winglet for example, the weights can be accentuated to better match the
        smoothness of the transition between the wing and the winglet.

        """
        [y_dihedral, z_dihedral] = self.dihedral_angles_to_y_z
        if not self.is_there_dihedral:
            return []
        if self.nb_pts_dihedral == 2:
            crv_dihedral_distrib = FittedCurve(
                [Point(y_dihedral[0], z_dihedral[0], 0), Point(y_dihedral[1], z_dihedral[1], 0)])
        else:

            pts_dihedral_lst = []
            FCPdi = Point(y_dihedral[0], z_dihedral[0], 0)
            VCPdi = Point(((y_dihedral[1] + y_dihedral[0]) / 2), (z_dihedral[1] + z_dihedral[0]) / 2, 0)
            pts_dihedral_lst.append(FCPdi)
            pts_dihedral_lst.append(VCPdi)

            for i in range(1, self.nb_pts_dihedral - 1):
                VCPd1 = Point(((y_dihedral[i - 1] + y_dihedral[i]) / 2), (z_dihedral[i - 1] + z_dihedral[i]) / 2, 0)
                FCPd = Point(y_dihedral[i], z_dihedral[i], 0)
                VCPd2 = Point(((y_dihedral[i + 1] + y_dihedral[i]) / 2),
                              (z_dihedral[i + 1] + z_dihedral[i]) / 2, 0)
                pts_dihedral_lst.append(VCPd1)
                pts_dihedral_lst.append(FCPd)
                pts_dihedral_lst.append(VCPd2)
            VCPdf = Point(((y_dihedral[-2] + y_dihedral[-1]) / 2), (z_dihedral[-2] + z_dihedral[-1]) / 2, 0)
            FCPdf = Point(y_dihedral[-1], z_dihedral[-1], 0)
            pts_dihedral_lst.append(VCPdf)
            pts_dihedral_lst.append(FCPdf)
            crv_dihedral_distrib = BezierCurve(control_points = pts_dihedral_lst,
                                               weights = self.weights_dihedral)

        eq_ss_lst = np.linspace(0, y_dihedral[-1], self.nb_eq_pts_rails)
        eq_pts_dihedral_distrib = [Point(eq_ss_lst[i], 0, 0) for i in range(self.nb_eq_pts_rails)]
        eq_pts_dihedral_distrib_lst = []

        for i in range(self.nb_eq_pts_rails):
            pt_dihedral = crv_dihedral_distrib.intersection_point(
                Plane(reference = eq_pts_dihedral_distrib[i], normal = Vector(1, 0, 0)))
            eq_pts_dihedral_distrib_lst.append(pt_dihedral)

        return eq_pts_dihedral_distrib_lst, crv_dihedral_distrib

    @Attribute
    def add_twist(self):
        """
        Returns the updated 'eq_pts_le_ss_lst' and 'eq_pts_te_ss_lst' points defining the leading and trailing edges of
        the sweeped surface by adding the twist distribution to the model with respect to the twist axis. If the surface
        is not twisted, the already sweeped points defined are not changed.

        """
        [eq_pts_le_ss_lst, eq_pts_te_ss_lst] = [self.sweep_crvs[0][:],self.sweep_crvs[1][:]]
        if not self.is_twisted:
            return eq_pts_le_ss_lst, eq_pts_te_ss_lst
        [eq_pts_twist_distrib_lst, eq_pts_twist_axis_lst] = self.calculate_twist_distrib

        for i in range(self.nb_eq_pts_rails):
            if self.nb_eq_pts_rails==2 :
                eq_pts_le_ss_lst[i] = eq_pts_le_ss_lst[i].rotate_around(eq_pts_twist_axis_lst[i],
                    Vector(0,np.cos(np.deg2rad(self.dihedral_angles_airfoils[i])),np.sin(np.deg2rad(self.dihedral_angles_airfoils[i]))), np.deg2rad(eq_pts_twist_distrib_lst[i][1]))
                eq_pts_te_ss_lst[i] = eq_pts_te_ss_lst[i].rotate_around(eq_pts_twist_axis_lst[i],
                                                                        Vector(0, np.cos(np.deg2rad(
                                                                            self.dihedral_angles_airfoils[i])), np.sin(
                                                                            np.deg2rad(
                                                                                self.dihedral_angles_airfoils[i]))),
                                                                        np.deg2rad(eq_pts_twist_distrib_lst[i][1]))
            else :
                eq_pts_le_ss_lst[i] = eq_pts_le_ss_lst[i].rotate_around(eq_pts_twist_axis_lst[i],
                                                                        Vector(0,1,0),
                                                                        np.deg2rad(eq_pts_twist_distrib_lst[i][1]))
                eq_pts_te_ss_lst[i] = eq_pts_te_ss_lst[i].rotate_around(eq_pts_twist_axis_lst[i],
                                                                        Vector(0,1,0),
                                                                        np.deg2rad(eq_pts_twist_distrib_lst[i][1]))

        return eq_pts_le_ss_lst, eq_pts_te_ss_lst

    @Attribute
    def add_dihedral(self):
        """
        Returns the updated 'eq_pts_le_ss_lst' and 'eq_pts_te_ss_lst' points defining the leading and trailing edges of
        the sweeped and twisted surface by adding the dihedral distribution to the model. The points are translated
        along the upward z-axis but also along the y-axis to be sure to match the condition in angle at the specific
        location defined by the distribution. The points coordinates are rounded for further processing purposes.

        """
        [eq_pts_le_ss_lst, eq_pts_te_ss_lst] = [self.add_twist[0][:],self.add_twist[1][:]]
        eq_pts_dihedral_distrib_lst = self.calculate_dihedral_distrib[0]
        if not self.is_there_dihedral:
            return eq_pts_le_ss_lst, eq_pts_te_ss_lst

        eq_pts_le_dihedral_ss_lst = eq_pts_le_ss_lst[:]
        eq_pts_te_dihedral_ss_lst = eq_pts_te_ss_lst[:]

        if self.nb_eq_pts_rails == 2:
            for i in range(self.nb_eq_pts_rails):
                eq_pts_le_dihedral_ss_lst[i] = eq_pts_le_ss_lst[i].translate(Vector(0, 0, 1),
                                                                             eq_pts_dihedral_distrib_lst[i][1],
                                                                             Vector(0, 1, 0), -eq_pts_le_ss_lst[i][1] +
                                                                             eq_pts_dihedral_distrib_lst[i][0])
                eq_pts_te_dihedral_ss_lst[i] = eq_pts_te_ss_lst[i].translate(Vector(0, 0, 1),
                                                                             eq_pts_dihedral_distrib_lst[i][1],
                                                                             Vector(0, 1, 0), -eq_pts_le_ss_lst[i][1] +
                                                                             eq_pts_dihedral_distrib_lst[i][0])
        else :
            for i in range(self.nb_eq_pts_rails):
                eq_pts_le_dihedral_ss_lst[i] = eq_pts_le_ss_lst[i].translate(Vector(0, 0, 1),
                                                                         eq_pts_dihedral_distrib_lst[i][1])
                eq_pts_te_dihedral_ss_lst[i] = eq_pts_te_ss_lst[i].translate(Vector(0, 0, 1),
                                                                         eq_pts_dihedral_distrib_lst[i][1])

        return [Point(round(eq_pts_le_dihedral_ss_lst[i][0], 7), round(eq_pts_le_dihedral_ss_lst[i][1], 7),
                      round(eq_pts_le_dihedral_ss_lst[i][2], 7)) for
                i in range(self.nb_eq_pts_rails)], [
                   Point(round(eq_pts_te_dihedral_ss_lst[i][0], 7), round(eq_pts_te_dihedral_ss_lst[i][1], 7),
                         round(eq_pts_te_dihedral_ss_lst[i][2], 7)) for
                   i in range(self.nb_eq_pts_rails)]

    @Attribute
    def determine_position_all_airfoils(self):
        """
        Returns the positions along the y-axis of all the airfoils by considering the defined airfoils but also the
        total number of airfoils given by the user that will serve as cross-sectional curves to fit the final external
        surface of the considered surface.
        A list 'length_defined_airfoils_along_dihedral' is created and contains the distances (lengths) between the
        leading edge root of the surface and the defined airfoils alongside the surface span (not along the y-axis, thus
        take into account the distance 'added' because of dihedral).
        The total number of airfoil is then equally distributed along these distances (considered as sections) to be
        avoid points concentration in the region of low dihedral and points scattering in the high dihedral regions (of
        a winglet for example).
        The dihedral distribution is trimmed into curves defining the sections between the defined airfoils. The y-axis
        coordinates of the airfoils is thus obtained by equally distributing the curve with the number of airfoils
        corresponding to the section.

        """
        crv_dihedral_distrib = self.calculate_dihedral_distrib[1]
        crv_dihedral_distrib = ScaledCurve(crv_dihedral_distrib,Point(0,0,0),self.y_defined_airfoils[-1]/crv_dihedral_distrib.end[0])
        if self.nb_airfoils_ss == self.nb_defined_airfoils:
            return self.y_defined_airfoils
        length_defined_airfoils_along_dihedral = []
        length_defined_airfoils_along_dihedral.append(0)
        crvs_btween_defined_airfoils_lst = []
        length_along_dihedral = 0
        for i in range(1, self.nb_defined_airfoils):
            crv = TrimmedCurve(basis_curve = crv_dihedral_distrib,
                               limit1 = Plane(Point(self.y_defined_airfoils[i - 1], 0, 0), Vector(1, 0, 0)),
                               limit2 = Plane(Point(self.y_defined_airfoils[i], 0, 0), Vector(1, 0, 0)),
                               tolerance = 1e-2)
            crvs_btween_defined_airfoils_lst.append(crv)
            length_along_dihedral = length_along_dihedral + crv.length
            length_defined_airfoils_along_dihedral.append(length_along_dihedral)
        nb_airfoils_sec1 = int(
            (length_defined_airfoils_along_dihedral[1] - length_defined_airfoils_along_dihedral[
                0]) * self.nb_airfoils_ss / (
                    length_defined_airfoils_along_dihedral[-1] - length_defined_airfoils_along_dihedral[0]))
        nb_airfoils_remain = self.nb_airfoils_ss - nb_airfoils_sec1
        y_airfoils_ss = [crvs_btween_defined_airfoils_lst[0].equispaced_points(nb_airfoils_sec1)[k][0] for k in
                         range(nb_airfoils_sec1)]
        # print('sec_1=', y_airfoils_ss)
        for i in range(1, self.nb_defined_airfoils - 2):
            nb_airfoils_sec = int(
                (length_defined_airfoils_along_dihedral[i + 1] - length_defined_airfoils_along_dihedral[
                    i]) * nb_airfoils_remain / (length_defined_airfoils_along_dihedral[-1] -
                                                length_defined_airfoils_along_dihedral[i]))
            # print('nb_airfoils_sec =',(length_defined_airfoils_along_dihedral[i+1]-length_defined_airfoils_along_dihedral[i])*nb_airfoils_remain/(length_defined_airfoils_along_dihedral[-1]-length_defined_airfoils_along_dihedral[i]))
            # print('nb_airf_sec_i=', nb_airfoils_sec)
            nb_airfoils_remain = nb_airfoils_remain - nb_airfoils_sec
            # print('remain_i=',nb_airfoils_remain)
            y_airfoils_ss = y_airfoils_ss[:-1] + [
                crvs_btween_defined_airfoils_lst[i].equispaced_points(nb_airfoils_sec + 1)[k][0] for k in
                range(nb_airfoils_sec + 1)]
            # print('new_sec =',
            #      [crvs_btween_defined_airfoils_lst[i].equispaced_points(nb_airfoils_sec + 1)[k][0] for k in
            #       range(nb_airfoils_sec + 1)])
            # print('y_airfoil_i=', y_airfoils_ss)
        y_airfoils_ss = y_airfoils_ss[:-1] + [
            crvs_btween_defined_airfoils_lst[-1].equispaced_points(nb_airfoils_remain + 1)[k][0] for k in
            range(nb_airfoils_remain + 1)]
        #print('y_airfoil=', y_airfoils_ss)

        return [round(y_airfoils_ss[i], 7) for i in range(self.nb_airfoils_ss)]

    @Attribute
    def determine_position_intermediate_airfoils(self):
        """
        Returns the positions along the y-axis of the intermediate airfoils created between each defined airfoils. It
        basically uses the 'determine_position_all_airfoils' and suppress the y values that correspond to the defined
        airfoils given by the user. This intermediary function is useful to better defined the position of the
        thereafter defined subsections part.

        """
        y_subsections_lst = []
        j_lst = []
        y_airfoil_ss = self.determine_position_all_airfoils
        y_defined_airfoils = [float(self.y_defined_airfoils[i]) for i in range(self.nb_defined_airfoils)]
        # print(y_defined_airfoils)
        j = 1
        i = 1
        while j < self.nb_defined_airfoils:
            while y_airfoil_ss[i] < y_defined_airfoils[j] - 1e-5 or y_airfoil_ss[i] > y_defined_airfoils[j] + 1e-5:
                y_subsections_lst.append(y_airfoil_ss[i])
                j_lst.append(j - 1)
                i += 1
            j += 1
            i += 1
        # print("y_subsections_lst=", y_subsections_lst)
        # print("j_lst=", j_lst)
        return y_subsections_lst, j_lst

    @Part
    def sections(self):
        """"
        Returns the cross-sectional sections corresponding to the defined airfoils given by the user. For further
        information, see 'Section' subclass.

        """
        return Section(quantify = self.nb_defined_airfoils,
                       airfoil_name = self.airfoils_names[child.index],
                       y_airfoil = self.y_defined_airfoils[child.index],
                       cant_angle = self.cant_angles_airfoils[child.index],
                       dihedral_angle = self.dihedral_angles_airfoils[child.index],
                       eq_pts_le_ss_lst = self.add_dihedral[0],
                       eq_pts_te_ss_lst = self.add_dihedral[1],
                       height_factor = self.heights_factors[child.index],
                       nb_eq_pts_rails = self.nb_eq_pts_rails,
                       nb_eq_pts_airfoils = self.nb_eq_pts_airfoils)

    @Part
    def subsections(self):
        """
        Returns the cross-sectional subsections corresponding to the intermediate airfoils (between the defined airfoil
        given by the user created to fit the final external surface. For further information, see 'SubSection' subclass.

        """
        return SubSection(quantify = self.nb_airfoils_ss - self.nb_defined_airfoils,
                          y_airfoil = self.determine_position_intermediate_airfoils[0][child.index],
                          nb_eq_pts_rails = self.nb_eq_pts_rails,
                          nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                          eq_pts_le_ss_lst = self.add_dihedral[0],
                          eq_pts_te_ss_lst = self.add_dihedral[1],
                          crv_defined_airfoils_ss_1 =
                          self.sections[self.determine_position_intermediate_airfoils[1][child.index]].curve[0],
                          crv_defined_airfoils_ss_2 =
                          self.sections[self.determine_position_intermediate_airfoils[1][child.index] + 1].curve[0],
                          y_defined_airfoils = self.y_defined_airfoils,
                          cant_angles_airfoils = self.cant_angles_airfoils,
                          dihedral_angles_airfoils = self.dihedral_angles_airfoils,
                          chords_defined_airfoils = self.chords_defined_airfoils,
                          j = self.determine_position_intermediate_airfoils[1][child.index])

    @Attribute
    def chords_defined_airfoils(self):
        """
        Returns the chords of the defined airfoils given by the user and defining the 'sections' part.

        """
        chords_defined_airfoils_lst = []
        for i in range(self.nb_defined_airfoils):
            chords_defined_airfoils_lst.append(self.sections[i].curve[1])
        # print("chords_defined_airfoils_lst=",chords_defined_airfoils_lst)
        return chords_defined_airfoils_lst

    @Attribute
    def crvs_surface_defined_airfoils(self):
        """
        Returns the chords of the intermediate airfoils defining the 'subsections' part.

        """
        if self.nb_defined_airfoils == 2:
            [eq_pts_le_ss_lst, eq_pts_te_ss_lst] = self.add_dihedral
            crvs_airfoils_vb_lst = []
            crvs_airfoils_vb_lst.append(
                Airfoil(airfoil_name = self.airfoils_names[0], pt_le = eq_pts_le_ss_lst[0], pt_te = eq_pts_te_ss_lst[0],
                        dihedral_angle = self.dihedral_angles_airfoils[0], height_factor = self.heights_factors[0],
                        nb_eq_pts_airfoils = self.nb_eq_pts_airfoils).airfoil)
            crvs_airfoils_vb_lst.append(
                Airfoil(airfoil_name = self.airfoils_names[1], pt_le = eq_pts_le_ss_lst[1], pt_te = eq_pts_te_ss_lst[1],
                        dihedral_angle = self.dihedral_angles_airfoils[1], height_factor = self.heights_factors[1],
                        nb_eq_pts_airfoils = self.nb_eq_pts_airfoils).airfoil)

            return crvs_airfoils_vb_lst
        crvs_surface_lst = []
        for j in range(self.nb_defined_airfoils):
            crvs_surface_lst.append(self.sections[j].curve)
        return crvs_surface_lst

    @Attribute
    def chords_all_airfoils(self):
        """
        Returns the chords of all the airfoils along the surface span. Useful to evaluate the total area of the
        planform.

        """
        if self.nb_airfoils_ss == self.nb_defined_airfoils:
            return self.chords_defined_airfoils
        tol = 1e-5
        y_defined_airfoils = self.y_defined_airfoils
        y_subsections_lst = self.determine_position_intermediate_airfoils[0]
        chords_lst = []
        i = 0
        j = 0
        while j < self.nb_defined_airfoils - 1:
            chords_lst.append(self.sections[j].curve[1])
            while i < self.nb_airfoils_ss - self.nb_defined_airfoils and y_defined_airfoils[j] + tol < \
                    y_subsections_lst[i] < y_defined_airfoils[
                j + 1] - tol:
                chords_lst.append(self.subsections[i].curve[3])
                i += 1
            j += 1
        chords_lst.append(self.sections[-1].curve[1])
        #print('chords_lst=', chords_lst)
        return chords_lst

    @Attribute
    def crvs_surface_all_airfoils(self):
        """
        Returns the curves of all the airfoils distributed along the surface.
        For j varying from 0 to the number of defined airfoils-1, the section j curve is added to the list. Then, all the
        intermediate curves between section j and section j+1 are added (the 'subsections' curves). A tolerance is added
        to avoid the rounding error to interfere with the correct detection of the sections curves. When all the
        intermediate airfoils (between section j and j+1) have been added, the section j+1 curve is in turn added.

        """
        if self.nb_airfoils_ss == self.nb_defined_airfoils:
            return self.crvs_surface_defined_airfoils
        crvs_airfoils_surface_lst = []
        tol = 1e-5
        y_defined_airfoils = self.y_defined_airfoils
        y_subsections_lst = self.determine_position_intermediate_airfoils[0]
        crvs_surface_lst = []
        pts = []
        # crvs_airfoils_lst = []
        i = 0
        j = 0
        while j < self.nb_defined_airfoils - 1:
            crvs_surface_lst.append(self.sections[j].curve[0])
            while i < self.nb_airfoils_ss - self.nb_defined_airfoils and y_defined_airfoils[j] + tol < \
                    y_subsections_lst[i] < y_defined_airfoils[
                j + 1] - tol:
                crvs_surface_lst.append(self.subsections[i].curve[0])
                crvs_airfoils_surface_lst.append(self.subsections[i].curve[1])
                pts.append(self.subsections[i].curve[2])
                # pts.append(self.subsections[i].curve[1])
                # pts.append(self.subsections[i].curve[2])
                # crvs_airfoils_lst.append(self.subsections[i].curve[3])
                i += 1
            j += 1
            # crvs_airfoils_lst.append(self.subsections[i-1].curve[3])
        crvs_surface_lst.append(self.sections[-1].curve[0])
        return crvs_surface_lst


    @Attribute
    def crvs_surface_to_points(self):
        """
        Returns the points equally distributed along the curves of all the airfoils distributed along the surface. This
        function is useful to defined the mirrored surface of the surface by using the 'mirror_points' function.

        """
        crvs_surface_lst = self.crvs_surface_all_airfoils
        pts_surface_lst = []
        for i in range(self.nb_airfoils_ss):
            pts_airfoil_lst = crvs_surface_lst[i].equispaced_points(self.nb_eq_pts_airfoils)
            pts_surface_lst.append(pts_airfoil_lst)
        return pts_surface_lst

    @Attribute
    def mirror_points(self):
        """
        Returns the mirrored points of all the airfoils curves distributed along the surface. The mirror transformation
        is made with respect to the plane (OXY).

        """
        pts_surface_ss_lst = self.crvs_surface_to_points
        pts_surface_ss_mirored_lst = []
        for i in range(len(pts_surface_ss_lst)):
            pts_airfoil_mirored_lst = []
            for k in range(self.nb_eq_pts_airfoils):
                pt_mirored = Point(pts_surface_ss_lst[i][k][0], -pts_surface_ss_lst[i][k][1],
                                   pts_surface_ss_lst[i][k][2])
                pts_airfoil_mirored_lst.append(pt_mirored)
            pts_surface_ss_mirored_lst.append(pts_airfoil_mirored_lst)
        return pts_surface_ss_mirored_lst

    @Part
    def surface_ss(self):
        """
        Returns the external surface of the surface by fitting all the points defining all the airfoils distributed
        along the surface. A mininmum degree of 8 is applied to closely fit the points, which is proving useful when
        high dihedral changes characterize the surface.

        """
        return FittedSurface(points = self.crvs_surface_to_points, min_degree = 8)

    @Part
    def mirorred_surface(self):
        """
        Returns the mirrored external surface corresponding to the points of the airfoils of the surface previously
        mirrored with the 'mirror_points' function.

        """
        return FittedSurface(points = self.mirror_points, min_degree = 8)

    @Part
    def surface(self):
        """
        Returns the surface of the surface and of its mirror surface with respect to the (OXY) plane. A better approach
        would have been to directly fit all along the span of the two joined surfaces the external surface. The
        difficulties experienced by the fit function to obtain a feasible result (see commented @Part surface below) has
        prompted to directly fused two fitted surfaces at their roots even thought a discontinuity can appear at this
        location.

        """
        return Fused(self.surface_ss, FittedSurface(points = self.mirror_points, min_degree = 8))

    # @Part
    # def surface(self):
    #    return FittedSurface(points=self.crvs_surface_to_points+self.mirror_points, min_degree = 8)

    @Attribute
    def determine_position_sections_avl(self):
        """
        Returns the y-positions of the AVL sections. These are equally distributed along the span of the surface given
        a number of sections chosen to conduct the AVL analyse.

        """
        y_sections_AVL = [0]
        for i in range(1,self.nb_sections_AVL):
            y_sections_AVL.append(self.sspan/self.nb_sections_AVL*(i+1))
        return y_sections_AVL

    @Attribute
    def crvs_sections_AVL(self):
        """
        Returns the curves of the AVL sections.
        For each y-position of an AVL section, the nearest airfoil curve among all the airfoils curves distributed
        along the surface is assigned to be the AVL section curve.

        """
        crvs_AVL_lst = []
        crvs_surface_lst = self.crvs_surface_all_airfoils
        y_airfoil_ss_array = np.asarray(self.determine_position_all_airfoils)
        y_airfoils_AVL = self.determine_position_sections_avl
        for i in range(self.nb_sections_AVL):
            idx = (np.abs(y_airfoil_ss_array - y_airfoils_AVL[i])).argmin()
            crvs_AVL_lst.append(avl.SectionFromCurve(curve_in =crvs_surface_lst[idx]))
        print("Surface '" + self.surface_type + "' placed")
        return crvs_AVL_lst

    @Part
    def avl_surface(self):
        """
        Returns the AVL surface to be analysed.

        """
        return avl.Surface(name = self.name,
                           n_chordwise = self.nb_chordwise_vortices,
                           chord_spacing = avl.Spacing.cosine,
                           n_spanwise = self.nb_spanwise_vortices,
                           span_spacing = avl.Spacing.cosine,
                           y_duplicate = self.position.point[1] if self.is_mirrored else None,
                           sections = self.crvs_sections_AVL)


if __name__ == '__main__':
    from parapy.gui import display

    obj = Surface_BWB()
    display(obj)
