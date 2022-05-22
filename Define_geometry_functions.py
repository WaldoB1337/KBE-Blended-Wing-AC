from copy import deepcopy
from math import *
import numpy as np
from parapy.gui import display
from parapy.geom import *
from parapy.core import *
from parapy.core.validate import AdaptedValidator

import kbeutils.avl as avl
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve


# Airfoil definition
def airfoil(airfoil_name, pt_le, pt_te, dihedral_angle):
    # airfoil_name = "whitcomb"
    # pt_le = Point(0,0,0)
    # pt_te = Point(1,0,2)
    chord = -pt_le.distance(pt_te)

    if len(airfoil_name) == 4 and (type(int(char) == int) for char in airfoil_name):
        print('NACA4')
        crv_airfoil = DynamicType(type = Naca4AirfoilCurve,
                                  designation = airfoil_name,
                                  mesh_deflection = 0.00001,
                                  hidden = True)

        crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = OXY, to_position = Position(pt_le))
        crv_airfoil = ScaledCurve(crv_airfoil,
                                  pt_le,
                                  chord,
                                  mesh_deflection = 0.00001)

    if len(airfoil_name) == 5 and (type(int(char) == int) for char in airfoil_name):
        print('NACA5')
        crv_airfoil = DynamicType(type = Naca5AirfoilCurve,
                                  designation = airfoil_name,
                                  mesh_deflection = 0.00001,
                                  hidden = True)

        crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = OXY, to_position = Position(pt_le))
        crv_airfoil = ScaledCurve(crv_airfoil,
                                  pt_le,
                                  chord,
                                  mesh_deflection = 0.00001
                                  )

    if not ((len(airfoil_name) == 4 or len(airfoil_name) == 5) and (type(int(char) == int) for char in airfoil_name)):
        print('Airfoil specified')
        with open(airfoil_name + ".dat", 'r') as f:
            pts_airfoil_lst = []
            for line in f:
                x, z = line.split(' ', 1)  # the cartesian coordinates are directly interpreted as X and Z coordinates
                pts_airfoil_lst.append(pt_le.translate(
                    "x", float(x) * chord,  # the x points are scaled according to the airfoil chord length
                    "z", float(z) * chord))  # the y points are scaled according to the
        # Creation of a curve te->le(up)+le->te(lo)

        idx = pts_airfoil_lst.index(pt_le)
        upper_pts_lst = pts_airfoil_lst[:idx + 1]
        lower_pts_lst = pts_airfoil_lst[idx:]

        crv_upper = FittedCurve(points = upper_pts_lst, max_degree = 8)
        crv_lower = FittedCurve(points = lower_pts_lst, max_degree = 8)
        crv_airfoil = ComposedCurve(built_from = [crv_upper, crv_lower], allow_multiple = 'True')

    rot_angle_z = np.arcsin((pt_te[1] - pt_le[1]) / chord)
    crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                               vector = Vector(0, 0, 1), angle = rot_angle_z)
    rot_angle_y = np.arcsin((pt_te[2] - pt_le[2]) / chord)
    crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                               vector = Vector(0, 1, 0), angle = -rot_angle_y)
    crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                               vector = Vector(pt_te[0] - pt_le[0], pt_te[1] - pt_le[1], pt_te[2] - pt_le[2]),
                               angle = -np.deg2rad(dihedral_angle))

    # display(crv_airfoil)
    return crv_airfoil


def sweep_crvs(nb_eq_pts_rails, nb_fixed_control_points, y_fcp, chords, sweeps_le, sweeps_te=[], rho_le=[], rho_te=[],
               surface_type='vb', x_fcp=[]):
    # Creation of the sweep curves
    lambda_le = np.deg2rad(sweeps_le)
    if surface_type == 'vb':
        x_fcp = [0, -tan(lambda_le[0]) * span_vb]
        crv_le = LineSegment(start = Point(x_fcp_vb[0], y_fcp_vb[0], 0), end = Point(x_fcp_vb[1], y_fcp_vb[1], 0))
        crv_te = LineSegment(start = Point(x_fcp_vb[0] - chords_vb[0], y_fcp_vb[0], 0),
                             end = Point(x_fcp_vb[1] - chords_vb[1], y_fcp_vb[1], 0))

    else:
        lambda_te = np.deg2rad(sweeps_te)
        crvs_le_lst = []
        crvs_te_lst = []

        for i in range(nb_fixed_control_points - 1):
            xc_1_le = x_fcp[i] + rho_le[i] * -sin(lambda_le[i])
            yc_1_le = y_fcp[i] + rho_le[i] * cos(lambda_le[i])
            xc_1_te = x_fcp[i] - chords[i] + rho_te[i] * -sin(lambda_te[i])
            yc_1_te = y_fcp[i] + rho_te[i] * cos(lambda_te[i])

            xc_2_le = x_fcp[i + 1] + rho_le[i + 1] * -sin(lambda_le[i + 1] + pi)
            yc_2_le = y_fcp[i + 1] + rho_le[i + 1] * cos(lambda_le[i + 1] + pi)
            xc_2_te = x_fcp[i + 1] - chords[i + 1] + rho_te[i + 1] * -sin(lambda_te[i + 1] + pi)
            yc_2_te = y_fcp[i + 1] + rho_te[i + 1] * cos(lambda_te[i + 1] + pi)

            FCP0 = Point(x_fcp[i], y_fcp[i], 0)
            VCP0 = Point(xc_1_le, yc_1_le, 0)
            VCP1 = Point(xc_2_le, yc_2_le, 0)
            FCP1 = Point(x_fcp[i + 1], y_fcp[i + 1], 0)
            pts = [FCP0, VCP0, VCP1, FCP1]
            crv = BezierCurve(pts)

            crvs_le_lst.append(crv)

            FCP0 = Point(x_fcp[i] - chords[i], y_fcp[i], 0)
            VCP0 = Point(xc_1_te, yc_1_te, 0)
            VCP1 = Point(xc_2_te, yc_2_te, 0)
            FCP1 = Point(x_fcp[i + 1] - chords[i + 1], y_fcp[i + 1], 0)
            pts = [FCP0, VCP0, VCP1, FCP1]
            crv = BezierCurve(pts)

            crvs_te_lst.append(crv)

        crv_le = ComposedCurve(built_from = crvs_le_lst, allow_multiple = 'True')
        # display(crv_le)
        crv_te = ComposedCurve(built_from = crvs_te_lst, allow_multiple = 'True')
        # display(crv_te)

    eq_ss_lst = np.linspace(0, y_fcp[-1], nb_eq_pts_rails)
    eq_pts_rails_lst = [Point(0, eq_ss_lst[i], 0) for i in range(nb_eq_pts_rails)]
    eq_pts_le_ss_lst = []
    eq_pts_te_ss_lst = []
    for i in range(nb_eq_pts_rails):
        pt_le = crv_le.intersection_point(Plane(reference = eq_pts_rails_lst[i], normal = Vector(0, 1, 0)))
        eq_pts_le_ss_lst.append(pt_le)
        pt_te = crv_te.intersection_point(Plane(reference = eq_pts_rails_lst[i], normal = Vector(0, 1, 0)))
        eq_pts_te_ss_lst.append(pt_te)

    return eq_pts_le_ss_lst, eq_pts_te_ss_lst


def calculate_twist_distrib(nb_eq_pts_rails, nb_pts_twist=0, angles_twist=[], y_twist=[], nb_pts_twist_axis=0,
                            x_twist_axis=[], y_twist_axis=[], is_twisted=False):
    if not is_twisted:
        return []
    # Twist distribution
    pts_twist_lst = []
    for i in range(nb_pts_twist):
        pt = Point(y_twist[i], angles_twist[i], 0)
        pts_twist_lst.append(pt)
    crv_twist_distrib = FittedCurve(points = pts_twist_lst)

    # Twist axis definition
    crvs_twist_axis = []
    for i in range(nb_pts_twist_axis - 1):
        crv = LineSegment(Point(x_twist_axis[i], y_twist_axis[i], 0),
                          Point(x_twist_axis[i + 1], y_twist_axis[i + 1], 0))
        crvs_twist_axis.append(crv)

    crv_twist_axis = ComposedCurve(built_from = crvs_twist_axis)

    eq_ss_lst = np.linspace(0, y_twist[-1], nb_eq_pts_rails)
    eq_pts_twist_distrib = [Point(eq_ss_lst[i], 0, 0) for i in range(nb_eq_pts_rails)]
    eq_pts_twist_axis = [Point(0, eq_ss_lst[i], 0) for i in range(nb_eq_pts_rails)]

    eq_pts_twist_distrib_lst = []
    eq_pts_twist_axis_lst = []

    for i in range(nb_eq_pts_rails):
        pt_twist_distrib = crv_twist_distrib.intersection_point(
            Plane(reference = eq_pts_twist_distrib[i], normal = Vector(1, 0, 0)))
        eq_pts_twist_distrib_lst.append(pt_twist_distrib)
        pt_twist_axis = crv_twist_axis.intersection_point(
            Plane(reference = eq_pts_twist_axis[i], normal = Vector(0, 1, 0)))
        eq_pts_twist_axis_lst.append(pt_twist_axis)

    return eq_pts_twist_distrib_lst, eq_pts_twist_axis_lst


def calculate_dihedral_distrib(nb_eq_pts_rails, nb_pts_dihedral, z_dihedral, y_dihedral, weights_dihedral=[],
                               is_dihedral=True):
    if not is_dihedral:
        return []
    if nb_pts_dihedral == 2:
        crv_dihedral_distrib = FittedCurve(
            [Point(y_dihedral[0], z_dihedral[0], 0), Point(y_dihedral[1], z_dihedral[1], 0)])

    else:
        pts_dihedral_lst = []
        FCPdi = Point(y_dihedral[0], z_dihedral[0], 0)
        VCPdi = Point(((y_dihedral[1] + y_dihedral[0]) / 2), (z_dihedral[1] + z_dihedral[0]) / 2, 0)
        pts_dihedral_lst.append(FCPdi)
        pts_dihedral_lst.append(VCPdi)
        for i in range(1, nb_pts_dihedral - 1):
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
                                           weights = weights_dihedral)
        # display(crv_dihedral_distrib)

    eq_ss_lst = np.linspace(0, y_dihedral[-1], nb_eq_pts_rails)
    eq_pts_dihedral_distrib = [Point(eq_ss_lst[i], 0, 0) for i in range(nb_eq_pts_rails)]
    eq_pts_dihedral_distrib_lst = []

    for i in range(nb_eq_pts_rails):
        pt_dihedral = crv_dihedral_distrib.intersection_point(
            Plane(reference = eq_pts_dihedral_distrib[i], normal = Vector(1, 0, 0)))
        eq_pts_dihedral_distrib_lst.append(pt_dihedral)

    return eq_pts_dihedral_distrib_lst, crv_dihedral_distrib


def add_twist(nb_eq_pts_rails, eq_pts_le_ss_lst, eq_pts_te_ss_lst, eq_pts_twist_distrib_lst=[],
              eq_pts_twist_axis_lst=[], is_twisted=False):
    if not is_twisted:
        return eq_pts_le_ss_lst, eq_pts_te_ss_lst

    for i in range(nb_eq_pts_rails):
        height_added_le = (eq_pts_le_ss_lst[i][0] - eq_pts_twist_axis_lst[i][0]) * tan(
            np.deg2rad(eq_pts_twist_distrib_lst[i][1]))
        height_added_te = (eq_pts_te_ss_lst[i][0] - eq_pts_twist_axis_lst[i][0]) * tan(
            np.deg2rad(eq_pts_twist_distrib_lst[i][1]))

        eq_pts_le_ss_lst[i] = eq_pts_le_ss_lst[i].translate(Vector(0, 0, 1), height_added_le)
        eq_pts_te_ss_lst[i] = eq_pts_te_ss_lst[i].translate(Vector(0, 0, 1), - height_added_te)

    return eq_pts_le_ss_lst, eq_pts_te_ss_lst


def add_dihedral(nb_eq_pts_rails, eq_pts_le_ss_lst, eq_pts_te_ss_lst, eq_pts_dihedral_distrib_lst=[],
                 is_dihedral=False):
    if not is_dihedral:
        return eq_pts_le_ss_lst, eq_pts_te_ss_lst

    eq_pts_le_dihedral_ss_lst = deepcopy(eq_pts_le_ss_lst)
    eq_pts_te_dihedral_ss_lst = deepcopy(eq_pts_te_ss_lst)
    for i in range(nb_eq_pts_rails):
        eq_pts_le_dihedral_ss_lst[i] = eq_pts_le_ss_lst[i].translate(Vector(0, 0, 1), eq_pts_dihedral_distrib_lst[i][1],
                                                                     Vector(0, 1, 0), -eq_pts_le_ss_lst[i][1] +
                                                                     eq_pts_dihedral_distrib_lst[i][0])
        eq_pts_te_dihedral_ss_lst[i] = eq_pts_te_ss_lst[i].translate(Vector(0, 0, 1), eq_pts_dihedral_distrib_lst[i][1],
                                                                     Vector(0, 1, 0), -eq_pts_le_ss_lst[i][1] +
                                                                     eq_pts_dihedral_distrib_lst[i][0])

    return eq_pts_le_dihedral_ss_lst, eq_pts_te_dihedral_ss_lst


# Mirror the points

# def mirror_rails(eq_pts_le_ss_lst,eq_pts_te_ss_lst):
#    eq_pts_le_ss_mirored_lst = []
#    eq_pts_te_ss_mirored_lst = []

#    for i in range(nb_eq_pts_rails):
#        eq_pts_le_ss_mirored_lst.append(Point(eq_pts_le_ss_lst[i][0], -eq_pts_le_ss_lst[i][1], eq_pts_le_ss_lst[i][2]))
#        eq_pts_te_ss_mirored_lst.append(Point(eq_pts_te_ss_lst[i][0], -eq_pts_te_ss_lst[i][1], eq_pts_te_ss_lst[i][2]))

#    return eq_pts_le_ss_mirored_lst,eq_pts_te_ss_mirored_lst

## End of the definition of rails

# display(eq_pts_le_ss_lst + eq_pts_le_ss_mirored_lst + eq_pts_te_ss_lst + eq_pts_te_ss_mirored_lst)

def place_defined_airfoils(nb_eq_pts_rails, nb_defined_airfoils, airfoils_names, y_defined_airfoils,
                           dihedral_angles_airfoils, eq_pts_le_ss_lst, eq_pts_te_ss_lst, cant_angles_airfoils=[]):
    if nb_defined_airfoils == 2:
        crvs_airfoils_vb_lst = []
        crvs_airfoils_vb_lst.append(
            airfoil(airfoils_names[0], eq_pts_le_ss_lst[0], eq_pts_te_ss_lst[0], dihedral_angles_airfoils[0]))
        crvs_airfoils_vb_lst.append(
            airfoil(airfoils_names[1], eq_pts_le_ss_lst[-1], eq_pts_te_ss_lst[-1], dihedral_angles_airfoils[-1]))

        return crvs_airfoils_vb_lst

    crvs_defined_airfoils_ss = []

    # Define le and te of each airfoil and place it in space
    for i in range(nb_defined_airfoils):

        pln_le = Plane(reference = Point(0, y_defined_airfoils[i], 0), normal = Vector(0, 1, 0))
        p = 0
        while (p < nb_eq_pts_rails - 1) and not (
                eq_pts_le_ss_lst[p][1] <= y_defined_airfoils[i] <= eq_pts_le_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_le = LineSegment(start = eq_pts_le_ss_lst[p], end = eq_pts_le_ss_lst[p + 1]).surface_intersections(
            pln_le)
        pt_le = lst_inter_le[0]['point']
        p = 0
        while (p < nb_eq_pts_rails - 1) and not (
                eq_pts_te_ss_lst[p][1] <= y_defined_airfoils[i] <= eq_pts_te_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_te_no_cant = LineSegment(start = eq_pts_te_ss_lst[p],
                                           end = eq_pts_te_ss_lst[p + 1]).surface_intersections(pln_le)
        pt_te_no_cant = lst_inter_te_no_cant[0]['point']
        pt_te_projected = pt_te_no_cant.translate(Vector(0, 1, 0),
                                                  np.sin(np.deg2rad(cant_angles_airfoils[i])) * pt_le.distance(
                                                      pt_te_no_cant))

        pln_te = Plane(reference = pt_te_projected,
                       normal = Vector(pt_le[2] - pt_te_projected[2] - (pt_le[1] - pt_te_projected[1]),
                                       pt_le[0] - pt_te_projected[0], 0))

        p = 0
        while (p < nb_eq_pts_rails - 1) and not (
                eq_pts_te_ss_lst[p][1] <= pt_te_projected[1] <= eq_pts_te_ss_lst[p + 1][1]):
            p = p + 1
        lst_inter_te = LineSegment(start = eq_pts_te_ss_lst[p], end = eq_pts_te_ss_lst[p + 1]).surface_intersections(
            pln_te)
        pt_te = lst_inter_te[0]['point']

        crvs_defined_airfoils_ss.append(airfoil(airfoils_names[i], pt_le, pt_te, dihedral_angles_airfoils[i]))

    return crvs_defined_airfoils_ss


# display(crv_te)
# display([crv_le,crv_te,crvs_airfoils_ss])
# display([crv_le_ss,crv_te_ss,pts])

# Mirror the airfoils

# def mirror_defined_airfoils(crvs_airfoils_ss):
#    crvs_defined_airfoils_ss_mirored = []
#    for i in range(len(crvs_airfoils_ss)):
#        crvs_defined_airfoils_ss_mirored.append(MirroredCurve(curve_in = crvs_airfoils_ss[i], reference_point = Point(0, 0, 0),
#                                                      vector1 = Vector(1, 0, 0),
#                                                      vector2 = Vector(0, 0, 1)))
#    return crvs_defined_airfoils_ss_mirored

# crvs_airfoils = crvs_airfoils_ss + crvs_airfoils_ss_mirored
# display([eq_pts_le_ss_lst,crvs_airfoils])

# display([crv_le,crv_te]+crvs_airfoils+pts)

# Define intermediate airfoils location


def determine_position_intermediate_airfoils(nb_airfoils_ss, nb_defined_airfoils, y_defined_airfoils,
                                             crv_dihedral_distrib):
    if nb_airfoils_ss == nb_defined_airfoils:
        return y_defined_airfoils
    length_defined_airfoils_along_dihedral = []
    length_defined_airfoils_along_dihedral.append(0)
    crvs_btween_defined_airfoils_lst = []
    length_along_dihedral = 0
    for i in range(1, nb_defined_airfoils):
        print('trim_nb=', i)
        crv = TrimmedCurve(basis_curve = crv_dihedral_distrib,
                           limit1 = Plane(Point(y_defined_airfoils[i - 1], 0, 0), Vector(1, 0, 0)),
                           limit2 = Plane(Point(y_defined_airfoils[i], 0, 0), Vector(1, 0, 0)))
        crvs_btween_defined_airfoils_lst.append(crv)
        length_along_dihedral = length_along_dihedral + crv.length
        length_defined_airfoils_along_dihedral.append(length_along_dihedral)
    print('length_defined_airfoils_along_dihedral=', length_defined_airfoils_along_dihedral)
    # display(crvs_btween_defined_airfoils_lst)

    nb_airfoils_sec1 = int(
        (length_defined_airfoils_along_dihedral[1] - length_defined_airfoils_along_dihedral[0]) * nb_airfoils_ss / (
                length_defined_airfoils_along_dihedral[-1] - length_defined_airfoils_along_dihedral[0]))
    nb_airfoils_remain = nb_airfoils_ss - nb_airfoils_sec1
    y_airfoils_ss = [crvs_btween_defined_airfoils_lst[0].equispaced_points(nb_airfoils_sec1)[k][0] for k in
                     range(nb_airfoils_sec1)]
    print('sec_1=', y_airfoils_ss)
    for i in range(1, nb_defined_airfoils - 2):
        nb_airfoils_sec = int((length_defined_airfoils_along_dihedral[i + 1] - length_defined_airfoils_along_dihedral[
            i]) * nb_airfoils_remain / (length_defined_airfoils_along_dihedral[-1] -
                                        length_defined_airfoils_along_dihedral[i]))
        # print('nb_airfoils_sec =',(length_defined_airfoils_along_dihedral[i+1]-length_defined_airfoils_along_dihedral[i])*nb_airfoils_remain/(length_defined_airfoils_along_dihedral[-1]-length_defined_airfoils_along_dihedral[i]))
        print('nb_airf_sec_i=', nb_airfoils_sec)
        nb_airfoils_remain = nb_airfoils_remain - nb_airfoils_sec
        # print('remain_i=',nb_airfoils_remain)
        y_airfoils_ss = y_airfoils_ss[:-1] + [
            crvs_btween_defined_airfoils_lst[i].equispaced_points(nb_airfoils_sec + 1)[k][0] for k in
            range(nb_airfoils_sec + 1)]
        print('new_sec =', [crvs_btween_defined_airfoils_lst[i].equispaced_points(nb_airfoils_sec + 1)[k][0] for k in
                            range(nb_airfoils_sec + 1)])
        print('y_airfoil_i=', y_airfoils_ss)
    y_airfoils_ss = y_airfoils_ss[:-1] + [
        crvs_btween_defined_airfoils_lst[-1].equispaced_points(nb_airfoils_remain + 1)[k][0] for k in
        range(nb_airfoils_remain + 1)]
    print('y_airfoil=', y_airfoils_ss)

    return y_airfoils_ss


def place_intermediate_airfoils(object_type, crvs_defined_airfoils_ss, nb_eq_pts_rails=0, nb_eq_pts_airfoils=0,
                                nb_defined_airfoils=0, y_defined_airfoils=[],
                                cant_angles_airfoils=[], dihedral_angles_airfoils=[], y_airfoils_ss=[],
                                eq_pts_le_ss_lst=[], eq_pts_te_ss_lst=[], tol=0, round_nb=0):
    # better : if y_airfoils_ss == y_defined_airfoils, return crvs.defi...
    if object_type == 'vb':
        return crvs_defined_airfoils_ss

    y_defined_airfoils = [round(y_defined_airfoils[i], round_nb) for i in range(nb_defined_airfoils)]
    eq_pts_te_ss_lst = [
        Point(round(eq_pts_te_ss_lst[i][0], round_nb), round(eq_pts_te_ss_lst[i][1], round_nb),
              round(eq_pts_te_ss_lst[i][2], round_nb)) for
        i in range(nb_eq_pts_rails)]
    eq_pts_le_ss_lst = [
        Point(round(eq_pts_le_ss_lst[i][0], round_nb), round(eq_pts_le_ss_lst[i][1], round_nb),
              round(eq_pts_le_ss_lst[i][2], round_nb)) for
        i in range(nb_eq_pts_rails)]

    crvs_surface_lst = []
    i = 1
    j = 0
    while j < nb_defined_airfoils - 1:
        print('j =', j)
        pts_airfoil_1_lst = crvs_defined_airfoils_ss[j].equispaced_points(nb_eq_pts_airfoils)
        pts_airfoil_2_lst = crvs_defined_airfoils_ss[j + 1].equispaced_points(nb_eq_pts_airfoils)

        crvs_surface_lst.append(crvs_defined_airfoils_ss[j])
        print('i =', i)

        while i < nb_airfoils_ss and y_defined_airfoils[j] + tol < y_airfoils_ss[i] < y_defined_airfoils[j + 1] - tol:
            pts_airfoil_lst = []
            cant_angle = cant_angles_airfoils[j] + ((y_airfoils_ss[i] - y_defined_airfoils[j]) / (
                    y_defined_airfoils[j + 1] - y_defined_airfoils[j])) * (
                                 cant_angles_airfoils[j + 1] - cant_angles_airfoils[j])
            print('cant_angle' + str(i) + '=', cant_angle)
            dihedral_angle = dihedral_angles_airfoils[j] + ((y_airfoils_ss[i] - y_defined_airfoils[j]) / (
                    y_defined_airfoils[j + 1] - y_defined_airfoils[j])) * (
                                     dihedral_angles_airfoils[j + 1] - dihedral_angles_airfoils[j])

            pln_le = Plane(reference = Point(0, y_airfoils_ss[i], 0), normal = Vector(0, 1, 0))
            p = 0
            # print("eq_pts_le_ss_lst=",eq_pts_le_ss_lst)
            while (p < nb_eq_pts_rails - 1) and not (
                    eq_pts_le_ss_lst[p][1] <= y_airfoils_ss[i] <= eq_pts_le_ss_lst[p + 1][1]):
                p = p + 1
            lst_inter_le = LineSegment(start = eq_pts_le_ss_lst[p],
                                       end = eq_pts_le_ss_lst[p + 1]).surface_intersections(
                pln_le)
            pt_le = lst_inter_le[0]['point']

            p = 0
            while (p < nb_eq_pts_rails - 1) and not (
                    eq_pts_te_ss_lst[p][1] <= y_airfoils_ss[i] <= eq_pts_te_ss_lst[p + 1][1]):
                p = p + 1
            lst_inter_te_no_cant = LineSegment(start = eq_pts_te_ss_lst[p],
                                               end = eq_pts_te_ss_lst[p + 1]).surface_intersections(pln_le)
            pt_te_no_cant = lst_inter_te_no_cant[0]['point']
            pt_te_projected = pt_te_no_cant.translate(Vector(0, 1, 0),
                                                      np.sin(np.deg2rad(cant_angle)) * pt_le.distance(
                                                          pt_te_no_cant))
            pln_te = Plane(reference = pt_te_projected,
                           normal = Vector(pt_le[2] - pt_te_projected[2] - (pt_le[1] - pt_te_projected[1]),
                                           pt_le[0] - pt_te_projected[0], 0))

            p = 0
            while (p < nb_eq_pts_rails - 1) and not (
                    eq_pts_te_ss_lst[p][1] <= pt_te_projected[1] <= eq_pts_te_ss_lst[p + 1][1]):
                p = p + 1

            lst_inter_te = LineSegment(start = eq_pts_te_ss_lst[p],
                                       end = eq_pts_te_ss_lst[p + 1]).surface_intersections(
                pln_te)
            if not lst_inter_te == []:
                pt_te = lst_inter_te[0]['point']
            else:
                print('Intersection with plane not found')
                display([LineSegment(start = eq_pts_te_ss_lst[p - 1], end = eq_pts_te_ss_lst[p]),
                         LineSegment(start = eq_pts_te_ss_lst[p], end = eq_pts_te_ss_lst[p + 1]),
                         LineSegment(start = eq_pts_te_ss_lst[p + 1], end = eq_pts_te_ss_lst[p + 2]), pln_te])
                break
            # display([crv_le_ss,crv_te_ss,pts])

            chord = -pt_le.distance(pt_te)

            chord_interp = -pt_le.distance(pt_te_no_cant)

            # te->le(up)+le->te(lo) for x
            x_upper = np.linspace(pt_te[0], pt_le[0], int(nb_eq_pts_airfoils / 2))
            x_lower = np.linspace(pt_le[0], pt_te[0], int(nb_eq_pts_airfoils / 2))
            x_airfoil = np.concatenate((x_upper, x_lower))

            for k in range(nb_eq_pts_airfoils):
                z_airfoil_interp = pts_airfoil_1_lst[k][2] + (
                        (y_airfoils_ss[i] - pts_airfoil_1_lst[k][1]) / (
                        pts_airfoil_2_lst[k][1] - pts_airfoil_1_lst[k][1])) * (
                                           pts_airfoil_2_lst[k][2] -
                                           pts_airfoil_1_lst[k][2])
                z_airfoil_norm = z_airfoil_interp / np.abs(chord_interp)
                z_airfoil_linear = z_airfoil_norm * np.abs(chord)

                pt = Point(x_airfoil[k], y_airfoils_ss[i], z_airfoil_linear)
                if k == nb_eq_pts_airfoils / 2:
                    pt_le_interp = pt

                pts_airfoil_lst.append(pt)

            crv_airfoil = FittedCurve(pts_airfoil_lst)
            crv_airfoil = TransformedCurve(curve_in = crv_airfoil, from_position = Position(pt_le_interp),
                                           to_position = Position(pt_le))
            rot_angle_z = np.arcsin((pt_te[1] - pt_le[1]) / chord)
            crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                       vector = Vector(0, 0, 1), angle = rot_angle_z)
            crv_airfoil = RotatedCurve(curve_in = crv_airfoil, rotation_point = pt_le,
                                       vector = Vector(pt_te[0] - pt_le[0], pt_te[1] - pt_le[1], pt_te[2] - pt_le[2]),
                                       angle = -np.deg2rad(dihedral_angle))

            crvs_surface_lst.append(crv_airfoil)
            print('Airfoil placed nb=', i)
            i = i + 1
        print('Next section considered, i=', i)
        j = j + 1
        i = i + 1

    crvs_surface_lst.append(crvs_defined_airfoils_ss[-1])

    return crvs_surface_lst


# display([pts_surface_ss_lst,pts_y_airfoils_ss_lst])
# Mirror the points around the symmetrical plane

def crvs_surface_to_points(nb_airfoils_ss, nb_eq_pts_airfoils, crvs_surface_lst):
    pts_surface_lst = []
    print(crvs_surface_lst)
    for i in range(nb_airfoils_ss):
        pts_airfoil_lst = crvs_surface_lst[i].equispaced_points(nb_eq_pts_airfoils)
        pts_surface_lst.append(pts_airfoil_lst)
    return pts_surface_lst


def mirror_object(nb_eq_pts_airfoils=0, pts_surface_ss_lst=[], crvs_engine_l=[], object_type='surface'):
    if object_type == 'engine':
        crvs_engine_r = []
        for i in range(len(crvs_engine_l)):
            crvs_engine_r.append(MirroredCurve(curve_in = crvs_engine_l[i], reference_point = Point(0, 0, 0),
                                               vector1 = Vector(1, 0, 0),
                                               vector2 = Vector(0, 0, 1)))
        return crvs_engine_r

    pts_surface_ss_mirored_lst = []
    for i in range(len(pts_surface_ss_lst)):
        pts_airfoil_mirored_lst = []
        for k in range(nb_eq_pts_airfoils):
            pt_mirored = Point(pts_surface_ss_lst[i][k][0], -pts_surface_ss_lst[i][k][1], pts_surface_ss_lst[i][k][2])
            pts_airfoil_mirored_lst.append(pt_mirored)
        pts_surface_ss_mirored_lst.append(pts_airfoil_mirored_lst)
    return pts_surface_ss_mirored_lst


# pts_surface_lst = pts_surface_ss_lst + pts_surface_ss_mirored_lst
# display(pts_surface_lst)

def create_object_surface(object_type='surface', pts_surface_lst=[], pts_surface_mirored_lst=[], crv_engine_l=[]):
    if object_type == 'engine':
        surf1_l = LoftedSurface(profiles = crv_engine_l[0:3], min_degree = 25)
        surf2_l = LoftedSurface(profiles = crv_engine_l[2:], min_degree = 25)
        return Fused(surf1_l, surf2_l)
    surface_l = FittedSurface(points = pts_surface_lst, min_degree = 8)
    if not pts_surface_ss_mirored_lst == []:
        surface_r = FittedSurface(points = pts_surface_mirored_lst, min_degree = 8)
        surface = Fused(surface_l, surface_r)
        return surface_l, surface_r, surface
    else:
        return surface_ss_l


# wing_surface = FittedSurface(points=pts_surface_lst, tolerance=0.0001)
# display(wing_surface)

def translate_crvs_surface(nb_airfoils_vt, crvs_vt_surface_lst, x_vt, y_vt):
    for i in range(nb_airfoils_vt):
        crvs_vt_surface_lst[i] = TranslatedCurve(crvs_vt_surface_lst[i], Vector(x_vt, y_vt, 10))
    return crvs_vt_surface_lst


def define_attachment_pt_object(crvs_surface_vt_lst, type_object='surface'):
    if type_object == 'engine':
        return crvs_surface_vt_lst[4].midpoint
    else:
        return crvs_surface_vt_lst[0].midpoint


def position_surf1_on_surf2(crvs_surface_vt_lst, surface_ss_l, nb_eq_pts_airfoils_vt=None, nb_airfoils_vt=None,
                            dtheta1=None, dtheta2=None, nb_max_rotations=None):
    # intermediate_crvs1 = []
    # intermediate_crvs2 = []
    if nb_airfoils_vt == None:
        pt_attachment_vt = define_attachment_pt_object(crvs_surface_vt_lst, 'engine')
    else:
        pt_attachment_vt = define_attachment_pt_object(crvs_surface_vt_lst)
    # display([crvs_surface_vt_lst,pt_attachment_vt,surface_ss_l])
    lst_inter_le = LineSegment(start = pt_attachment_vt,
                               end = Point(pt_attachment_vt[0], pt_attachment_vt[1], -1e3)).surface_intersections(
        surface_ss_l)
    pt_attachment_le_vt_projected = lst_inter_le[0]['point']
    vect_displ_vt = Vector(pt_attachment_le_vt_projected[0] - pt_attachment_vt[0],
                           pt_attachment_le_vt_projected[1] - pt_attachment_vt[1],
                           pt_attachment_le_vt_projected[2] - pt_attachment_vt[2])
    # surface_vt = TranslatedSurface(surface_vt,vect_displ_vt)
    # display([surface_wing,surface_vt,surf_engine_l,surf_engine_r])

    # Positioning on the surface of the wing with the attachment point
    for i in range(len(crvs_surface_vt_lst)):
        crvs_surface_vt_lst[i] = TranslatedCurve(crvs_surface_vt_lst[i], vect_displ_vt)
        # intermediate_crvs1.append(crvs_surface_vt_lst[i])
    # display([crvs_surface_vt_lst,surface_ss_l])

    if not nb_airfoils_vt == None:
        # Placement trailing edge
        count = 0
        crvs_vt_rotation_1_lst = []

        lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].start,
                                   end = Point(crvs_surface_vt_lst[0].start[0], crvs_surface_vt_lst[0].start[1],
                                               -1e3)).surface_intersections(
            surface_ss_l)
        pt_attachment_te_vt_projected = lst_inter_te[0]['point']
        dz_old = 1e3
        dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]

        while dz_new < dz_old and count < nb_max_rotations:
            count = count + 1
            print('rotation1 increment nb=', count)
            for i in range(nb_airfoils_vt):
                crvs_surface_vt_lst[i] = RotatedCurve(crvs_surface_vt_lst[i], pt_attachment_le_vt_projected,
                                                      Vector(0, 1, 0), -np.deg2rad(dtheta1))
                # intermediate_crvs1.append(crvs_surface_vt_lst[i])
                crvs_vt_rotation_1_lst.append(crvs_surface_vt_lst[i])
            print('rotation1 increment nb=', count, 'is finished')
            lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].start,
                                       end = Point(crvs_surface_vt_lst[0].start[0], crvs_surface_vt_lst[0].start[1],
                                                   -1e3)).surface_intersections(
                surface_ss_l)
            if lst_inter_te == []:
                break
            pt_attachment_te_vt_projected = lst_inter_te[0]['point']
            dz_old = dz_new
            dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]
        # display([crvs_surface_vt_lst,surface_ss_l,crvs_vt_rotation_1_lst])
        AoA_vt = dtheta1 * count

        # Placement 3rd point non colinear on the based airfoil
        count = 0
        crvs_vt_rotation_2_lst = []

        lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].sample_points[int(115)],
                                   end = Point(crvs_surface_vt_lst[0].sample_points[int(115)][0],
                                               crvs_surface_vt_lst[0].sample_points[int(115)][1],
                                               -1e3)).surface_intersections(
            surface_ss_l)
        pt_attachment_te_vt_projected = lst_inter_te[0]['point']
        dz_old = 1e3
        dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]

        while dz_new < dz_old and count < nb_max_rotations:
            count = count + 1
            print('rotation2 increment nb=', count)
            for i in range(nb_airfoils_vt):
                crvs_surface_vt_lst[i] = RotatedCurve(crvs_surface_vt_lst[i], pt_attachment_le_vt_projected,
                                                      Vector(1, 0, 0), -np.deg2rad(dtheta2))
                # intermediate_crvs2.append(crvs_surface_vt_lst[i])
                crvs_vt_rotation_2_lst.append(crvs_surface_vt_lst[i])
            print('rotation2 increment nb=', count, 'is finished')
            lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].sample_points[int(23 / 2)],
                                       end = Point(crvs_surface_vt_lst[0].sample_points[int(23 / 2)][0],
                                                   crvs_surface_vt_lst[0].sample_points[int(23 / 2)][1],
                                                   -1e3)).surface_intersections(
                surface_ss_l)
            if lst_inter_te == []:
                break
            pt_attachment_te_vt_projected = lst_inter_te[0]['point']
            dz_old = dz_new
            dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]
        # display([crvs_surface_vt_lst,surface_ss_l,crvs_vt_rotation_2_lst])

        return crvs_surface_vt_lst, AoA_vt
    return crvs_surface_vt_lst


def create_engine_l_casing(x_engine, y_engine, fact_engine, AoA_engine):
    with open("coord_engine.dat", 'r') as f:
        crvs_engine_rot = []
        center_engine = Point(120.847 / (2 * fact_engine) + x_engine, y_engine / fact_engine, 0)
        for line in f:
            x, y, z = line.split(' ', 2)
            x = float(x) / fact_engine + x_engine
            y = float(y) / fact_engine + y_engine
            z = float(z.strip()) / fact_engine
            crv = Circle(z - 40 / fact_engine, Point(x, 0, 0))
            crv = RotatedCurve(crv, Point(x, 0, 0), Vector(0, 1, 0), -np.pi / 2)
            crv = RotatedCurve(crv, center_engine, Vector(0, 1, 0), np.deg2rad(AoA_engine))
            crv = TranslatedCurve(curve_in = crv, displacement = Vector(0, 0, 10))
            crv = TranslatedCurve(curve_in = crv, displacement = Vector(0, y, 0))
            crvs_engine_rot.append(crv)

        pt_attachment_engine = crvs_engine_rot[4].midpoint
        crvs_engine_l = crvs_engine_rot
        for i in range(len(crvs_engine_l)):
            crvs_engine_l[i] = RotatedCurve(crvs_engine_rot[i], pt_attachment_engine, Vector(0, 0, 1), -np.pi)

    return crvs_engine_l


def create_base_projection_engine(nb_eq_pts_airfoils_vb, pts_surface_vb_lst, x_vb, y_vb, chords_vb, span_vb, AoA_vb):
    pt_center_rect_placement_engine = Point(x_vb - chords_vb[0] / 2, y_vb,
                                            pts_surface_vb_lst[-1][int(nb_eq_pts_airfoils_vb / 2)][2] - span_vb / 5)
    surface_base_rect_placement_engine = RectangularSurface(width = chords_vb[0] * 2, length = chords_vb[0],
                                                            position = Position(pt_center_rect_placement_engine))
    surface_base_rect_placement_engine = RotatedSurface(surface_base_rect_placement_engine,
                                                        pt_center_rect_placement_engine, Vector(0, 1, 0),
                                                        -np.deg2rad(AoA_vb))
    return surface_base_rect_placement_engine


def position_engine_on_vertical_base(x_vb, y_vb, span_vb, pts_surface_vb_lst, chords_vb, nb_eq_pts_airfoils_vb,
                                     pt_attachment_engine, crvs_engine_l, AoA_vb):
    pt_center_rect_placement_engine = Point(x_vb - chords_vb[0] / 2, y_vb,
                                            pts_surface_vb_lst[-1][int(nb_eq_pts_airfoils_vb / 2)][2] - span_vb / 2)
    surface_base_rect_placement_engine = RectangularSurface(width = chords_vb[0] * 2, length = chords_vb[0],
                                                            position = Position(pt_center_rect_placement_engine))
    surface_base_rect_placement_engine = RotatedSurface(surface_base_rect_placement_engine,
                                                        pt_center_rect_placement_engine, Vector(0, 1, 0),
                                                        -np.deg2rad(AoA_vb))
    lst_inter_te = LineSegment(start = pt_attachment_engine,
                               end = Point(pt_attachment_engine[0], pt_attachment_engine[1], -5)).surface_intersections(
        surface_base_rect_placement_engine)
    pt_attachment_engine_projected = lst_inter_te[0]['point']
    vect_displ_engine = Vector(pt_attachment_engine_projected[0] - pt_attachment_engine[0],
                               pt_attachment_engine_projected[1] - pt_attachment_engine[1],
                               pt_attachment_engine_projected[2] - pt_attachment_engine[2])

    for i in range(len(crvs_engine_l)):
        crvs_engine_l[i] = TranslatedCurve(curve_in = crvs_engine_l[i], displacement = vect_displ_engine)

    return crvs_engine_l

def dihedral_angles_to_z(nb_pts_dihedral,angles_dihedral,spans_dihedral):
    y_dihedral = [0, ]  # location of the definition of dihedral. 1st and last points must be defined at root and tip
    z_dihedral = [0, ]
    print(y_dihedral)
    for i in range(nb_pts_dihedral - 1):
        print([y_dihedral[-1] + cos(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]])
        y_dihedral += [y_dihedral[-1] + cos(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]]
        z_dihedral += [z_dihedral[-1] + sin(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]]
    return y_dihedral,z_dihedral



# Data BWB
nb_eq_pts_rails = 200
nb_eq_pts_airfoils = 200
round_nb = 7

nb_fixed_control_points = 3
sweeps_le = [0, 40, 0]
sweeps_te = [40, 45, 45]
sspan = 6
x_fcp = [0, -.5, -1]
y_fcp = [0, 4, sspan]
chords = [2, 1, .5]
rho_le = [.5, .5, .5, .5]
rho_te = [.5, .5, .5, .5]

nb_pts_twist = 3
angles_twist = [0, 0, 0]  # degrees of twist
y_twist = [0, 3, 6]  # location of the definition of twist
nb_pts_twist_axis = 3
x_twist_axis = [0, -1, -1]
y_twist_axis = [0, 3, sspan]  # 1st and last points must be at the root and tip

nb_pts_dihedral = 4
angles_dihedral = [10, 10, 10, 89]  # angle defined on the right of the section
spans_dihedral = [2.75, 2.75, 0.5]

y_dihedral = [0, ]  # location of the definition of dihedral. 1st and last points must be defined at root and tip
z_dihedral = [0, ]
print(y_dihedral)
for i in range(nb_pts_dihedral - 1):
    print([y_dihedral[-1] + cos(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]])
    y_dihedral += [y_dihedral[-1] + cos(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]]
    z_dihedral += [z_dihedral[-1] + sin(np.deg2rad(angles_dihedral[i + 1])) * spans_dihedral[i]]

print("y_dihedral=", y_dihedral)
print("z_dihedral=", z_dihedral)
# z_dihedral = [2, 2.5, 3, 3.5]
# y_dihedral = [0, 3, 5.99,6]

# location of the definition of dihedral. 1st and last points must be defined at root and tip
weights_dihedral = [.5, .5, .5, .5, .5, .5, 1000, 1000, 1000, 1000]

nb_defined_airfoils = 5
airfoils_names = ["whitcomb", "24012", "whitcomb", "24012", "whitcomb"]
y_defined_airfoils = [0, 3, 4, 5, y_dihedral[-1]]  # 1st and last correspond to root and tip
cant_angles_airfoils = [0, 0, 0, 0, 0]  # degree
dihedral_angles_airfoils = [0, 0, 0, 10, 90]

nb_airfoils_ss = 50

tol = 1e-5

nb_fixed_control_points_vt = 2
sweeps_vt_le = [50, 50]
sweeps_vt_te = [10, 10]
x_vt = -1.1
y_vt = 1.5
span_vt = 0.4
chords_vt = [0.5, 0.4]
x_fcp_vt = [0, -tan(np.deg2rad(sweeps_vt_le[0])) * span_vt]
y_fcp_vt = [0, span_vt]
rho_vt_le = [.01, .01]
rho_vt_te = [.01, .01]

nb_pts_twist_vt = 2
angles_twist_vt = [0, 0]  # degrees of twist
y_twist_vt = [0, span_vt]  # location of the definition of twist
nb_eq_pts_rails_vt = 30
nb_eq_pts_airfoils_vt = 200
nb_pts_twist_axis_vt = 2
x_twist_axis_vt = [0, -tan(np.deg2rad(sweeps_vt_le[0])) * span_vt - chords_vt[1]]
y_twist_axis_vt = y_fcp_vt  # 1st and last points must be at the root and tip

nb_pts_dihedral_vt = 2
angle_dihedral_vt = [80, 80]  # degree
y_dihedral_vt = [0, cos(np.deg2rad(angle_dihedral_vt[
                                       -1])) * span_vt]  # location of the definition of dihedral. 1st and last points must be defined at root and tip
z_dihedral_vt = [0.0, sin(np.deg2rad(angle_dihedral_vt[
                                         -1])) * span_vt]
weights_dihedral_vt = [.1, .1, .1, .1]

nb_defined_airfoils_vt = 2
airfoils_names_vt = ["24012", "24012"]
y_defined_airfoils_vt = y_dihedral_vt  # 1st and last correspond to root and tip
cant_angles_airfoils_vt = [0, 0]  # degree
dihedral_angles_airfoils_vt = [90, 90]

nb_airfoils_vt = 2

dtheta1 = 0.1  # degree
dtheta2 = 0.001  # degree
nb_max_rotations = 500

nb_fixed_control_points_vb = 2
sweeps_vb_le = [72, 72]
lambda_vb_le = np.deg2rad(sweeps_vb_le)
x_vb = -1.1
y_vb = 0.6
span_vb = 0.075
y_fcp_vb = [0, span_vb]
chords_vb = [0.5, 0.2]

x_fcp_vb = [0, -tan(lambda_vb_le[0]) * span_vb]
y_fcp_vb = [0, span_vb]

nb_pts_dihedral_vb = 2
angle_dihedral_vb = [110, 110]  # degree
y_dihedral_vb = [0, cos(np.deg2rad(angle_dihedral_vb[
                                       -1])) * span_vb]  # location of the definition of dihedral. 1st and last points must be defined at root and tip
z_dihedral_vb = [0.0, sin(np.deg2rad(angle_dihedral_vb[
                                         -1])) * span_vb]
print(y_dihedral_vb)
nb_eq_pts_rails_vb = 2
nb_eq_pts_airfoils_vb = 200

nb_defined_airfoils_vb = 2
nb_airfoils_vb = 2
airfoils_names_vb = ["24030", "24030"]
y_defined_airfoils_vb = y_dihedral_vb  # 1st and last correspond to root and tip
dihedral_angles_airfoils_vb = [90, 90]
crvs_airfoils_vb_lst = []
pts_vb = []
chords_defined_airfoils_vb = []

fact_engine = 180
y_engine = y_vb + y_fcp_vb[0]
x_engine = x_vb + x_fcp_vb[0] - chords_vb[0]
AoA_engine = 10  # degree

# End Data BWB

# Define wing
[eq_pts_le_ss_lst, eq_pts_te_ss_lst] = sweep_crvs(nb_eq_pts_rails, nb_fixed_control_points, y_fcp, chords, sweeps_le,
                                                  sweeps_te, rho_le, rho_te, 'wing', x_fcp)
[eq_pts_twist_distrib_lst, eq_pts_twist_axis_lst] = calculate_twist_distrib(nb_eq_pts_rails, nb_pts_twist, angles_twist,
                                                                            y_twist, nb_pts_twist_axis, x_twist_axis,
                                                                            y_twist_axis, True)
[eq_pts_dihedral_distrib_lst, crv_dihedral_distrib] = calculate_dihedral_distrib(nb_eq_pts_rails, nb_pts_dihedral,
                                                                                 z_dihedral, y_dihedral,
                                                                                 weights_dihedral, True)
[eq_pts_le_ss_lst, eq_pts_te_ss_lst] = add_twist(nb_eq_pts_rails, eq_pts_le_ss_lst, eq_pts_te_ss_lst,
                                                 eq_pts_twist_distrib_lst, eq_pts_twist_axis_lst, True)
[eq_pts_le_ss_lst, eq_pts_te_ss_lst] = add_dihedral(nb_eq_pts_rails, eq_pts_le_ss_lst, eq_pts_te_ss_lst,
                                                    eq_pts_dihedral_distrib_lst, True)
crvs_defined_airfoils_ss = place_defined_airfoils(nb_eq_pts_rails, nb_defined_airfoils, airfoils_names,
                                                  y_defined_airfoils, dihedral_angles_airfoils,
                                                  eq_pts_le_ss_lst, eq_pts_te_ss_lst, cant_angles_airfoils)
# display([crvs_defined_airfoils_ss,eq_pts_le_ss_lst, eq_pts_te_ss_lst])
y_airfoils_ss = determine_position_intermediate_airfoils(nb_airfoils_ss, nb_defined_airfoils, y_defined_airfoils,
                                                         crv_dihedral_distrib)
crvs_surface_lst = place_intermediate_airfoils('surface', crvs_defined_airfoils_ss, nb_eq_pts_rails, nb_eq_pts_airfoils,
                                               nb_defined_airfoils,
                                               y_defined_airfoils,
                                               cant_angles_airfoils, dihedral_angles_airfoils,
                                               y_airfoils_ss,
                                               eq_pts_le_ss_lst, eq_pts_te_ss_lst, tol, round_nb)
# display([crvs_surface_lst,eq_pts_le_ss_lst, eq_pts_te_ss_lst])
pts_surface_ss_lst = crvs_surface_to_points(nb_airfoils_ss, nb_eq_pts_airfoils, crvs_surface_lst)
pts_surface_ss_mirored_lst = mirror_object(nb_eq_pts_airfoils, pts_surface_ss_lst)
[surface_ss_l, surface_ss_r, surface_wing] = create_object_surface('surface', pts_surface_ss_lst,
                                                                   pts_surface_ss_mirored_lst)
# display(surface_wing)

# Define VT

[eq_pts_le_vt_lst, eq_pts_te_vt_lst] = sweep_crvs(nb_eq_pts_rails_vt, nb_fixed_control_points_vt, y_fcp_vt, chords_vt,
                                                  sweeps_vt_le,
                                                  sweeps_vt_te, rho_vt_le, rho_vt_te, 'vt', x_fcp_vt)
display([eq_pts_le_vt_lst, eq_pts_te_vt_lst])
[eq_pts_twist_distrib_vt_lst, eq_pts_twist_axis_vt_lst] = calculate_twist_distrib(nb_eq_pts_rails_vt, nb_pts_twist_vt,
                                                                                  angles_twist_vt,
                                                                                  y_twist_vt, nb_pts_twist_axis_vt,
                                                                                  x_twist_axis_vt,
                                                                                  y_twist_axis_vt, True)
[eq_pts_dihedral_distrib_vt_lst, crv_dihedral_distrib_vt] = calculate_dihedral_distrib(nb_eq_pts_rails_vt,
                                                                                       nb_pts_dihedral_vt,
                                                                                       z_dihedral_vt, y_dihedral_vt,
                                                                                       weights_dihedral_vt, True)
[eq_pts_le_vt_lst, eq_pts_te_vt_lst] = add_twist(nb_eq_pts_rails_vt, eq_pts_le_vt_lst, eq_pts_te_vt_lst,
                                                 eq_pts_twist_distrib_vt_lst, eq_pts_twist_axis_vt_lst, True)
[eq_pts_le_vt_lst, eq_pts_te_vt_lst] = add_dihedral(nb_eq_pts_rails_vt, eq_pts_le_vt_lst, eq_pts_te_vt_lst,
                                                    eq_pts_dihedral_distrib_vt_lst, True)
display([eq_pts_le_vt_lst, eq_pts_te_vt_lst])
crvs_defined_airfoils_vt = place_defined_airfoils(nb_eq_pts_rails_vt, nb_defined_airfoils_vt, airfoils_names_vt,
                                                  y_defined_airfoils_vt, dihedral_angles_airfoils_vt,
                                                  eq_pts_le_vt_lst, eq_pts_te_vt_lst, cant_angles_airfoils_vt)
y_airfoils_vt = determine_position_intermediate_airfoils(nb_airfoils_vt, nb_defined_airfoils_vt, y_defined_airfoils_vt,
                                                         crv_dihedral_distrib_vt)
crvs_surface_vt_lst = place_intermediate_airfoils('surface', crvs_defined_airfoils_vt, nb_eq_pts_rails_vt,
                                                  nb_eq_pts_airfoils_vt, nb_defined_airfoils_vt,
                                                  y_defined_airfoils_vt,
                                                  cant_angles_airfoils_vt, dihedral_angles_airfoils_vt,
                                                  y_airfoils_vt,
                                                  eq_pts_le_vt_lst, eq_pts_te_vt_lst, tol, round_nb)
crvs_vt_surface_lst = translate_crvs_surface(nb_airfoils_vt, crvs_surface_vt_lst, x_vt, y_vt)
pt_attachment_vt = define_attachment_pt_object(crvs_surface_vt_lst, type_object = 'surface')
[crvs_surface_vt_lst, AoA_vt] = position_surf1_on_surf2(crvs_surface_vt_lst, surface_ss_l, nb_eq_pts_airfoils_vt,
                                                        nb_airfoils_vt, dtheta1, dtheta2, nb_max_rotations)
pts_surface_vt_lst = crvs_surface_to_points(nb_airfoils_vt, nb_eq_pts_airfoils_vt, crvs_surface_vt_lst)
surface_vt_l = create_object_surface('surface', pts_surface_vt_lst)
pts_surface_vt_mirorred_lst = mirror_object(nb_eq_pts_airfoils_vt, pts_surface_vt_lst)
surface_vt_r = create_object_surface('surface', pts_surface_vt_mirorred_lst)

# Define VB

[eq_pts_le_vb_lst, eq_pts_te_vb_lst] = sweep_crvs(nb_eq_pts_rails_vb, nb_fixed_control_points_vb, y_fcp_vb, chords_vb,
                                                  sweeps_vb_le, surface_type = 'vb')
[] = calculate_twist_distrib(nb_eq_pts_rails_vb, False)
[eq_pts_dihedral_distrib_vb_lst, crv_dihedral_distrib_vb] = calculate_dihedral_distrib(nb_eq_pts_rails_vb,
                                                                                       nb_pts_dihedral_vb,
                                                                                       z_dihedral_vb, y_dihedral_vb,
                                                                                       True)
[eq_pts_le_vb_lst, eq_pts_te_vb_lst] = add_twist(nb_eq_pts_rails_vb, eq_pts_le_vb_lst, eq_pts_te_vb_lst, False)
[eq_pts_le_vb_lst, eq_pts_te_vb_lst] = add_dihedral(nb_eq_pts_rails_vb, eq_pts_le_vb_lst, eq_pts_te_vb_lst,
                                                    eq_pts_dihedral_distrib_vb_lst, True)
crvs_defined_airfoils_vb = place_defined_airfoils(nb_eq_pts_rails_vb, nb_defined_airfoils_vb, airfoils_names_vb,
                                                  y_defined_airfoils_vb, dihedral_angles_airfoils_vb,
                                                  eq_pts_le_vb_lst, eq_pts_te_vb_lst)
y_airfoils_vb = determine_position_intermediate_airfoils(nb_airfoils_vb, nb_defined_airfoils_vb, y_defined_airfoils_vb,
                                                         crv_dihedral_distrib_vb)
crvs_surface_vb_lst = place_intermediate_airfoils('vb', crvs_defined_airfoils_vb)
crvs_surface_vb_lst = translate_crvs_surface(nb_airfoils_vb, crvs_surface_vb_lst, x_vb, y_vb)
[crvs_surface_vb_lst, AoA_vb] = position_surf1_on_surf2(crvs_surface_vb_lst, surface_ss_l, nb_airfoils_vb, dtheta1,
                                                        dtheta2, nb_max_rotations)
pts_surface_vb_lst = crvs_surface_to_points(nb_airfoils_vb, nb_eq_pts_airfoils_vb, crvs_surface_vb_lst)
surface_vb_l = create_object_surface('surface', pts_surface_vb_lst)

# Define engine

crvs_engine_l = create_engine_l_casing(x_engine, y_engine, fact_engine, AoA_engine)
surface_engine_l = create_object_surface('engine', [], [], crvs_engine_l)
surface_base_rect_placement_engine = create_base_projection_engine(nb_eq_pts_airfoils_vb, pts_surface_vb_lst, x_vb,
                                                                   y_vb, chords_vb, span_vb, AoA_vb)
crvs_engine_l = position_surf1_on_surf2(crvs_engine_l, surface_base_rect_placement_engine)
surface_engine_l = create_object_surface('engine', [], [], crvs_engine_l)
display([surface_wing, surface_vt_l, surface_vb_l, surface_base_rect_placement_engine, surface_engine_l])
