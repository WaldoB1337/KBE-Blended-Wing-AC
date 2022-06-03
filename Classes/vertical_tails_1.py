from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from copy import deepcopy
from parapy.gui import display

import numpy as np
from math import *

from surface_1 import Surface_BWB


class VerticalTail(GeomBase):
    nb_eq_pts_rails = Input()
    nb_eq_pts_airfoils = Input()

    nb_fixed_control_points = Input()
    sweeps_le = Input()
    span = Input()
    x_fcp_except_last = Input()
    y_fcp_except_last = Input()
    chords = Input()

    nb_pts_dihedral = Input()
    angles_dihedral = Input()  # angle defined on the right of the section
    spans_dihedral_except_last = Input()

    nb_defined_airfoils = Input()
    airfoils_names = Input()
    y_defined_airfoils_except_last = Input()  # 1st and last correspond to root and tip
    cant_angles_airfoils = Input()  # degree
    dihedral_angles_airfoils = Input()
    heights_factors = Input()

    nb_airfoils = Input()

    surface_type = Input()
    is_twisted = Input()
    is_there_dihedral = Input()

    nb_eq_pts_airfoils_vt = Input()
    x_vt = Input()
    y_vt = Input()

    dtheta1 = Input()
    dtheta2 = Input()
    nb_max_rotations = Input()

    surf2_type = Input()
    surfaces = Input()

    is_there_vt = Input()

    name = Input()
    nb_chordwise_vortices = Input()
    nb_spanwise_vortices = Input()
    is_mirrored = Input()


    @Attribute
    def initial_surface(self):
        print("Vertical Tail placement started")
        return Surface_BWB(nb_eq_pts_rails = self.nb_eq_pts_rails,
                           nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                           nb_fixed_control_points = self.nb_fixed_control_points,
                           sweeps_le = self.sweeps_le,
                           sweeps_te = [10, 10],
                           sspan = self.span,
                           x_fcp_except_last = self.x_fcp_except_last,
                           y_fcp_except_last = self.y_fcp_except_last,
                           chords = self.chords,
                           rho_le = [.01, .01],
                           rho_te = [.01, .01],
                           nb_pts_twist = 2,
                           angles_twist = [0, 0],
                           y_twist_except_last = [0, 0],
                           nb_pts_twist_axis = 2,
                           x_twist_axis = [0, -0.876701437037684],
                           y_twist_axis_except_last = [0, 0.4],
                           nb_pts_dihedral = self.nb_pts_dihedral,
                           angles_dihedral = self.angles_dihedral,
                           spans_dihedral_except_last = self.spans_dihedral,
                           weights_dihedral = [.1, .1, .1, .1],
                           nb_defined_airfoils = self.nb_defined_airfoils,
                           airfoils_names = self.airfoils_names,
                           y_defined_airfoils_except_last = self.y_defined_airfoils_except_last,
                           cant_angles_airfoils = self.cant_angles_airfoils,
                           dihedral_angles_airfoils = self.dihedral_angles_airfoils,
                           heights_factors = self.heights_factors,
                           nb_airfoils_ss = self.nb_airfoils,
                           surface_type = self.surface_type,
                           is_twisted = self.is_twisted,
                           is_there_dihedral = self.is_there_dihedral,
                           nb_sections_AVL = 2)

    @Attribute
    def x_fcp(self):
        return self.initial_surface.x_fcp

    @Attribute
    def spans_dihedral(self):
        return self.spans_dihedral_except_last + [self.span]

    @Attribute
    def y_defined_airfoils(self):
        return self.initial_surface.y_defined_airfoils

    @Attribute
    def y_fcp(self):
        return self.y_fcp_except_last + [self.span]

    @Attribute
    def crvs_vt_surface_lst(self):
        return self.initial_surface.crvs_surface_all_airfoils

    @Attribute
    def translate_crvs_surface(self):
        crvs_vt_surface_lst = self.crvs_vt_surface_lst[:]
        for i in range(self.nb_airfoils):
            crvs_vt_surface_lst[i] = TranslatedCurve(crvs_vt_surface_lst[i], Vector(self.x_vt, self.y_vt, 10))
        return crvs_vt_surface_lst

    @Attribute
    def define_attachment_pt_object(self):
        return self.translate_crvs_surface[0].midpoint

    @Attribute
    def position_surf1_on_surf2(self):
        pt_attachment_vt = self.define_attachment_pt_object
        crvs_surface_vt_lst = self.translate_crvs_surface[:]
        # display([crvs_surface_vt_lst,pt_attachment_vt,surface_ss_l])
        if self.surf2_type == 'wing':
            surf2 = self.surfaces[0]
        elif self.surf2_type == 'engine':
            surf2 = self.surfaces[1]
        lst_inter_le = LineSegment(start = pt_attachment_vt,
                                   end = Point(pt_attachment_vt[0], pt_attachment_vt[1],
                                               -1e3)).surface_intersections(
            surf2)
        if lst_inter_le == []:
            display([LineSegment(start = pt_attachment_vt,
                                 end = Point(pt_attachment_vt[0], pt_attachment_vt[1], -1e3)), surf2])
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

        # Placement trailing edge
        count = 0
        crvs_vt_rotation_1_lst = []

        lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].start,
                                   end = Point(crvs_surface_vt_lst[0].start[0], crvs_surface_vt_lst[0].start[1],
                                               -1e3)).surface_intersections(
            surf2)
        pt_attachment_te_vt_projected = lst_inter_te[0]['point']
        dz_old = 1e3
        dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]

        while dz_new < dz_old and count < self.nb_max_rotations:
            count = count + 1
            #print('rotation1 increment nb=', count)
            for i in range(self.nb_airfoils):
                crvs_surface_vt_lst[i] = RotatedCurve(crvs_surface_vt_lst[i], pt_attachment_le_vt_projected,
                                                      Vector(0, 1, 0), -np.deg2rad(self.dtheta1))
                # intermediate_crvs1.append(crvs_surface_vt_lst[i])
                crvs_vt_rotation_1_lst.append(crvs_surface_vt_lst[i])
            #print('rotation1 increment nb=', count, 'is finished')
            lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].start,
                                       end = Point(crvs_surface_vt_lst[0].start[0], crvs_surface_vt_lst[0].start[1],
                                                   -1e3)).surface_intersections(
                surf2)
            if lst_inter_te == []:
                break
            pt_attachment_te_vt_projected = lst_inter_te[0]['point']
            dz_old = dz_new
            dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]
            # display([crvs_surface_vt_lst,surface_ss_l,crvs_vt_rotation_1_lst])
        AoA_vt = self.dtheta1 * count

        # Placement 3rd point non colinear on the based airfoil
        count = 0
        crvs_vt_rotation_2_lst = []

        lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].sample_points[int(115)],
                                   end = Point(crvs_surface_vt_lst[0].sample_points[int(115)][0],
                                               crvs_surface_vt_lst[0].sample_points[int(115)][1],
                                               -1e3)).surface_intersections(
            surf2)
        pt_attachment_te_vt_projected = lst_inter_te[0]['point']
        dz_old = 1e3
        dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]

        while dz_new < dz_old and count < self.nb_max_rotations:
            count = count + 1
            #print('rotation2 increment nb=', count)
            for i in range(self.nb_airfoils):
                crvs_surface_vt_lst[i] = RotatedCurve(crvs_surface_vt_lst[i], pt_attachment_le_vt_projected,
                                                      Vector(1, 0, 0), -np.deg2rad(self.dtheta2))
                # intermediate_crvs2.append(crvs_surface_vt_lst[i])
                crvs_vt_rotation_2_lst.append(crvs_surface_vt_lst[i])
            #print('rotation2 increment nb=', count, 'is finished')
            lst_inter_te = LineSegment(start = crvs_surface_vt_lst[0].sample_points[int(23 / 2)],
                                       end = Point(crvs_surface_vt_lst[0].sample_points[int(23 / 2)][0],
                                                   crvs_surface_vt_lst[0].sample_points[int(23 / 2)][1],
                                                   -1e3)).surface_intersections(
                surf2)
            if lst_inter_te == []:
                break
            pt_attachment_te_vt_projected = lst_inter_te[0]['point']
            dz_old = dz_new
            dz_new = crvs_surface_vt_lst[0].start[2] - pt_attachment_te_vt_projected[2]
            # display([crvs_surface_vt_lst,surface_ss_l,crvs_vt_rotation_2_lst])

        return crvs_surface_vt_lst, AoA_vt

    @Attribute
    def crvs_surface_to_points(self):
        crvs_surface_lst = self.position_surf1_on_surf2[0]
        pts_surface_lst = []
        for i in range(self.nb_airfoils):
            pts_airfoil_lst = crvs_surface_lst[i].equispaced_points(self.nb_eq_pts_airfoils_vt)
            pts_surface_lst.append(pts_airfoil_lst)
        return pts_surface_lst

    @Attribute
    def mirror_points(self):
        pts_surface_lst = self.crvs_surface_to_points
        pts_surface_mirored_lst = []
        for i in range(len(pts_surface_lst)):
            pts_airfoil_mirored_lst = []
            for k in range(self.nb_eq_pts_airfoils):
                pt_mirored = Point(pts_surface_lst[i][k][0], -pts_surface_lst[i][k][1],
                                   pts_surface_lst[i][k][2])
                pts_airfoil_mirored_lst.append(pt_mirored)
            pts_surface_mirored_lst.append(pts_airfoil_mirored_lst)
        print("Vertical tail placed")
        return pts_surface_mirored_lst

    @Part
    def surface(self):
        return FittedSurface(points = self.crvs_surface_to_points, min_degree = 8)

    @Part
    def surface_mirorred(self):
        return FittedSurface(points = self.mirror_points, min_degree = 8)

    @Part
    def avl_surface(self):
        return avl.Surface(name = self.name,
                           n_chordwise = self.nb_chordwise_vortices,
                           chord_spacing = avl.Spacing.cosine,
                           n_spanwise = self.nb_spanwise_vortices,
                           span_spacing = avl.Spacing.cosine,
                           y_duplicate = self.position.point[1] if self.is_mirrored else None,
                           sections = self.position_surf1_on_surf2[0] if self.is_there_vt else None)

if __name__ == '__main__':
    from parapy.gui import display

    obj = VerticalTail()
    display(obj)
