from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from copy import deepcopy
from parapy.gui import display

import numpy as np
from math import *

from surface_1 import Surface_BWB
from engine_1 import Engine
from vertical_tails_1 import VerticalTail
from vertical_base_1 import VerticalBase


class BWB(GeomBase):
    nb_eq_pts_rails = Input(200)
    nb_eq_pts_airfoils = Input(200)

    nb_fixed_control_points_wing = Input(4)
    sweeps_le_wing = Input([25, 50, 55, 40])
    sweeps_te_wing = Input([-40, -40, 0, 30])
    sspan_wing = Input(7)
    x_fcp_wing = Input([0, -0.75, -2.5, -6])
    y_fcp_wing_except_last = Input([0, 1, 2.5])
    chords_wing = Input([7, 7.1, 3.5, 0.8])
    rho_le_wing = Input([.5, .5, .5, .5])
    rho_te_wing = Input([.5, .5, .5, .5])

    nb_pts_twist_wing = Input(3)
    angles_twist_wing = Input([0, 0, 0])  # degrees of twist
    y_twist_wing_except_last = Input([0, 3])  # location of the definition of twist
    nb_pts_twist_axis_wing = Input(3)
    x_twist_axis_wing = Input([0, -1, -1])
    y_twist_axis_wing_except_last = Input([0, 3])  # 1st and last points must be at the root and tip

    nb_pts_dihedral_wing = Input(4)
    angles_dihedral_wing = Input([0, 0, 0, 0])  # angle defined on the right of the section
    spans_dihedral_wing_except_last = Input([2.75, 2.75])

    weights_dihedral_wing = Input([.5, .5, .5, .5, .5, .5, 1000, 1000, 1000, 1000])

    nb_defined_airfoils_wing = Input(5)
    airfoils_names_wing = Input(["23012", "23012", "23012", "23012", "0003"])
    y_defined_airfoils_wing_except_last = Input([0, 3, 4, 5])  # 1st and last correspond to root and tip
    cant_angles_airfoils_wing = Input([0, 0, 0, 0, 0])  # degree
    dihedral_angles_airfoils_wing = Input([0, 0, 0, 0, 0])
    heights_factors_wing = Input([1.5, 1.2, 1, 1, 1])

    nb_airfoils_ss_wing = Input(50)

    surface_type_wing = Input('wing')
    is_twisted_wing = Input(True)
    is_there_dihedral_wing = Input(True)

    name_wing = Input('Main_wing')
    nb_chordwise_vortices_wing = Input(30)
    nb_spanwise_vortices_wing = Input(50)
    is_mirrored_wing = Input(True)
    nb_sections_AVL_wing = Input(5)

#----------------------------------------------

    nb_eq_pts_rails_vt = Input(2)
    nb_fixed_control_points_vt = Input(2)
    sweeps_le_vt = Input([50, 50])
    span_vt= Input(0.4)
    x_fcp_vt_except_last = Input([0])
    y_fcp_vt_except_last = Input([0])
    chords_vt = Input([0.5, 0.4])

    nb_pts_dihedral_vt = Input(2)
    angles_dihedral_vt = Input([80, 80])  # angle defined on the right of the section
    spans_dihedral_vt_except_last = Input([])

    nb_defined_airfoils_vt = Input(2)
    airfoils_names_vt = Input(["24012", "24012"])
    y_defined_airfoils_except_last_vt = Input([0])  # 1st and last correspond to root and tip
    cant_angles_airfoils_vt = Input([0, 0])  # degree
    dihedral_angles_airfoils_vt = Input([90, 90])
    heights_factors_vt = Input([1, 1])

    nb_airfoils_vt = Input(2)

    surface_type_vt = Input('vt')
    is_twisted_vt = Input(False)
    is_there_dihedral_vt = Input(True)

    x_vt = Input(0)
    y_vt = Input(0)

    dtheta1 = Input(0.1)
    dtheta2 = Input(0.01)
    nb_max_rotations = Input(200)

    surf2_type_vt = Input('wing')

    is_there_vt = Input(False)

    name_vt = Input('tail')
    nb_chordwise_vortices_vt = Input(30)
    nb_spanwise_vortices_vt = Input(50)

#---------------------------------------------------

    nb_eq_pts_rails_vb = Input(2)
    nb_fixed_control_points_vb = Input(2)
    sweeps_le_vb = Input([72, 72])
    span_vb = Input(0.25)
    x_fcp_vb_except_last = Input([0])
    y_fcp_vb_except_last = Input([0])
    chords_vb = Input([0.75, 0.3])

    nb_pts_dihedral_vb = Input(2)
    angles_dihedral_vb = Input([100, 100])  # angle defined on the right of the section
    spans_dihedral_vb_except_last = Input([])

    nb_defined_airfoils_vb = Input(2)
    airfoils_names_vb = Input(["24030", "24030"])
    y_defined_airfoils_except_last_vb = Input([0])  # 1st and last correspond to root and tip
    cant_angles_airfoils_vb = Input([0, 0])  # degree
    dihedral_angles_airfoils_vb = Input([90, 90])
    heights_factors_vb = Input([1, 1])

    nb_airfoils_vb = Input(2)

    surface_type_vb = Input('vb')
    is_twisted_vb = Input(False)
    is_there_dihedral_vb = Input(True)

    x_vb = Input(-6.23)
    y_vb = Input(0.9)

    surf2_type_vb = Input('wing')

    name_vb = Input('base')
    nb_chordwise_vortices_vb = Input(30)
    nb_spanwise_vortices_vb = Input(50)

#--------------------------------

    is_there_engine = Input(True)
    fact_engine = Input(110)
    AoA_engine = Input(20)

#--------------------------------

    Mach = Input(0.11)

# --------------------------------

    @Part
    def wing(self):
        return Surface_BWB(nb_eq_pts_rails = self.nb_eq_pts_rails,
                          nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                          nb_fixed_control_points = self.nb_fixed_control_points_wing,
                          sweeps_le = self.sweeps_le_wing,
                          sweeps_te = self.sweeps_te_wing,
                          sspan = self.sspan_wing,
                          x_fcp_except_last = self.x_fcp_wing,
                          y_fcp_except_last = self.y_fcp_wing_except_last,
                          chords = self.chords_wing,
                          rho_le = self.rho_le_wing,
                          rho_te = self.rho_te_wing,
                          nb_pts_twist = self.nb_pts_twist_wing,
                          angles_twist = self.angles_twist_wing,
                          y_twist_except_last = self.y_twist_wing_except_last,
                          nb_pts_twist_axis = self.nb_pts_twist_axis_wing,
                          x_twist_axis = self.x_twist_axis_wing,
                          y_twist_axis_except_last = self.y_twist_axis_wing_except_last,
                          nb_pts_dihedral = self.nb_pts_dihedral_wing,
                          angles_dihedral = self.angles_dihedral_wing,
                          spans_dihedral_except_last = self.spans_dihedral_wing_except_last,
                          weights_dihedral = self.weights_dihedral_wing,
                          nb_defined_airfoils = self.nb_defined_airfoils_wing,
                          airfoils_names = self.airfoils_names_wing,
                          y_defined_airfoils_except_last = self.y_defined_airfoils_wing_except_last,
                          cant_angles_airfoils = self.cant_angles_airfoils_wing,
                          dihedral_angles_airfoils = self.dihedral_angles_airfoils_wing,
                          heights_factors = self.heights_factors_wing,
                          nb_airfoils_ss = self.nb_airfoils_ss_wing,
                          surface_type = self.surface_type_wing,
                          is_twisted = self.is_twisted_wing,
                          is_there_dihedral = self.is_there_dihedral_wing,
                          name = self.name_wing,
                          nb_chordwise_vortices = self.nb_chordwise_vortices_wing,
                          nb_spanwise_vortices = self.nb_spanwise_vortices_wing,
                          is_mirrored = self.is_mirrored_wing,
                          nb_sections_AVL = self.nb_sections_AVL_wing)

    @Part
    def vertical_bases(self):
            return VerticalBase(nb_eq_pts_rails = self.nb_eq_pts_rails_vb,
                                nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                                nb_fixed_control_points = self.nb_fixed_control_points_vb,
                                sweeps_le = self.sweeps_le_vb,
                                span = self.span_vb,
                                x_fcp_except_last = self.x_fcp_vb_except_last,
                                y_fcp_except_last = self.y_fcp_vb_except_last,
                                chords = self.chords_vb,
                                nb_pts_dihedral = self.nb_pts_dihedral_vb,
                                angles_dihedral = self.angles_dihedral_vb,
                                spans_dihedral_except_last = self.spans_dihedral_vb_except_last,
                                nb_defined_airfoils = self.nb_defined_airfoils_vb,
                                airfoils_names = self.airfoils_names_vb,
                                y_defined_airfoils_except_last = self.y_defined_airfoils_except_last_vb,
                                cant_angles_airfoils = self.cant_angles_airfoils_vb,
                                dihedral_angles_airfoils = self.dihedral_angles_airfoils_vb,
                                heights_factors = self.heights_factors_vb,
                                nb_airfoils = self.nb_airfoils_vb,
                                surface_type = self.surface_type_vb,
                                is_twisted = self.is_twisted_vb,
                                is_there_dihedral = self.is_there_dihedral_vb,
                                nb_eq_pts_airfoils_vb = self.nb_eq_pts_airfoils,
                                x_vb = self.x_vb,
                                y_vb = self.y_vb,
                                dtheta1 = self.dtheta1,
                                dtheta2 = self.dtheta2,
                                nb_max_rotations = self.nb_max_rotations,
                                surf2_type = self.surf2_type_vb,
                                surfaces = [self.wing.surface_ss, self.engines.surface],
                                is_there_vb = self.is_there_engine,
                                name = self.name_vb,
                                nb_chordwise_vortices = self.nb_chordwise_vortices_vb,
                                nb_spanwise_vortices = self.nb_spanwise_vortices_vb,
                                is_mirrored = False)

    @Part
    def engines(self):
        return Engine(fact_engine = self.fact_engine, AoA_engine = self.AoA_engine, x_vb = self.x_vb, y_vb = self.y_vb,
                          chords_vb = self.chords_vb, span_vb = self.span_vb,
                          crvs_surface_to_points_vb = self.vertical_bases.crvs_surface_to_points,
                          position_surf1_on_surf2_vb = self.vertical_bases.position_surf1_on_surf2,
                          nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                          is_there_engine=self.is_there_engine)

    # @Part
    # def vertical_tails(self):
    #     return VerticalTail(nb_eq_pts_rails = self.nb_eq_pts_rails_vt,
    #                           nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
    #                           nb_fixed_control_points = self.nb_fixed_control_points_vt,
    #                           sweeps_le = self.sweeps_le_vt,
    #                           span = self.span_vt,
    #                           x_fcp_except_last = self.x_fcp_vt_except_last,
    #                           y_fcp_except_last = self.y_fcp_vt_except_last,
    #                           chords = self.chords_vt,
    #                           nb_pts_dihedral = self.nb_pts_dihedral_vt,
    #                           angles_dihedral = self.angles_dihedral_vt,
    #                           spans_dihedral_except_last=self.spans_dihedral_vt_except_last,
    #                           nb_defined_airfoils = self.nb_defined_airfoils_vt,
    #                           airfoils_names = self.airfoils_names_vt,
    #                           y_defined_airfoils_except_last = self.y_defined_airfoils_except_last_vt,
    #                           cant_angles_airfoils = self.cant_angles_airfoils_vt,
    #                           dihedral_angles_airfoils = self.dihedral_angles_airfoils_vt,
    #                           heights_factors = self.heights_factors_vt,
    #                           nb_airfoils = self.nb_airfoils_vt,
    #                           surface_type = self.surface_type_vt,
    #                           is_twisted = self.is_twisted_vt,
    #                           is_there_dihedral = self.is_there_dihedral_vt,
    #                           nb_eq_pts_airfoils_vt = self.nb_eq_pts_airfoils,
    #                           x_vt=self.x_vt,
    #                           y_vt=self.y_vt,
    #                           dtheta1=self.dtheta1,
    #                           dtheta2 = self.dtheta2,
    #                           nb_max_rotations=self.nb_max_rotations,
    #                           surf2_type=self.surf2_type_vt,
    #                           surfaces=[self.wing.surface_ss,self.engines.surface],
    #                           is_there_vt=self.is_there_vt,
    #                           name = self.name_vt,
    #                           nb_chordwise_vortices = self.nb_chordwise_vortices_vt,
    #                           nb_spanwise_vortices = self.nb_spanwise_vortices_vt,
    #                           is_mirrored = False)


    @Attribute
    def avl_surfaces(self):  # this scans the product tree and collect all instances of the avl.Surface class
        return self.find_children(lambda o: isinstance(o, avl.Surface))

    @Part
    def avl_configuration(self):
        return avl.Configuration(name='bwb',
                                 reference_area=self.wing.planform_area,
                                 reference_span=self.wing.sspan*2,
                                 reference_chord=self.wing.mac,
                                 reference_point=self.position.point,
                                 surfaces=self.avl_surfaces,
                                 mach=self.Mach)

if __name__ == '__main__':
    from parapy.gui import display
    obj = BWB(label="BWB_2")
    display(obj)