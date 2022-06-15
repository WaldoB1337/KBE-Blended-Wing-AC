from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from copy import deepcopy
from parapy.gui import display

import numpy as np
from math import *

from surface_2 import Surface_BWB
from engine_1 import Engine
from vertical_tail_1 import VerticalTail
from BWB_CentralFuselage import CentralFuselage
from BWB_Fuel import FuelTank


def is_point_inside_surf(airfoils_pts, pt):
    """
    Returns True (False) if the point 'pt' is inside (not) the surface fitted with the points 'airfoils_pts'. The
    surface that will be considered is the left surface of the wing 'surface_ss'. It is an opened surface at the root
    and tip. Note that a method .is_point_inside() is available in Parapy's functions, but it does not work with opened
    surfaces.
    The airfoil which is the closest to the point (along the y-axis) is firstly found. Then, the two closest x-ordinates
    (upper and lower surfaces) of this airfoil to the point are identified. If the point is (not) between the z_lower
    and z_upper ordinates associated to this x-ordinate of the considered airfoil, the point is considered (not) to be
    enclosed by the surface and the function returns True (False).

    This function will be very useful to ensure to place the 'Central_Fuselage' parts (Cargo, Cabin and Fuel tanks) in
    the external surface of the fuselage.

    """
    y_of_airfoil = [airfoils_pts[i][0][1] for i in range(len(airfoils_pts))]
    idx_airfoil = (np.abs(np.asarray(y_of_airfoil) - pt[1])).argmin()
    x_of_airfoil = [airfoils_pts[idx_airfoil][i][0] for i in range(len(airfoils_pts[idx_airfoil]))]
    idxs_x_in_airfoil = np.argpartition(np.abs(np.asarray(x_of_airfoil) - pt[0]), 2)
    # print('z_low=',airfoils_pts[idx_airfoil][idxs_x_in_airfoil[1]][2],'z_up=',airfoils_pts[idx_airfoil][idxs_x_in_airfoil[0]][2])
    z_1 = airfoils_pts[idx_airfoil][idxs_x_in_airfoil[1]][2]
    z_2 = airfoils_pts[idx_airfoil][idxs_x_in_airfoil[0]][2]
    if z_1 <= z_2:
        return airfoils_pts[idx_airfoil][idxs_x_in_airfoil[1]][2] <= pt[2] <= \
               airfoils_pts[idx_airfoil][idxs_x_in_airfoil[0]][2], airfoils_pts[idx_airfoil]
    return airfoils_pts[idx_airfoil][idxs_x_in_airfoil[0]][2] <= pt[2] <= \
           airfoils_pts[idx_airfoil][idxs_x_in_airfoil[1]][2], airfoils_pts[idx_airfoil]


class BWB(GeomBase):
    nb_eq_pts_rails = Input(200)
    nb_eq_pts_airfoils = Input(200)

    nb_fixed_control_points_wing = Input(5)
    sweeps_le_wing = Input([25, 50, 55, 40, 80])
    sweeps_te_wing = Input([-40, -40, 0, 30, 55])
    sspan_wing = Input(7.5)
    x_fcp_wing = Input([0, -0.75, -2.5, -6, -7])
    y_fcp_wing_except_last = Input([0, 1, 2.5, 7])
    chords_wing = Input([7, 7.1, 3.5, 0.8, 0.3])
    rho_le_wing = Input([.5, .5, .5, .5, .5])
    rho_te_wing = Input([.5, .5, .5, .5, .5])

    nb_pts_twist_wing = Input(3)
    angles_twist_wing = Input([0, 0, 0])  # degrees of twist
    y_twist_wing_except_last = Input([0, 7])  # location of the definition of twist
    nb_pts_twist_axis_wing = Input(3)
    x_twist_axis_wing = Input([0, -1, -1])
    y_twist_axis_wing_except_last = Input([0, 3])  # 1st and last points must be at the root and tip

    nb_pts_dihedral_wing = Input(5)
    angles_dihedral_wing = Input([5, 5, 5, 5, 89])  # angle defined on the right of the section
    spans_dihedral_wing_except_last = Input([2.75, 2.75, 1.4])

    weights_dihedral_wing = Input([.5, .5, .5, .5, .5, .5, .5, .5, .5, 1000, 1000, 1000, 1000])

    nb_defined_airfoils_wing = Input(6)
    airfoils_names_wing = Input(["23012", "23012", "23012", "23012", "0003", "0003"])
    y_defined_airfoils_wing_except_last = Input([0, 3, 4, 5, 7])  # 1st and last correspond to root and tip
    cant_angles_airfoils_wing = Input([0, 0, 0, 0, 0, 0])  # degree
    dihedral_angles_airfoils_wing = Input([0, 0, 0, 0, 0, 89])
    heights_factors_wing = Input([1.5, 1.2, 1, 1, 1, 1])

    nb_airfoils_ss_wing = Input(50)

    surface_type_wing = Input('wing')
    is_twisted_wing = Input(True)
    is_there_dihedral_wing = Input(True)

    name_wing = Input('Main_wing')
    nb_chordwise_vortices_wing = Input(30)
    nb_spanwise_vortices_wing = Input(50)
    is_mirrored_wing = Input(True)
    nb_sections_AVL_wing = Input(5)

    # ----------------------------------------------

    nb_eq_pts_rails_vt = Input(2)
    nb_fixed_control_points_vt = Input(2)
    sweeps_le_vt = Input([50, 50])
    span_vt = Input(1)
    x_fcp_vt_except_last = Input([0])
    y_fcp_vt_except_last = Input([0])
    chords_vt = Input([1.2, 0.4])

    nb_pts_dihedral_vt = Input(2)
    angles_dihedral_vt = Input([85, 85])  # angle defined on the right of the section
    spans_dihedral_vt_except_last = Input([])

    nb_pts_twist_vt = Input(2)
    angles_twist_vt = Input([30, 30])  # degrees of twist
    y_twist_vt_except_last = Input([0])  # location of the definition of twist
    nb_pts_twist_axis_vt = Input(2)
    x_twist_axis_vt = Input([-np.tan(10)*0.2, -np.tan(10)*0.2-0.5])
    y_twist_axis_vt_except_last = Input([0.2, 0.2])

    nb_defined_airfoils_vt = Input(2)
    airfoils_names_vt = Input(["24012", "24012"])
    y_defined_airfoils_except_last_vt = Input([0])  # 1st and last correspond to root and tip
    cant_angles_airfoils_vt = Input([0, 0])  # degree
    dihedral_angles_airfoils_vt = Input([90, 90])
    heights_factors_vt = Input([1, 1])

    nb_airfoils_vt = Input(2)

    surface_type_vt = Input('vt')
    is_twisted_vt = Input(True)
    is_there_dihedral_vt = Input(True)

    x_vt = Input(-6.75)
    y_vt = Input(1.4)

    dtheta1 = Input(0.1)
    dtheta2 = Input(0.01)
    nb_max_rotations = Input(200)

    surf2_type_vt = Input('wing')

    is_there_vt = Input(False)

    name_vt = Input('tail')
    nb_chordwise_vortices_vt = Input(30)
    nb_spanwise_vortices_vt = Input(50)

    # ---------------------------------------------------

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
    y_vb = Input(0.5)

    surf2_type_vb = Input('wing')

    name_vb = Input('base')
    nb_chordwise_vortices_vb = Input(30)
    nb_spanwise_vortices_vb = Input(50)

    # --------------------------------

    is_there_engine = Input(True)
    fact_engine = Input(110)
    AoA_engine = Input(20)

    # --------------------------------

    Mach = Input(0.11)

    # --------------------------------

    N_pax = Input(400)
    N_col = Input(3)
    nb_max_iter_adjust_fuselage = Input(30)
    scaling_external_fus_increment = Input(1.1)
    tank_c_root = Input(15)
    sweep_increment = Input(3)
    chord_factor_increment = Input(4 / 5)

    @Part
    def central_fuselage(self):
        """
        Returns the central fuselage defined with the given number of passengers 'N_pax' and the number of economy
        columns 'N_col'.

        """
        return CentralFuselage(N_pax = self.N_pax, N_col = self.N_col)

    @Part
    def wing(self):
        """
        Returns the wing of the original BWB defined with the data given by the user. This wing will endorse the AVL
        calculation even if it is not the scaled one at the central fuselage size. The data given by the user and
        especially the lengths are given for reference and to create a harmonic fuselage. The real lengths of the BWB
        are scaled afterward to include the internal geometries defined by the user.

        """
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
        """
        Returns the vertical base of the original BWB defined with the data given by the user.
        """
        return VerticalTail(nb_eq_pts_rails = self.nb_eq_pts_rails_vb,
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
                            nb_eq_pts_airfoils_vt = self.nb_eq_pts_airfoils,
                            x_vt = self.x_vb,
                            y_vt = self.y_vb,
                            dtheta1 = self.dtheta1,
                            dtheta2 = self.dtheta2,
                            nb_max_rotations = self.nb_max_rotations,
                            surf2_type = self.surf2_type_vb,
                            surfaces = [self.wing.surface_ss, self.engines.surface],
                            is_there_vt = self.is_there_engine,
                            name = self.name_vb,
                            nb_chordwise_vortices = self.nb_chordwise_vortices_vb,
                            nb_spanwise_vortices = self.nb_spanwise_vortices_vb,
                            is_mirrored = False)

    @Part
    def engines(self):
        """
        Returns the engine of the original BWB defined with the data given by the user.
        """
        return Engine(fact_engine = self.fact_engine, AoA_engine = self.AoA_engine, x_vb = self.x_vb, y_vb = self.y_vb,
                      chords_vb = self.chords_vb, span_vb = self.span_vb,
                      crvs_surface_to_points_vb = self.vertical_bases.crvs_surface_to_points,
                      position_surf1_on_surf2_vb = self.vertical_bases.position_surf1_on_surf2,
                      nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                      is_there_engine = self.is_there_engine)

    @Part
    def vertical_tails(self):
        return VerticalTail(nb_eq_pts_rails = self.nb_eq_pts_rails_vt,
                              nb_eq_pts_airfoils = self.nb_eq_pts_airfoils,
                              nb_fixed_control_points = self.nb_fixed_control_points_vt,
                              sweeps_le = self.sweeps_le_vt,
                              span = self.span_vt,
                              x_fcp_except_last = self.x_fcp_vt_except_last,
                              y_fcp_except_last = self.y_fcp_vt_except_last,
                              chords = self.chords_vt,
                              nb_pts_dihedral = self.nb_pts_dihedral_vt,
                              angles_dihedral = self.angles_dihedral_vt,
                              spans_dihedral_except_last=self.spans_dihedral_vt_except_last,
                              nb_pts_twist = self.nb_pts_twist_vt,
                              angles_twist = self.angles_twist_vt,
                              y_twist_except_last = self.y_twist_vt_except_last,
                              nb_pts_twist_axis = self.nb_pts_twist_axis_vt,
                              x_twist_axis = self.x_twist_axis_vt,
                              y_twist_axis_except_last = self.y_twist_axis_vt_except_last,
                              nb_defined_airfoils = self.nb_defined_airfoils_vt,
                              airfoils_names = self.airfoils_names_vt,
                              y_defined_airfoils_except_last = self.y_defined_airfoils_except_last_vt,
                              cant_angles_airfoils = self.cant_angles_airfoils_vt,
                              dihedral_angles_airfoils = self.dihedral_angles_airfoils_vt,
                              heights_factors = self.heights_factors_vt,
                              nb_airfoils = self.nb_airfoils_vt,
                              surface_type = self.surface_type_vt,
                              is_twisted = self.is_twisted_vt,
                              is_there_dihedral = self.is_there_dihedral_vt,
                              nb_eq_pts_airfoils_vt = self.nb_eq_pts_airfoils,
                              x_vt=self.x_vt,
                              y_vt=self.y_vt,
                              dtheta1=self.dtheta1,
                              dtheta2 = self.dtheta2,
                              nb_max_rotations=self.nb_max_rotations,
                              surf2_type=self.surf2_type_vt,
                              surfaces=[self.wing.surface_ss,self.engines.surface],
                              is_there_vt=self.is_there_vt,
                              name = self.name_vt,
                              nb_chordwise_vortices = self.nb_chordwise_vortices_vt,
                              nb_spanwise_vortices = self.nb_spanwise_vortices_vt,
                              is_mirrored = False)

    @Attribute
    def adjust_wing_central_fuselage(self):
        """
        Returns the adjusted fuselage curves points defining the external wing scaled with the interior central fuselage.
        The external surface is initiated with the original one. Then, the external surface is scaled with an increment
        given by the user until all the corners points of the left side of the internal geometries are enclosed by the
        left external surface of the wing. The total scaling factor is also returned by the function. The scaling
        increment can be modified by the user.

        """
        wing = Surface_BWB(nb_eq_pts_rails = self.nb_eq_pts_rails,
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

        is_inside = False
        total_scale_factor = 1
        iter = 1
        interior_corners = self.central_fuselage.symm_corners
        crvs = []
        pts = []
        airfoils_crvs = wing.crvs_surface_all_airfoils
        vect_displacement = airfoils_crvs[0].midpoint.translate("x",
                                                                self.central_fuselage.cabin.l_cabin * 1.25 + self.central_fuselage.cabin.l_cockpit / 2).vector_from(
            airfoils_crvs[0].midpoint)

        airfoils_crvs = [
            TranslatedCurve(airfoils_crvs[i], vect_displacement) for i in
            range(self.nb_airfoils_ss_wing)]
        iter = 1
        while not is_inside and iter < self.nb_max_iter_adjust_fuselage:
            airfoils_crvs = [
                ScaledCurve(airfoils_crvs[i], airfoils_crvs[0].midpoint, self.scaling_external_fus_increment) for i in
                range(self.nb_airfoils_ss_wing)]
            pts_surface_lst = []
            crvs.append(airfoils_crvs[0])
            for i in range(self.nb_airfoils_ss_wing):
                pts_airfoil_lst = airfoils_crvs[i].equispaced_points(self.nb_eq_pts_airfoils)
                pts_surface_lst.append(pts_airfoil_lst)
            pts.append(pts_surface_lst)
            total_scale_factor *= self.scaling_external_fus_increment
            is_inside = True
            for i in range(ceil(len(interior_corners[0]))):
                statement = is_point_inside_surf(pts_surface_lst,
                                                 interior_corners[0][i])[0]
                is_inside *= statement
            iter += 1
        if iter < self.nb_max_rotations:
            return pts_surface_lst, total_scale_factor

    @Attribute
    def adjust_fuel_tank(self):
        """
        Returns the tank root and tip chord as well as the sweep angle of the tanks in order to ensure that the tank is
        contained in the external surface of the wing.

        The tank is initialized with the pre-filled or parameters given by the user and defining the tank geometry.
        The tip of the tank is firstly modified to be contained at the 80% span of the wing location. The leading and
        trailing edge of the tank tip are respectively defined at 80% and 45% of the tank tip chord of the external
        airfoil. The tank tip height is modified according to the height of the airfoil at 80% span of the wing. The
        dihedral angle of the tank is modified so that the tank tip is aligned with the 80% span of the wing airfoil.

        Then, the tank tip is place inside the wing by playing with the sweep angle until the corners of the tank tip
        are all placed in the wing surface.

        Finally, the tank root chord is modified to ensure that the tank leading and trailing edge are inside the wing
        all along the span of the tank. Points are defined along these sides edges. If all these points are inside the
        wing surface, the tank is finally considered as enclosed by the wing and the configuration is retained.

        All the increments (sweep, tank root chord) can be modified by the user.

        """
        pts_surface_lst = self.adjust_wing_central_fuselage[0]
        total_scale_factor = self.adjust_wing_central_fuselage[1]
        airfoil_tank_tip = pts_surface_lst[
            int(self.sspan_wing * 0.8 / self.sspan_wing * self.nb_airfoils_ss_wing)]

        tank_tip_x_le = airfoil_tank_tip[0][0] + 0.8 * (
                airfoil_tank_tip[int(self.nb_eq_pts_airfoils / 2)][0] - airfoil_tank_tip[0][0])
        tank_tip_x_te = airfoil_tank_tip[0][0] + 0.45 * (
                airfoil_tank_tip[int(self.nb_eq_pts_airfoils / 2)][0] - airfoil_tank_tip[0][0])
        tank_c_tip = tank_tip_x_le - tank_tip_x_te

        x_of_airfoil = [airfoil_tank_tip[i][0] for i in range(len(airfoil_tank_tip))]
        idxs_x_in_airfoil_le = np.argpartition(np.abs(np.asarray(x_of_airfoil) - tank_tip_x_le), 2)
        idxs_x_in_airfoil_te = np.argpartition(np.abs(np.asarray(x_of_airfoil) - tank_tip_x_te), 2)
        z_le_1 = airfoil_tank_tip[idxs_x_in_airfoil_le[1]][2]
        z_le_2 = airfoil_tank_tip[idxs_x_in_airfoil_le[0]][2]
        z_te_1 = airfoil_tank_tip[idxs_x_in_airfoil_te[1]][2]
        z_te_2 = airfoil_tank_tip[idxs_x_in_airfoil_te[0]][2]
        if np.abs(z_le_1 - z_le_2) <= np.abs(z_te_1 - z_te_2):
            tank_tip_height = np.abs(z_le_1 - z_le_2)
        else:
            tank_tip_height = np.abs(z_te_1 - z_te_2)

        tank = FuelTank(w_cabin = self.central_fuselage.cabin.w_cabin,
                        wingspan = self.sspan_wing * 2 * total_scale_factor, taper = tank_c_tip / self.tank_c_root,
                        loc = self.central_fuselage.cog_interior.rotate("z", -90, deg = True),
                        h_tank_tip = tank_tip_height)

        dihedral_height = airfoil_tank_tip[int(0.33 * self.nb_eq_pts_airfoils)][2] + (
                airfoil_tank_tip[int(0.66 * self.nb_eq_pts_airfoils)][2] -
                airfoil_tank_tip[int(0.33 * self.nb_eq_pts_airfoils)][2]) / 2 - tank.fuel_tank1.profile1.position[2]

        sweep = self.sweep_increment
        are_inside = False
        iter = 1
        while not are_inside and iter < self.nb_max_iter_adjust_fuselage:
            sweep -= self.sweep_increment
            tank = FuelTank(w_cabin = self.central_fuselage.cabin.w_cabin,
                            wingspan = self.sspan_wing * 2 * total_scale_factor, taper = tank_c_tip / 15, sweep = sweep,
                            loc = self.central_fuselage.cog_interior.rotate("z", -90, deg = True),
                            h_tank_tip = tank_tip_height, dihedral_height = dihedral_height)
            are_inside = is_point_inside_surf(pts_surface_lst, tank.corners[0][-1])[0] * is_point_inside_surf(
                pts_surface_lst, tank.corners[0][-2])[0] * is_point_inside_surf(
                pts_surface_lst, tank.corners[0][-3])[0] * is_point_inside_surf(pts_surface_lst, tank.corners[0][-4])[0]
            iter += 1
        are_inside = False
        iter = 1
        tank_c_root = tank.tank_c_root / self.chord_factor_increment
        while not are_inside and iter < self.nb_max_iter_adjust_fuselage:
            are_inside = True
            tank_c_root = tank_c_root * self.chord_factor_increment
            tank = FuelTank(w_cabin = self.central_fuselage.cabin.w_cabin,
                            wingspan = self.sspan_wing * 2 * total_scale_factor, taper = tank_c_tip / tank_c_root,
                            sweep = sweep, tank_c_root = tank_c_root,
                            loc = self.central_fuselage.cog_interior.rotate("z", -90, deg = True),
                            h_tank_tip = tank_tip_height, dihedral_height = dihedral_height)
            for i in range(int(len(tank.tank_border) / 2)):
                statement = is_point_inside_surf(pts_surface_lst, tank.tank_border[i].position)[0]
                are_inside *= statement
            iter +=1
        if iter < self.nb_max_iter_adjust_fuselage:
            return tank
        return "The tank cannot fit in the external surface because part of the side edges are out of it"

    @Attribute
    def mirror_points(self):
        pts_surface_ss_lst = self.adjust_wing_central_fuselage[0]
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
    def wing_surface_scaled(self):
        """
        Returns the surface of the wing scaled with the internal central fuselage.
        """
        return Fused(FittedSurface(points = self.adjust_wing_central_fuselage[0], min_degree = 8)
                     , FittedSurface(points = self.mirror_points, min_degree = 8))

    @Part
    def tanks_scaled(self):
        """
        Returns the corresponding tank of the BWB scaled with the internal central fuselage.
        """
        return FuelTank(wingspan = self.sspan_wing * 2, taper = self.adjust_fuel_tank[1] / self.adjust_fuel_tank[0],
                        sweep = self.adjust_fuel_tank[2], tank_c_root = self.adjust_fuel_tank[0])

    @Attribute
    def avl_surfaces(self):
        """
        Returns the AVL surfaces of the original wing to be analyzed by AVL.
        """
        return self.find_children(lambda o: isinstance(o, avl.Surface))

    @Part
    def avl_configuration(self):
        """
        Returns the AVL configuration of the original wing to be analyzed by AVL.
        """
        return avl.Configuration(name = 'bwb',
                                 reference_area = self.wing.planform_area,
                                 reference_span = self.wing.sspan * 2,
                                 reference_chord = self.wing.mac,
                                 reference_point = self.position.point,
                                 surfaces = self.avl_surfaces,
                                 mach = self.Mach)


if __name__ == '__main__':
    from parapy.gui import display

    obj = BWB(label = "BWB_2")
    display(obj)
