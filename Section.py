from math import radians, tan
from parapy.geom import *
from parapy.core import *
from ref_frame import Frame


class Wing (Base,GeomBase):

    n_sections = Input(4)
    s_c_root = Input([3,4,5,6])  # section root chord
    s_c_tip = Input([4,5,6,7])  # section tip chord
    s_span = Input([8,9,2,4])  # section span
    sweep_locs = Input([0.2,0.2,0.2,0.2])  # location of sweep
    sweep_angles = Input([30,30,30,30])  # sweep angle
    dihedral_angles = Input([30,30,30,30])
    twist_locs = Input([0.2,0.2,0.2,0.2])
    twist_angles = Input([30,30,30,30])
    #: :type: str
    airfoils_root = Input(['withcomb','withcomb','withcomb','withcomb'])
    airfoils_tip = Input(['withcomb','withcomb','withcomb','withcomb'])

    @Part
    def wing_frame(self):
        return Frame(pos=self.position)  # this helps visualizing the aircraft reference frame, /
        # which, in this case, is the same as the global reference frame XOY)

    @Part
    def sections(self):
        return Section(quantify=self.n_sections,
                       c_root=self.s_c_root[child.index % len(self.s_c_root)],
                       c_tip=self.s_c_tip[child.index % len(self.s_c_tip)],
                       span=self.s_span[child.index % len(self.s_span)],
                       sweep_loc=self.sweep_locs[child.index % len(self.sweep_locs)],
                       sweep_angle=self.sweep_angles[child.index % len(self.sweep_angles)],
                       dihedral_angle=self.sweep_angles[child.index % len(self.dihedral_angles)],
                       twist_loc=self.twist_locs[child.index % len(self.twist_locs)],
                       twist_angle=self.twist_angles[child.index % len(self.twist_angles)],
                       airfoil_root=self.airfoils_root[child.index % len(self.airfoils_root)],
                       airfoil_tip=self.airfoils_tip[child.index % len(self.airfoils_root)])


class Section(GeomBase):
    """Basic wing geometry: a loft between root and tip airfoil"""
    c_root = Input()  # section root chord
    c_tip = Input()  # section tip chord
    span = Input()  # section span
    sweep_loc = Input()  # location of sweep
    sweep_angle = Input()  # sweep angle
    dihedral_angle = Input()
    twist_loc = Input()
    twist_angle = Input()
    #: :type: str
    airfoil_root = Input()
    airfoil_tip = Input()

    @Attribute
    def pts_root(self):
        """ Extract airfoil coordinates from a data file and create a list of 3D points"""
        with open(self.airfoils_root[child.index]+'.dat', 'r') as f:
            points = []
            for line in f:
                x, y = line.split(' ', 1)  # separator = " "; max number of split = 1
                # Convert the strings to numbers and make 3D points for the FittedCurve class
                points.append(Point(float(x), float(y)))
        return points

    def pts_tip(self):
        """ Extract airfoil coordinates from a data file and create a list of 3D points"""
        with open(self.airfoils_tip[child.index]+'.dat', 'r') as f:
            points = []
            for line in f:
                x, y = line.split(' ', 1)
                # separator = " "; max number of split = 1
                # Convert the strings to numbers and make 3D points for the FittedCurve class
                points.append(Point(float(x), float(y)))
        return points

    @Part
    def section_frame(self):
        return Frame(pos=self.position)  # this helps visualizing the wing local reference frame

    @Part
    def airfoil_from_3D_points_root(self):  # this curve is on the X-Y plane, with TE = (1, 0, 0) and LE = (0, 0, 0)
        return FittedCurve(self.pts_root,
                           mesh_deflection=0.0001)

    @Part
    def airfoil_from_3D_points_tip(self):  # this curve is on the X-Y plane, with TE = (1, 0, 0) and LE = (0, 0, 0)
        return FittedCurve(self.pts_tip,
                           mesh_deflection=0.0001)


    @Part  # TransformedCurve is making a carbon copy of the fitted curve, which can be moved (rotated and translated) /
    # from one position to another. /
    # In this case we want to position the fitted curve copy in the x-z plane of the wing reference system, with its /
    # TE in the origin (location) of this reference system. This requires a rotation and a few translations.
    def root_section_unscaled(self):
        return TransformedCurve(self.airfoil_from_3D_points_root,  # the curve to be repositioned
                                rotate(translate(XOY, 'x', 1), 'x', -90, deg=True),  # from_position
                                self.position,  # to_position. The wing relative reference system
                                hidden=False)

    @Part  # for the wing tip we use the same type of airfoil used for the wing root. We use again TransformedCurve
    def tip_section_unscaled(self):
        return TransformedCurve(self.airfoil_from_3D_points_root,  # the curve to be repositioned
    translate(self.root_section_unscaled.position, 'x', -self.sweep_loc),  # from_position
    rotate(translate(translate(self.root_section_unscaled.position,  # to_position, i.e. the wing tip section
                               'y', self.s_span,
                               'x', self.s_span * tan(radians(self.sweep_angle))), 'z',
                     self.s_span * tan(radians(self.dihedral_angle))), Vector(-self.twist_loc, 1, 0), self.twist_angle,
           deg=True),  # the sweep is applied
    hidden = True)

    '''curve_in = self.airfoil_from_3D_points_tip,
    displacement = rotate(
        translate(translate(self.root_section_unscaled.position,  # to_position, i.e. the wing tip section
                            'y', self.s_span,
                            'x', self.s_span * tan(radians(self.sweep_angle))), 'z',
                  self.s_span * tan(radians(self.dihedral_angle))), Vector(-self.twist_loc, 1, 0), self.twist_angle,
        deg=True)'''

    @Part
    def root_section(self):  # the ScaledCurve primitive allows scaling a given curve. Here it is used to scale /
        # the unit chord airfoil generated from the .dat file according to their actual chord length
        return ScaledCurve(curve_in=self.root_section_unscaled,
                           reference_point=self.root_section_unscaled.start,  # this point (the curve TE in this case) /
                           # is kept fixed during scaling
                           factor=self.s_c_root,  # uniform scaling
                           mesh_deflection=0.0001)

    @Part
    def tip_section(self):
        return ScaledCurve(self.tip_section_unscaled, self.tip_section_unscaled.start, self.s_c_tip)

    @Part
    def wing_loft_surf(self):  # generate a surface
        return LoftedSurface([self.root_section, self.tip_section],
                           mesh_deflection=0.0001)

    @Part
    def wing_loft_solid(self):  # generate a solid
        return LoftedSolid([self.root_section, self.tip_section],
                           mesh_deflection=0.0001)


if __name__ == '__main__':
    from parapy.gui import display

    obj = Aircraft(label="aircraft")
    display(obj)
