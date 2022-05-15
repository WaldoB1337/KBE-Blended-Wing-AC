from math import sqrt
from parapy.core import *
from parapy.geom import *

from aircraft import Fuselage, Wing
import kbeutils.avl as avl


def calc_span(area, aspect_ratio):
    return sqrt(area * aspect_ratio)


class Aircraft(GeomBase):
    length = Input()  # fuselage length, m
    slenderness = Input()  # ratio fuselage length / max fuselage diameter

    wing_location = Input()  # longitudinal wing location, as % of fuselage length
    wing_area = Input()  # planform area of the total wing (right + left wing)
    wing_aspect_ratio = Input()  # square root (wing span**2/ wing area)
    wing_taper_ratio = Input()  # chord_tip / chord_root
    wing_le_sweep = Input()  # sweep angle measured at leading edge, in degrees
    twist = Input()
    wing_airfoil = Input()  # Name of the NACA airfoil (4 or 5 digits according to designation)
    elevator_hinge = Input()

    tail_area = Input()
    tail_aspect_ratio = Input()
    tail_taper_ratio = Input()
    tail_airfoil = Input()
    rudder_hinge = Input()

    mach = Input(0.2)

    @Input
    def tail_le_sweep(self):
        return self.wing_le_sweep + 5  # engineering rule: to allow teh tail a higher critical mach than for the wing

    @Part
    def fuselage(self):
        return Fuselage(length=self.length,
                        diameter=self.length/self.slenderness)

    @Attribute
    def wing_position(self):  # wing reference system. Same orientation as fuselage's, but with origin on wing LE @ root
        return self.position.translate('x', self.length*self.wing_location)

    @Attribute
    def wing_te_position(self):
        return self.wing.position.translate('x', self.wing.chord_root)

    @Attribute  # Intersection point between the fuselage crown curve and a plane passing through the wing TE edge
    def fuselage_crown_at_te(self):
        plane = Plane(reference=self.wing_te_position.point,
                      normal=self.position.Vx)
        return self.fuselage.crown_curve.intersection_point(plane)

    @Attribute
    def tail_te_position(self): # position of vertical tailplane
        rotated = self.wing_te_position.rotate90(self.position.Vx)
        vec = self.wing_te_position.point.vector_to(self.fuselage_crown_at_te)
        return rotated.translate(vec, vec.magnitude)

    @Attribute
    def avl_surfaces(self):  # this scans the product tree and collect all instances of the avl.Surface class
        return self.find_children(lambda o: isinstance(o, avl.Surface))

    @Part
    def wing(self):
        return Wing(name="wing",
                    span=calc_span(self.wing_area, self.wing_aspect_ratio),  # call to function calc_span
                    # span=sqrt(self.wing_area * self.wing_aspect_ratio), # worse alternative to the use of the function
                    aspect_ratio=self.wing_aspect_ratio,
                    taper_ratio=self.wing_taper_ratio,
                    le_sweep=self.wing_le_sweep,
                    twist=self.twist,
                    airfoil=self.wing_airfoil,

                    control_name='elevator',
                    control_hinge_loc=self.elevator_hinge,
                    duplicate_sign=1,

                    position=self.wing_position)

    @Part
    def tail(self):
        return Wing(name="tail",
                    span=calc_span(self.tail_area, self.tail_aspect_ratio),
                    aspect_ratio=self.tail_aspect_ratio,
                    taper_ratio=self.tail_taper_ratio,
                    le_sweep=self.tail_le_sweep,
                    twist=0,
                    airfoil=self.tail_airfoil,
                    control_name='rudder',
                    control_hinge_loc=self.rudder_hinge,

                    is_mirrored=False,
                    # position=self.wing_position.translate('x', -self.tail.chord_root)
                    position=self.tail_te_position.translate('x', -1*child.chord_root)
                    )

    @Part
    def avl_configuration(self):
        return avl.Configuration(name='aircraft',
                                 reference_area=self.wing.planform_area,
                                 reference_span=self.wing.span,
                                 reference_chord=self.wing.mac,
                                 reference_point=self.position.point,
                                 surfaces=self.avl_surfaces,
                                 mach=self.mach)


if __name__ == '__main__':
    from parapy.gui import display
    obj = Aircraft(label="fast primi plane",
                   length=16,
                   slenderness=12,
                   wing_location=0.42,
                   wing_area=50,
                   wing_aspect_ratio=2.5,
                   wing_taper_ratio=0.2,
                   wing_le_sweep=46,
                   twist=-5,
                   wing_airfoil='23008',
                   tail_area=8,
                   tail_aspect_ratio=2,
                   tail_taper_ratio=0.2,
                   tail_airfoil='0010'
                   )
    display(obj)






























