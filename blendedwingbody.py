import os
from math import radians
from parapy.exchange import STEPWriter
from parapy.geom import *
from parapy.core import *
from liftingsurface import LiftingSurface  # note this is not teh same wing class defined in tutorial 5
from fuselage import Fuselage
from ref_frame import Frame

DIR = os.path.dirname(__file__)


class BlendedWingBody(GeomBase):
    n_section = Input(4)
    airfoils = Input(["whitcomb", "whitcomb", "whitcomb", "whitcomb", "whitcomb"])
    s_chords = Input([6., 6., 6., 6., 6.])
    t_factors = Input([1., 1., 1., 1., 1.])

    s_spans = Input([1., 1., 1., 1.])
    sweep_locs = Input([0.2, 0.2, 0.2, 0.2])  # location of sweep
    sweep_angles = Input([30, 30, 30, 30])  # sweep angle
    dihedral_angles = Input([0, 0, 0, 0])
    twist_locs = Input([0.2, 0.2, 0.2, 0.2, 0.2])
    twist_angles = Input([0, 0, 0, 0, 20])

    @Part
    def aircraft_frame(self):
        return Frame(pos=self.position)  # this helps visualizing the wing local reference frame


    @Part
    def right_wing(self):
        return LiftingSurface(pass_down="n_section, airfoils, s_chords, t_factors, s_spans,"
                                        "sweep_locs, sweep_angles, twist_locs, twist_angles"
                                        "dihedral_angles",
                              position=self.position,
                              mesh_deflection=0.0001)

    @Part
    def left_wing(self):
        return MirroredShape(shape_in=self.right_wing,
                             reference_point=self.position,
                             # Two vectors to define the mirror plane
                             vector1=self.position.Vz,
                             vector2=self.position.Vx,
                             mesh_deflection=0.0001)

if __name__ == '__main__':
    from parapy.gui import display
    obj = BlendedWingBody(label="BWB")
    display(obj)