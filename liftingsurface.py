from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil  # note this is a different class than n exe_17_wing.py
from ref_frame import Frame


class LiftingSurface(LoftedSolid):  # note use of loftedSolid as superclass
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

    @Attribute  # required input for the superclass LoftedSolid
    def profiles(self):
        return [A for A in self.airfoils_pos]

    @Part
    def lifting_surf_frame(self):  # to visualize the given lifting surface reference frame
        return Frame(pos=self.position,
                     hidden=False)

    @Attribute
    def airfoils_pos(self):
        lst_airfoils = []
        lst_positions = []
        for i in range(self.n_section+1):
            if i==0:
                P = rotate(self.position,Vector(self.twist_locs[i], 1, 0), self.twist_angles[i],
               deg=True)
                A = Airfoil(airfoil_name=self.airfoils[i],
                       chord=self.s_chords[i],
                       thickness_factor=self.t_factors[i],
                       mesh_deflection=0.0001)
                lst_airfoils.append(A)
                lst_positions.append(P)
            else:
                P = rotate(translate(lst_positions[-1], 'x', self.sweep_locs[i-1]*self.s_chords[i-1] - (self.s_chords[i]  - (
                                       tan(radians(self.sweep_angles[i-1])) * self.s_spans[i-1] + (1-self.sweep_locs[i-1])*self.s_chords[i])), 'y', self.s_spans[i-1], 'z', self.s_spans[i-1] * tan(radians(self.dihedral_angles[i-1]))),
               Vector(self.twist_locs[i], 1, 0), self.twist_angles[i],
               deg=True)
                A = Airfoil(airfoil_name=self.airfoils[i],
                       chord=self.s_chords[i],
                       thickness_factor=self.t_factors[i],
                       position=P,
                       mesh_deflection=0.0001)
                lst_airfoils.append(A)
                lst_positions.append(P)

        return lst_airfoils

    @Part
    def lofted_surf(self):
        return LoftedSurface(profiles=self.profiles,
                             hidden=not(__name__ == '__main__'))


if __name__ == '__main__':
    from parapy.gui import display
    obj = LiftingSurface(label="lifting surface",
                         mesh_deflection=0.0001
                        )
    display(obj)
