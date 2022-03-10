from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil  # note this is a different class than n exe_17_wing.py
from ref_frame import Frame


class LiftingSurface(LoftedSolid):  # note use of loftedSolid as superclass
    n_section = Input(4)
    airfoils = Input(["whitcomb", "whitcomb", "whitcomb", "whitcomb", "whitcomb"])
    s_chords = Input([20.0, 8.0, 5.0, 5.0, 5.0])
    t_factors = Input([2.0, 1.0, 1.0, 1.0, 1.0])

    s_spans = Input([8.0, 13.0, 1.0, 3])
    sweep_locs = Input([0.2, 0.2, 0.2, 0.2])  # location of sweep
    sweep_angles = Input([40, 50, 40, 40])  # sweep angle
    dihedral_angles = Input([0, -10, 45, 45])
    twist_locs = Input([0.2, 0.2, 0.2, 0.2, 0.2])
    twist_angles = Input([0, 0, 0, 0, 0])

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
                P = rotate(rotate(translate(lst_positions[-1], 'x', self.sweep_locs[i-1]*self.s_chords[i-1] - (self.s_chords[i]  - (
                                       tan(radians(self.sweep_angles[i-1])) * self.s_spans[i-1] + (1-self.sweep_locs[i-1])*self.s_chords[i])), 'y', self.s_spans[i-1], 'z', self.s_spans[i-1] * tan(radians(self.dihedral_angles[i-1]))),
               Vector(self.twist_locs[i], 1, 0), self.twist_angles[i],
               deg=True),Vector(self.twist_locs[i], 0, 0),self.dihedral_angles[i-1],deg=True)
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
        return LoftedSolid(profiles=self.profiles,
                           color='red',
                           hidden=not(__name__ == '__main__'))


if __name__ == '__main__':
    from parapy.gui import display
    obj = LiftingSurface(label="lifting surface",
                         mesh_deflection=0.0001
                        )
    display(obj)
