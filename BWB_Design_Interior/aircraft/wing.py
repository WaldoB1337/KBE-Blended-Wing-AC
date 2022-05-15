from math import radians, tan

from parapy.core import *
from parapy.geom import *
import kbeutils.avl as avl
from aircraft import Section


# from aircraft import Section


class Wing(GeomBase):
    name = Input()
    span = Input()
    aspect_ratio = Input()
    taper_ratio = Input()
    le_sweep = Input()
    twist = Input()
    airfoil = Input()
    is_mirrored = Input(True)

    # Controls
    control_name = Input(None)
    control_hinge_loc = Input(None)
    duplicate_sign = Input(1)

    @Attribute
    def planform_area(self):
        return self.span**2 / self.aspect_ratio

    @Attribute
    def half_span(self):
        return self.span /2 if self.is_mirrored else self.span

    @Attribute
    def chords(self):
        root = ((2 * self.planform_area) /
                (self.span * (1 + self.taper_ratio)))
        tip = root * self.taper_ratio
        return root, tip

    @Attribute
    def chord_root(self):
        return self.chords[0]

    @Attribute
    def mac(self):
        return 2/3*self.chord_root*(1+self.taper_ratio+self.taper_ratio**2)/(1+self.taper_ratio)

    @Attribute
    def section_positions(self):
        sweep = radians(self.le_sweep)
        root = self.position
        tip = rotate(self.position.translate('x', self.half_span * tan(sweep),'y', self.half_span),
                     'y', self.twist, deg=True)
        return root, tip

    @Part
    def sections(self):
        return Section(quantify=len(self.chords),
                       airfoil_name=self.airfoil,
                       chord=self.chords[child.index],
                       position=self.section_positions[child.index],

                       control_name=self.control_name,
                       control_hinge_loc=self.control_hinge_loc,
                       duplicate_sign=self.duplicate_sign
                       )

    @Part
    def surface(self):
        return LoftedShell(profiles=[section.curve for section in self.sections],
                           mesh_deflection=0.0001)

    @Part
    def mirrored(self):
        return MirroredSurface(surface_in=self.surface.faces[0],
                               reference_point=self.position.point,
                               vector1=self.position.Vx,
                               vector2=self.position.Vz,
                               suppress=not self.is_mirrored,
                               mesh_deflection=0.0001)

    @Part
    def avl_surface(self):
        return avl.Surface(name=self.name,
                           n_chordwise=12,
                           chord_spacing=avl.Spacing.cosine,
                           n_spanwise=20,
                           span_spacing=avl.Spacing.cosine,
                           y_duplicate=self.position.point[1] if self.is_mirrored else None,
                           sections=[section.avl_section
                                     for section in self.sections])


if __name__ == '__main__':
    from parapy.gui import display
    obj = Wing(name='wing',
               span=20,
               aspect_ratio=6,
               taper_ratio=0.2,
               le_sweep=30,
               twist=-5,
               airfoil='2024',
               is_mirrored=False)
    display(obj)




























