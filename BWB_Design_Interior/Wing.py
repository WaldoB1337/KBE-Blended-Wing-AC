from contextlib import suppress
from math import radians, tan

from parapy.core import Input, Attribute, Part, child
from parapy.geom import GeomBase, LoftedShell, MirroredSurface

from Section import Section ## Importing class from Section.py

class Wing(GeomBase):

    name = Input()
    span = Input()
    aspect_ratio = Input()
    taper_ratio = Input()
    le_sweep = Input()
    airfoil = Input()

    is_mirrored = Input(True)

    @Attribute
    def planform_area(self):
        return self.span**2/self.aspect_ratio

    @Attribute
    def half_span(self):
        return self.span / 2 if self.is_mirrored else self.span

    @Attribute
    def chords(self):
        root = ((2 * self.planform_area)/
                (self.span * (1 + self.taper_ratio)))
        tip = root * self.taper_ratio
        return root, tip
    
    @Attribute
    def section_positions(self):
        sweep = radians(self.le_sweep)
        root = self.position
        tip = self.position.translate("x", self.half_span * tan(sweep),
                                      "y", self.half_span)
        return root, tip
    
    @Part
    def sections(self):
        return Section(quantify=len(self.chords),
                       airfoil_name=self.airfoil,
                       chord=self.chords[child.index],
                       position=self.section_positions[child.index])
    
    @Part
    def surface(self):
        return LoftedShell(profiles=[section.curve for section in self.sections],
                           ruled=True)
    
    @Part
    def mirrored(self):
        return MirroredSurface(surface_in=self.surface.faces[0],
                               reference_point=self.positon.point,
                               vector1=self.position.Vx,
                               vector2=self.position.Vz,
                               suppress=not self.is_mirrored)
    
    @Attribute
    def root_chord(self):
        return self.chords(self)[0]

    @Attribute
    def mac(self):
        return (2*self.root_chord/3)*((1 + self.taper_ratio + self.taper_ratio**2)/(1 + self.taper_ratio))



if __name__ == '__main__':
    from parapy.gui import display
    obj = Wing(name='wing',
               span=20,
               aspect_ratio=6,
               taper_ratio=0.2,
               le_sweep=30,
               airfoil='2024',
               is_mirrored=False)
    display(obj)