from parapy.core import Input, Part
from parapy.geom import GeomBase, ScaledCurve

from kbeutils.geom import Naca4AirfoilCurve


class Section(GeomBase):
    airfoil_name = Input()
    chord = Input()

    @Part(in_tree=False)
    def airfoil(self):
        return Naca4AirfoilCurve(designation=self.airfoil_name)

    @Part
    def curve(self):
        return ScaledCurve(curve_in=self.airfoil,
                           reference_point=self.position.point,
                           factor=self.chord)

if __name__ == "__main__":
    from parapy.gui import display

    obj = Section(airfoil_name="0012",
                  chord = 2)
    display(obj)