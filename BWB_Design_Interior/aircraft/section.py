from parapy.geom import *
from parapy.core import *
from parapy.core.validate import AdaptedValidator

import kbeutils.avl as avl
# from kbeutils.geom import Naca4AirfoilCurve
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve

_len_4_or_5 = AdaptedValidator(lambda a: 4 <= len(a) <= 5)


class Section(GeomBase):
    airfoil_name = Input(validator=_len_4_or_5)
    chord = Input()

    control_name = Input(None)
    control_hinge_loc = Input(None)
    duplicate_sign = Input(1)

    @Attribute
    def avl_controls(self):
        if self.control_name is not None and self.control_hinge_loc is not None:
            return [avl.Control(name=self.control_name,
                                gain=1.0,
                                x_hinge=self.control_hinge_loc,
                                duplicate_sign=self.duplicate_sign)]
        else:
            return []

    @Part
    def airfoil(self):
        return DynamicType(type=Naca5AirfoilCurve if len(self.airfoil_name) == 5
                                                     else Naca4AirfoilCurve,
                           designation=self.airfoil_name,
                           mesh_deflection=0.00001,
                           hidden=True)

    @Part
    def curve(self):
        return ScaledCurve(self.airfoil,
                           self.position.point,
                           self.chord,
                           mesh_deflection=0.00001
                           )

    # @Part  # the camber of the airfoil is ignored
    # def avl_section_no_curvature(self):
    #     return avl.Section(chord=self.chord)

    # @Part  # the camber of the airfoil is accounted, but works only for NACA4 airfoil (AVL
    # limitation) def avl_section(self): return avl.Section(chord=self.chord,
    # airfoil=avl.NacaAirfoil(designation=self.airfoil_name))

    # @Part  # the camber of the airfoil is accounted. Any curve is allowed
    # def avl_section(self):
    #     return avl.Section(chord=self.chord,
    #                        airfoil=avl.DataAirfoil(self.curve.sample_points))

    @Part  # the camber of the airfoil is accounted. Any curve is allowed. It includes
    # avl_section_by_points
    def avl_section(self):
        return avl.SectionFromCurve(curve_in=self.curve,
                                    controls=self.avl_controls)


if __name__ == '__main__':
    from parapy.gui import display

    obj = Section(airfoil_name="2412",
                  chord=2,
                  label="wing section",
                  )
    display(obj)
