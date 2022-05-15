import profile
from parapy.core import *
from parapy.geom import *
import numpy as np

"""
TO DO LIST:
    # 
"""

"""
RULE VIOLATION LIST:
    
"""
"""
variable = Input(derived)
@variable.getter
def variable(self):
    ...
    return ...

Range = Attribute(XXX)

@Part(in_tree=False)
"""

class FuelTank(GeomBase):
    Range = Input(6000) # [km] (A320 Range)
    V_cr = Input(230) # [m/s]
    C_T = Input(1.8639e-4) # [N/Ns] (MDO fuel consumption)
    LD = Input(16) # [-] (To be replaced from aero analysis)
    cruise_frac = Input(1.35) # Start / End of Cruise
    W_TO_max = Input(50e3) # [kg] MTOW
    rho_fuel = Input(0.81715e3) #[kg/m3] Fuel Density

    l_tank = Input(6) # [m]
    w_tank = Input(3) # [m]
    w_cabin = Input(12) # [m]

    wingspan = Input(25) # [m]
    sweep = Input(30) # [deg]

    loc = Input(5)

    def cruise_frac(self):
        W_cruise_pre = (self.Range * self.C_T) * (self.V_cr * self.LD)
        W_cruise = np.exp(W_cruise_pre)
        return W_cruise
    
    def W_fuel(self):
        W_fuel = (1 -0.938*self.cruise_frac)*self.W_TO_max
        return W_fuel

    # Coor. Sys Origin:
    @Attribute
    def origin(self):
        return XOY.translate("y",self.loc)
    
    @Part
    def origin_marker(self):
        return Sphere(radius=0.25,position=self.origin,color="orange")

    @Attribute
    def tank_span(self):
        return 0.5*self.wingspan

    @Attribute
    def tank_sweep(self):
        return 0.5*self.tank_span*np.cos(self.sweep)
    

    ## Fuel Tank Left:
    @Attribute
    def tank_root1(self):
        return Rectangle(width=self.w_tank,length=self.l_tank,
                        position=self.origin.translate("x",self.w_cabin/2).rotate("y",90,deg=True))
    
    @Attribute
    def tank_root2(self):
        return Rectangle(width=self.w_tank,length=self.l_tank,
                        position=self.origin.translate("x",-1*self.w_cabin/2).rotate("y",90,deg=True))
    
    @Attribute
    def tank_tip1(self):
        X_span = self.tank_span + self.w_cabin/2
        Y_span = self.tank_sweep
        return Rectangle(width=self.w_tank,length=self.l_tank,
                        position=self.origin.translate("x",X_span).translate("y",Y_span).rotate("y",90,deg=True))
    
    @Attribute
    def tank_tip2(self):
        X_span = -1*(self.tank_span + self.w_cabin/2)
        Y_span = self.tank_sweep
        return Rectangle(width=self.w_tank,length=self.l_tank,
                        position=self.origin.translate("x",X_span).translate("y",Y_span).rotate("y",90,deg=True))

    @Part
    def fuel_tank1(self):
        return RuledSolid(profile1=self.tank_root1, profile2=self.tank_tip1,color="orange")
    
    @Part
    def fuel_tank2(self):
        return RuledSolid(profile1=self.tank_root2, profile2=self.tank_tip2, color="orange")


if __name__ == '__main__':
    from parapy.gui import display

    obj = FuelTank()
    display(obj)