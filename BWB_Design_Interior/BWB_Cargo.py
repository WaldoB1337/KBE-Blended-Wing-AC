from parapy.core import *
from parapy.geom import *

"""
TO DO LIST:
    # 
"""

class Cargo(GeomBase):
    w_cargo = Input(6) #6
    l_cargo = Input(15) #15
    h_cargo = Input(1.5) #1.5

    m_cargo = Input(10e3)# 10e3
    rho_cargo = Input(160) #160 

    # Coor. Sys Origin:
    @Part
    def origin(self):
        return Sphere(radius=0.1,position=XOY)

    @Attribute
    def cargo_vol(self):
        return self.m_cargo / self.rho_cargo
    
    @Attribute
    def cargo_length(self):
        return self.l_cargo
    
    @Attribute
    def cargo_width(self):
        return self.w_cargo

    #@Attribute
    def cargoSizeToPayload(self):
        V_cargo = self.l_cargo * self.w_cargo * self.h_cargo
        m_payload = V_cargo /self.rho_cargo
        return m_payload
    
    #@Attribute
    def payloadToCargoSize(self): ## Looking for relations...
        pass

    @Part
    def cargo_hold(self):
        return Box(length=self.l_cargo,width=self.w_cargo,height=self.h_cargo,
                   centered=False, position=XOY.rotate("y",180,deg=True).translate("x",-self.w_cargo/2))
    
if __name__ == '__main__':
    from parapy.gui import display

    obj = Cargo()
    display(obj)