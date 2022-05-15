from parapy.core import *
from parapy.geom import *
import numpy as np

"""
TO DO LIST:
    # Retrieve Fuse part vertices for fuselage ellipse radius
    # Update sizing laws ?
    # Search for Cockpit Size Regulations
    # Setup checks for integers
    # Establish rules for consistent geometry
"""

"""
RULE VIOLATION LIST:
    # Cockpit Width <= 0.5*Cabin Width
        - if not, looks wonky; 
"""

class Cabin(GeomBase):
    N_pax      = Input(250) # 200
    N_col      = Input(2) # 3
    h_cabin    = Input(2.1) # 2.5
    w_seat     = Input(0.508)  # 0.m ;20 inches (HAW Hamburg)
    w_aisle    = Input(0.49) # 0.49m; 19 inches (HAW Hamburg)
    l_row      = Input(0.89) # 0.89m ;35 inches (HAW Hamburg)
    w_cockpit = Input(2) #3
    l_cockpit  = Input(5) #6

    # Coor. Sys Origin:
    @Part
    def origin(self):
        return Sphere(radius=0.1,position=XOY)

    @Attribute
    def paxToCabinSize(self):
        # Under construction...
        pass

    @Attribute
    def cabinSizeToPax(self):
        # Under construction...
        pass

    @Attribute
    def l_cabin(self):
        N_pax_row = 6 * self.N_col
        N_row = int(np.ceil(self.N_pax / N_pax_row))
        return N_row * self.l_row
    
    @Attribute
    def w_cabin(self):
        self.w_col = 6*self.w_seat + self.w_aisle
        return self.w_col * self.N_col
    
    
    
    # @Attribute
    # def cabin_vol(self):
    #     return self.cabin_space

    @Part
    def cabin_main(self):
        return Box(length=self.l_cabin,width=self.w_cabin,height=self.h_cabin,centered=False,
                   position=XOY.translate("x",-self.w_cabin/2))

    @Part(parse=False)
    def cabin_1c(self): ## First Class front trapezoid
        self.x_min = (self.w_cabin/2) - (self.w_cockpit/2)
        self.x_max = (self.w_cabin/2) + (self.w_cockpit/2)
        return Wedge(dx=self.w_cabin, dy=self.l_cockpit, dz=self.h_cabin, xmin=self.x_min, xmax=self.x_max, zmin=0, zmax=self.h_cabin,
                    position=XOY.translate("x",self.w_cabin).rotate("z",180,deg=True).translate("x",self.w_cabin/2))

    @Part(parse=False)
    def aft_bay(self): ## First Class front trapezoid
        self.x_min = (self.w_cabin/2) - (self.w_cockpit/2)
        self.x_max = (self.w_cabin/2) + (self.w_cockpit/2)
        return Wedge(dx=self.w_cabin, dy=self.l_cockpit/2, dz=self.h_cabin, xmin=self.x_min, xmax=self.x_max, zmin=0, zmax=self.h_cabin,
                    position=XOY.translate("y",self.l_cabin).translate("x",-self.w_cabin/2))
    
    @Part
    def fused_cabin(self):
        return FusedSolid(color="blue",shape_in=self.cabin_main,tool=self.cabin_1c + self.aft_bay,
                            position=XOY.translate("x",-self.w_cabin/2))
    
    @Attribute
    def l_cabin_tot(self):
        return 1.5*self.l_cockpit + self.l_cabin

if __name__ == '__main__':
    from parapy.gui import display

    obj = Cabin()
    display(obj)