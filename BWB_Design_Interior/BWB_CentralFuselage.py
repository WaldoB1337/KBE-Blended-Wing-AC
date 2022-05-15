from parapy.core import *
from parapy.geom import *
#from parapy.geom.ref_frame import Frame

from BWB_Cabin import Cabin
from BWB_Cargo import Cargo
from BWB_Fuel import FuelTank
import numpy as np
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve


"""
TO DO LIST:
    # Retrieve fused geometry attributes such as:
        - Coordinates of vertices
        - Total Length
"""

class CentralFuselage(GeomBase):
    #Cabin Inputs
    N_pax      = Input(250) # 200
    N_col      = Input(3) # 3
    h_cabin    = Input(2.1) # 2.5
    w_seat     = Input(0.508)  # 0.508m ;20 inches (HAW Hamburg)
    w_aisle    = Input(0.49) # 0.49m; 19 inches (HAW Hamburg)
    l_row      = Input(0.89) # 0.89m ;35 inches (HAW Hamburg)
    w_cockpit = Input(3) #3
    l_cockpit  = Input(5) #6
    airfoil_name = Input(str(2312))
    # Cargo Hold Inputs:
    w_cargo = Input(6)
    l_cargo = Input(15)
    h_cargo = Input(1.5)

    m_cargo = Input(10e3)
    rho_cargo = Input(160)

    total_length = 22.5

    # Coor. Sys Origin:
    # @Part
    # def Frame(self):
    #     return Frame(pos=XOY, hidden=False)

    @Attribute
    def origin_points(self):
        XYZ_cabin = XOY
        y_tilde = ((-self.h_cargo/2)*(self.h_cargo*self.w_cargo)+\
                        (self.h_cabin/2)*(self.h_cabin*self.cabin.w_cabin))\
                 /((self.h_cargo*self.w_cargo)+(self.h_cabin*self.cabin.w_cabin))
        XYZ_front = XOY.translate("z",y_tilde,"y",-self.l_cockpit*1.1)
        return [XYZ_front,XYZ_cabin]
    
    # @Attribute
    # def cabin_dimensions(self):
    #     cab_dim = {"l_cabin":self.cabin.l_cabin,
    #                "w_cabin":self.cabin.w_cabin,
    #                "l_cabin_tot":self.cabin.l_cabin_tot,}
    #     return cab_dim
    
    @Part 
    def origin_front(self):
        """HELP!: How can I output all the markers from origin_points in one go ?"""
        return Sphere(radius=0.1,quantify=len(self.origin_points),
                      position=self.origin_points[0])
    @Part
    def cabin(self):
        return Cabin(N_pax=self.N_pax, N_col=self.N_col,
                    h_cabin=self.h_cabin, w_seat=self.w_seat, w_aisle=self.w_aisle,
                    w_cockpit=self.w_cockpit, l_cockpit=self.l_cockpit)
    
    @Part
    def cargo(self):
        return Cargo(w_cargo=self.w_cargo,l_cargo=self.cabin.l_cabin,
                     h_cargo=self.h_cargo,m_cargo=self.m_cargo,
                     rho_cargo=self.rho_cargo)
    # @Part
    # def fused_interior(self):
    #     return FusedSolid(color="blue",shape_in=self.cabin,tool=self.cargo,
    #                         position=self.origin_front)
    @Part 
    def fuel_tanks(self):
        return FuelTank(w_cabin=self.cabin.w_cabin,
                        loc=self.cabin.l_cabin_tot/2)
    
    ## Inner Shell:
    """
    @Attribute
    def shell_ellipses(self):
        #https://stackoverflow.com/questions/433371/ellipse-bounding-a-rectangle
        # Assume larger rectangle with dimensions dictated by Cargo + Cabin extremes
        Rh = self.h_cabin + self.h_cargo
        Rw = self.cabin.w_cabin \
            if self.cabin.w_cabin > self.w_cargo else self.w_cargo
        # 5% (?) margin between internal structure and outer shell
        Fuse_margin_a = 1.0  
        Fuse_margin_b = 1.0 
        # Ellipse Parameters:
        a = Fuse_margin_a * (Rw/np.sqrt(2)) # Major Axis
        b = Fuse_margin_b * (Rh/np.sqrt(2)) # Minor Axis
        #e1 = Ellipse(major_radius=a, minor_radius = b, 
        #            position=self.origin_points[0].rotate("x",90,deg=True))
        #e2 = Ellipse(major_radius=a, minor_radius = b, 
        #            position=self.origin_points[0].translate("y",22.5).rotate("x",90,deg=True))
        reso = 10
        shell_rings = []
        l_tot = self.total_length
        axis = np.linspace(0,l_tot,reso)
        for step in range(reso):
            ring = Ellipse(major_radius=a, minor_radius = b, 
                    position=self.origin_points[0].translate("y",axis[step]).rotate("x",90,deg=True))
            shell_rings.append(ring)
        return shell_rings#e1, e2

    @Part(parse=False)
    def InnerShell(self): # Inner Fuselage Skin
        #Returns a single surface through multiple curves
        return self.shell_ellipses
        #return LoftedShell(profiles=self.shell_ellipses)
    
    @Part(parse=False)
    def OuterShell(self):
        Airfoil = DynamicType(type=Naca5AirfoilCurve \
                           if len(self.airfoil_name) == 5 else Naca4AirfoilCurve,
                           designation=self.airfoil_name,
                           mesh_deflection=0.00001,
                           hidden=True)

        return ScaledCurve(Airfoil,
                           XOY,
                           self.total_length,
                           mesh_deflection=0.00001
                           )
        # return FittedCurve(curve_in=Airfoil, 
        #                  reference_point=XOY, 
        #                  factor=1)
"""
if __name__ == '__main__':
    from parapy.gui import display

    obj = CentralFuselage()
    display(obj)
    #display(obj.ellipses)