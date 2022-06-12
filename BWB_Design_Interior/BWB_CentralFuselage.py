import os
from parapy.core import *
from parapy.geom import *
#from parapy.geom.ref_frame import Frame

from BWB_Cabin import Cabin
from BWB_Cargo import Cargo
from BWB_Fuel import FuelTank
from BWB_Wing import Wing
#from BWB_Class_I import AircraftMission
#from bwb_1 import BWB
import numpy as np
from kbeutils.geom.curve import Naca5AirfoilCurve, Naca4AirfoilCurve
import kbeutils.avl as avl
from parapy.exchange.step import STEPWriter

DIR = os.path.dirname(__file__)

"""
TO DO LIST:
    # Retrieve fused geometry attributes such as:
        - Coordinates of vertices
        - Total Length
"""

class CentralFuselage(GeomBase,avl.Interface):
    ### User Inputs ###
    #mission = Input()

    # Cabin Inputs:
    N_pax  = Input(700) # 200
    N_crew = Input(12)
    pax_dist = Input("5-10-85")
   
    airfoil_name = Input(str(2312))
    # Cargo Hold Inputs:
    m_cargo = Input(50e3)

    # Fuel Tank Inputs:
    W_fuel = Input(250e3)
    

    #total_length = 22.5

    # Wing Inputs
    taper = Input(0.55)
    surface_scale=Input(65)
    wingspan = (60)
    Mach = Input(0.6)

    # AVL:
    case_settings = Input()

    ### Set Constants ###
    """ REMOVE !!!"""
    N_col      = Input(3) # 3
    h_cabin    = Input(2.1) # 2.5
    w_seat     = Input(0.508)  # 0.508m ;20 inches (HAW Hamburg)
    w_aisle    = Input(0.49) # 0.49m; 19 inches (HAW Hamburg)
    l_row      = Input(0.89) # 0.89m ;35 inches (HAW Hamburg)
    w_cockpit = Input(3) #3
    l_cockpit  = Input(5) #6
    rho_cargo = Input(160)
    w_cargo = Input(6)
    l_cargo = Input(15)
    h_cargo = Input(1.5)

    # Coor. Sys Origin:
    # @Part
    # def Frame(self):
    #     return Frame(pos=XOY, hidden=False)
    
    # @Attribute
    # def mission(self):
    #     return AircraftMission()

    frame = Input(XOY)

    @Attribute
    def origin_points(self):
        XYZ_cabin = self.frame.rotate("z",90,deg=True)
        y_tilde = ((-self.h_cargo/2)*(self.h_cargo*self.w_cargo)+\
                        (self.h_cabin/2)*(self.h_cabin*self.cabin.w_cabin))\
                 /((self.h_cargo*self.w_cargo)+(self.h_cabin*self.cabin.w_cabin))
        XYZ_front = self.frame.translate("z",y_tilde,"y",-self.l_cockpit*1.1).rotate("z",90,deg=True)
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
        return Cabin(N_pax=self.N_pax, N_col_ec=self.N_col,
                    h_cabin=self.h_cabin, w_seat_ec=self.w_seat, w_aisle_ec=self.w_aisle,
                    w_cockpit=self.w_cockpit, l_cockpit=self.l_cockpit,
                    loc=XOY.rotate("z",-90,deg=True))
    
    @Part
    def cargo(self):
        return Cargo(w_cabin=self.cabin.w_cabin,l_cabin=self.cabin.l_cabin,
                     h_cargo=self.h_cargo,m_cargo=self.m_cargo,
                     rho_cargo=self.rho_cargo,
                     loc=XOY)
    
    @Attribute
    def fuel_weight(self):
        return self.W_fuel
    
    @Attribute
    def cog_interior(self):
        # CG_cabin = self.cabin.fused_cabin.cog
        # CG_cargo = self.cargo.cargo_hold.cog
        # cg_int = CG_cargo - CG_cabin
        return Position(self.cabin.fused_cabin.cog)

    @Part 
    def fuel_tanks(self):
        return FuelTank(w_cabin=self.cabin.w_cabin,
                        l_cabin=self.cabin.l_cabin,
                        h_cabin=self.cabin.h_cabin,
                        W_fuel = self.W_fuel,
                        loc=self.cog_interior.rotate("z",-90,deg=True),
                        wingspan=self.wingspan)
                        #h_tank=self.cabin.h_cabin)
    
    @Attribute
    def Interior_Corners(self):
        #vrtx_coor_s = vrtx_coor[np.argsort(vrtx_coor[:,1])]
        ## Use argsort to sort along x axis instead
        cor_cab = self.cabin.corners[1][np.argsort(self.cabin.corners[1][:,0])]
        cor_carg = self.cargo.corners[1][np.argsort(self.cargo.corners[1][:,0])]
        fuel_cor = self.fuel_tanks.corners[1][np.argsort(self.fuel_tanks.corners[1][:,0])]
        #print(cor_cab)
        """To extract corners on left side, set if statement """
        # return [self.cabin.corners,
        #         self.cargo.corners,
        #         self.fuel_tanks.corners]
        return [cor_cab,cor_carg,fuel_cor]
    
    @Attribute
    def symm_corners(self):
        cor_l = []
        cor_r = []
        for part in self.Interior_Corners:
            for coor in part:
                point = Point(coor[0],coor[1],coor[2])
                if coor[1] > 0:    
                    cor_l.append(point)
                elif coor[1] < 0:
                    cor_r.append(point)
        #print(cor_l)
        return [cor_l,cor_r]

    ### Wing Surface
    """
    @Part(parse=False)
    def WingSurface(self):
        return Wing(origin=self.origin_points,
                    cabin=self.cabin,
                    fuel_tanks=self.fuel_tanks,
                    wingspan=self.wingspan)

    ### AVL Setup:

    @Attribute
    def avl_surfaces(self):
        return self.find_children(lambda o: isinstance(o, avl.Surface))

    @Part
    def avl_configuration(self):
        return avl.Configuration(name='aircraft',
                                 reference_area=self.WingSurface.wing_area,
                                 reference_span=self.WingSurface.wing_span+self.cabin.w_cabin,
                                 reference_chord=self.WingSurface.wing_mac,
                                 reference_point=self.WingSurface.wing_cog,
                                 surfaces=self.avl_surfaces,
                                 mach=self.Mach)
    
    @Attribute
    def configuration(self):
        return self.avl_configuration
    
    @Part
    def cases(self):
        return avl.Case(quantify=len(self.case_settings),
                        name=self.case_settings[child.index][0],
                        settings=self.case_settings[child.index][1])

    @Attribute
    def l_over_d(self):
        return {case_name: result['Totals']['CLtot'] / result['Totals']['CDtot']
                for case_name, result in self.results.items()}
    """
    ### Export Files to STP:
    # Uncomment .write() method
    # @Attribute
    # def CAD_Export(self):
    #     print(DIR)
    #     return STEPWriter(trees=[self.cabin, self.cabin,
    #                             self.fuel_tanks,self.WingSurface],
    #                       default_directory=DIR,
    #                       filename="BWB_CAD.stp")
    
if __name__ == '__main__':
    from parapy.gui import display

    cases = [('fixed_aoa', {'alpha': 3}), 
            ('fixed_cl', {'alpha': avl.Parameter(
                name='alpha', value=0.3, setting='CL')})]

    obj = CentralFuselage(case_settings=cases)

    
    
    #obj.CAD_Export.write()
    display(obj)
    #display(obj.cases)
    #display(obj.ellipses)