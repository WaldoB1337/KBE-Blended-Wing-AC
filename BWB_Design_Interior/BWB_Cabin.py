from turtle import position
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

    ### User Inputs: ###
    N_pax      = Input(700) # 200
    pax_dist   = Input("3-12-85") # Pax distribution [%]: Emirates B777
    # Retrieved from: https://www.pngkey.com/detail/u2e6y3q8t4r5e6q8_seatmap-777-300er-garuda-boeing-777-seat-map/


    ### Fine Tuning: ###
    # Standard Cabin values according to Raymer. Retrieved from:
    # https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_6_Fuselage.pdf
    h_cabin     = Input(2.1) # 2.5
    # First Class:
    N_col_1st   = Input(2) # 3
    w_seat_1st  = Input(0.7112) # 28 inches (HAW Hamburg)
    w_aisle_1st = Input(0.7112) # 28 inches (HAW Hamburg)
    l_row_1st   = Input(1.016)   # 40 inches (HAW Hamburg)
    
    #Business: (Dimensions chosen as midpoint from Economy and First)
    N_col_bss   = Input(2) # 3
    w_seat_bss  = Input(0.508) # 20 inches 
    w_aisle_bss = Input(0.6096)  # 24 inches 
    l_row_bss   = Input(0.9652)  # 38 inches 
    
    #Economy
    N_col_ec   = Input(2) # 3
    w_seat_ec  = Input(0.5588) # 22 inches (HAW Hamburg)
    w_aisle_ec = Input(0.5080)  # 20 inches (HAW Hamburg)
    l_row_ec   = Input(0.9144)  # 36 inches (HAW Hamburg)
    
    # Cockpit & Aft Bulkhead:
    w_cockpit  = Input(2) #3
    l_cockpit  = Input(5) #6
    
    
    # Coor. Sys Origin:
    @Part
    def origin(self):
        return Sphere(radius=0.1,position=XOY)

    # @Attribute
    # def paxToCabinSize(self):
    #     # Under construction...
    #     pass

    # @Attribute
    # def cabinSizeToPax(self):
    #     # Under construction...
    #     pass

    @Attribute
    def N_pax_dist(self):
        dist = self.pax_dist.split("-")
        dist = [(self.N_pax/100)*float(x) for x in dist]
        return dist
    
    @Attribute
    def l_cabin(self):
        ## First:
        N_pax_row_1st = 6 * self.N_col_1st
        N_row_1st = int(np.ceil(( self.N_pax_dist[0])/ N_pax_row_1st))
        l_cabin_1st = N_row_1st * self.l_row_1st

        ## Business:
        N_pax_row_bss = 6 * self.N_col_bss
        N_row_bss = int(np.ceil((self.N_pax_dist[1])/ N_pax_row_bss))
        l_cabin_bss = N_row_bss * self.l_row_bss

        ## Economy Seats:
        N_pax_row_ec = 6 * self.N_col_ec
        N_row_ec = int(np.ceil((self.N_pax_dist[2])/ N_pax_row_ec))
        l_cabin_ec = N_row_ec * self.l_row_ec

        return l_cabin_1st + l_cabin_bss +  l_cabin_ec 
    
    @Attribute
    def w_cabin(self):
        w_seat_max = max(self.w_seat_1st,self.w_seat_bss,self.w_aisle_ec)
        self.w_col = 6*w_seat_max + self.w_aisle_ec
        return self.w_col * self.N_col_ec

    @Part
    def cabin_main(self):
        return Box(length=self.l_cabin,width=self.w_cabin,height=self.h_cabin,centered=False,
                   position=XOY.translate("x",-self.w_cabin/2))

    @Part(parse=False)
    def cockpit(self): ## First Class front trapezoid
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
        return FusedSolid(color="blue",shape_in=self.cabin_main,tool=self.cockpit + self.aft_bay,
                            position=XOY.translate("x",-self.w_cabin/2))
    
    @Attribute
    def l_cabin_tot(self):
        return 1.5*self.l_cockpit + self.l_cabin
    
    @Part
    def cabin_cg(self):
        return Sphere(radius=0.5,position=self.fused_cabin.cog,color="black")
    
    @Attribute
    def corners(self):
        vertices = []
        corner_mark = []
        vrtx_coor = []
        component = self.fused_cabin.vertices
        #for component in components:
        for i in range(len(component)):
            vertex_cabin = component[i].point
            vertices.append(vertex_cabin)
            vrtx_coor.append([vertex_cabin.x,vertex_cabin.y,vertex_cabin.z])

        vrtx_coor = np.array(vrtx_coor)
        #print(vertices)
        #print(np.array(vrtx_coor))
        vrtx_coor_s = vrtx_coor[np.argsort(vrtx_coor[:,1])]

        return [vertices, vrtx_coor_s]
    
    @Part(parse=False)
    def cabin_vertices(self):
        #return self.corners[0]
        vertices = self.corners[0]
        corner_mark = []
        for corner in range(len(vertices)):
            marker = Sphere(radius=0.2,position=vertices[corner], color="red")
            corner_mark.append(marker)
        return corner_mark
        
    
    @Attribute
    def control_points(self):
        #print(self.corners[1])
        data = self.corners[1]
        ctrl_loc = []
        ctrl_points = []
        for i in range(0,int(0.5*len(data))):
            if data[i][2] == 0:
                point = Point(data[i][0],data[i][1],data[i][2])\
                    .translate("z",0.2*self.h_cabin,
                                "y",-0.5*np.sqrt(3)*self.h_cabin*1.5)
                #print(point)
                ctrl_loc.append(point)

        # for i in range(len(ctrl_loc)):
        #     ctrl_points.append(Sphere(radius=0.2,position=ctrl_loc[i], color="black"))

        return ctrl_loc
    
    @Part(parse=False)
    def control_spheres(self):
        ctrl_spheres = []
        for i in range(len(self.control_points)):
            ctrl_spheres.append(Sphere(radius=0.2,
            position=self.control_points[i], color="black"))

        return ctrl_spheres

if __name__ == '__main__':
    from parapy.gui import display

    obj = Cabin()
    #print(obj.fused_cabin.cog)
    display(obj)