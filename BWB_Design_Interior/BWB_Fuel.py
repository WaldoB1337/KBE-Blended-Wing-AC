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
    - Fuel Tank Chord <= 80% Main Cabin Length
    
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
def error_message(header,msg):
        from tkinter import Tk, mainloop, X, messagebox

        window = Tk()
        window.withdraw()

        messagebox.showwarning(header,msg)

class FuelTank(GeomBase):


    ### User Inputs: ###
    wingspan = Input(60) # [m]
    sweep = Input(30) # [deg]
    taper = Input(0.4)
    W_fuel = Input(250e3)

    ## Constant:###
    h_tank = Input(0.5) # [m]
    w_cabin = Input(12) # [m]
    l_cabin = Input(20)
    h_cabin = Input(2.1)
    rho_fuel = Input(0.81715e3) # [kg/m3]
    h_cargo = Input(1.5)
    tank_c_root = Input(15)
    
    # Location passed as Cabin COG in CentralFuselage
    loc = Input(XOY)

    # Coor. Sys Origin:
    @Attribute
    def origin(self):
        return self.loc
    
    @Part
    def origin_marker(self):
        return Sphere(radius=0.25,position=self.origin,color="orange")
    
    @Attribute
    def tank_span(self):
        return 0.5*self.wingspan*0.8 # Constraint tank to 80\% of Half wing span
    
    @Attribute
    def fuel_vol(self):
        #V_fuel = (0.5*self.wingspan*self.tank_c_r)*(1+self.taper)*self.h_tank
        V_fuel = self.W_fuel / self.rho_fuel
        return V_fuel
    
    @Attribute
    def tank_height(self):
        h_tank = self.h_tank
        if h_tank > self.h_cabin + 0.5*self.h_cargo:
            header = "Warning: Fuel Tank Exceeds Limit!"
            msg = "The depth of the fuel tank exceeds preset boundary.\
                Dimensions will be reset to match the cabin height.\
                    Consider entering a larger wingspan to compensate."
            error_message(header,msg)
            h_tank = self.h_cabin + 0.5*self.h_cargo#self.h_cabin 

        return h_tank
    
    @Attribute
    def tank_c_r(self):
        """In case Tank Root exceeds cabin, take minimum"""
        root_chord = self.tank_c_root

        if root_chord > 0.75*self.l_cabin:
            header = "Warning: Fuel Tank Exceeds Limit!"
            msg = "The root chord of the fuel tank exceeds that of the main cabin.\
                It will be reset to be 75% \of the main cabin length.\
                    Consider entering a larger wingspan to compensate."
            error_message(header,msg)
            root_chord = 0.75*self.l_cabin
        return root_chord
    
    

    @Attribute
    def tank_sweep(self):
        return self.tank_span*np.tan(np.radians(self.sweep)) # Y-Axis shift due to wing sweep
    

    ## Fuel Tank Left:
    @Attribute
    def tank_root1(self):
        return Rectangle(width=self.tank_height,length=self.tank_c_r,
                        position=self.origin.translate("x",self.w_cabin/2).rotate("y",90,deg=True))
    
    @Attribute
    def tank_root2(self):
        return Rectangle(width=self.tank_height,length=self.tank_c_r,
                        position=self.origin.translate("x",-1*self.w_cabin/2).rotate("y",90,deg=True))
    
    @Attribute
    def tank_tip1(self):
        X_span = self.tank_span + self.w_cabin/2
        Y_span = self.tank_sweep
        return Rectangle(width=self.tank_height,length=self.tank_c_r*self.taper,
                        position=self.origin.translate("x",X_span).translate("y",Y_span).rotate("y",90,deg=True))
    
    @Attribute
    def tank_tip2(self):
        X_span = -1*(self.tank_span + self.w_cabin/2)
        Y_span = self.tank_sweep
        return Rectangle(width=self.tank_height,length=self.tank_c_r*self.taper,
                        position=self.origin.translate("x",X_span).translate("y",Y_span).rotate("y",90,deg=True))

    @Part
    def fuel_tank1(self):
        return RuledSolid(profile1=self.tank_root1, profile2=self.tank_tip1,color="orange")
    
    @Part
    def fuel_tank2(self):
        return RuledSolid(profile1=self.tank_root2, profile2=self.tank_tip2, color="orange")

    @Attribute
    def corners(self):
        vertices = []
        corner_mark = []
        vrtx_coor = []
        components = [self.fuel_tank1.vertices,
                     self.fuel_tank2.vertices]
        for component in components:
            for i in range(len(component)):
                vertex_cargo = component[i].point
                vertices.append(vertex_cargo)
                vrtx_coor.append([vertex_cargo.x,vertex_cargo.y,vertex_cargo.z])

        vrtx_coor = np.array(vrtx_coor)
        #print(np.array(vrtx_coor))
        vrtx_coor_s = vrtx_coor[np.argsort(vrtx_coor[:,1])]

        return [vertices, vrtx_coor_s]
    
    @Attribute
    def midline(self):
        midline_1 = self.fuel_tank1.profile2.center - self.fuel_tank1.profile1.center
        midline_2 = self.fuel_tank2.profile2.center - self.fuel_tank2.profile1.center
        #print(midline_1)
        #print(midline_2)
        
        return [midline_1,midline_2]

    @Part(parse=False)
    def cargo_vertices(self):
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
        data = self.corners[0]
        ctrl_loc = []
        ## Add manually upper points
        ctrl_points_fr = [1,2,9,10] ## Pointing forward
        ctrl_points_rr = [4,5,12,13] ## Pointing Rearward


        for i in range(0,len(data)):
            if i in ctrl_points_fr:
            #if i==1 or i==2 or i==9 or i==10:
                point = Point(data[i][0],data[i][1],data[i][2])\
                    .translate("z",0.5*self.tank_height,
                                "x",-1.2*0.5*np.sqrt(3)*self.tank_height)
                ctrl_loc.append(point)
            elif i in ctrl_points_rr:
            #if i==1 or i==2 or i==9 or i==10:
                point = Point(data[i][0],data[i][1],data[i][2])\
                    .translate("z",0.5*self.tank_height,
                                "x",1.2*0.5*np.sqrt(3)*self.tank_height)
                ctrl_loc.append(point)
            
                
        # for i in range(len(ctrl_loc)):
        #     ctrl_points.append(Sphere(radius=0.2,position=ctrl_loc[i], color="black"))

        return ctrl_loc
    
    @Part(parse=False)
    def tank_border(self):
        ## LE and TE
        reso = 10
        data = self.control_points
        #print(len(data))
        vectors,points,ref_points,spheres = [],[],[],[]
        for i in range(int(0.5*len(data))):
            line_vec = (data[2*i+1] - data[2*i])/ reso
            vectors.append(line_vec)
            ref_points.append(data[2*i])

        for i in range(len(vectors)):
            for k in range(reso):
                point = ref_points[i] + (k+1)*vectors[i]
                #print(point)
                spheres.append(Sphere(radius=0.2,
                                position=point, color="black"))
        
        ## Tank Tips
        #ref_points = [2,5,10,13]
        reso = 5
        corners = self.corners[0]
        line_vec_l = (corners[5] - corners[2]) / reso
        line_vec_r = (corners[13] - corners[10]) / reso
        ref_points = [corners[2],corners[10]]
        vectors = [line_vec_l,line_vec_r]
        for ref in range(len(ref_points)):
            for i in range(reso+1):
                if ref==0:
                    buffer = -1.2*0.5*np.sqrt(3)*self.tank_height
                else:
                    buffer = 1.2*0.5*np.sqrt(3)*self.tank_height
                point = ref_points[ref] + i*vectors[ref] + Vector(0,buffer,0.5*self.tank_height)
                # point_l = corners[2] + i*line_vec_l
                # point_r = corners[10] + i*line_vec_r
                marker = Sphere(radius=0.2,position=point, color="black")
                spheres.append(marker)
        return spheres
    
    @Part(parse=False)
    def control_spheres(self):
        ctrl_spheres = []
        for i in range(len(self.control_points)):
            ctrl_spheres.append(Sphere(radius=0.2,
            position=self.control_points[i], color="black"))

        return ctrl_spheres


if __name__ == '__main__':
    from parapy.gui import display

    obj = FuelTank()
    display(obj)