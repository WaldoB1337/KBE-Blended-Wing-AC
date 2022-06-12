from parapy.core import *
from parapy.geom import *
import numpy as np

"""
TO DO LIST:
    # 
"""
"""
RULE VIOLATION LIST:
    - Cargo Length == main Cabin Length
    - Cargo Width <= 85% Cabin Width
    
"""
def error_message(header,msg):
        from tkinter import Tk, mainloop, X, messagebox

        window = Tk()
        window.withdraw()

        messagebox.showwarning(header,msg)

class Cargo(GeomBase):

    ### User Inputs : ###
    m_cargo = Input(10e3)# 10e3
    rho_cargo = Input(160) #160

    w_cargo = Input()
    h_cargo = Input(2.1) # Standrd LD3 Container Height = 162.6 cm
    l_cargo = Input()

    ### Constants: ### 
    w_cabin = Input(6) #6
    l_cabin = Input(15) 
    #h_cargo = Input(1.626) 

    loc = Input(XOY.rotate("z",90,deg=True))

    # Coor. Sys Origin:
    @Part
    def origin(self):
        return Sphere(radius=0.1,position=self.loc)

    @Attribute
    def cargo_vol(self):
        return self.m_cargo / self.rho_cargo
    
    @Attribute
    def cargo_length(self):
        # Set to Central Cabin Length in Central Fuselage
        l_cargo = self.l_cabin
        if l_cargo >= 0.558*self.l_cabin:
            header = "Warning: Cargo Bay Exceeds Airfoil!"
            msg = "The length of the cargo bay exceeds the wing surface.\
                Length will be reset to ??% \that of the Cabin.\
                    Consider entering... ."
            error_message(header,msg)
            l_cargo = 0.558*self.l_cabin        
        return l_cargo
    
    @Attribute
    def cargo_height(self):
        #h_cargo = self.cargo_vol / (self.cargo_width + self.cargo_length)
        return self.h_cargo
    
    @Attribute
    def cargo_width(self):
        w_cargo = self.cargo_vol / (self.cargo_height + self.cargo_length)
        
        if w_cargo >= 0.85*self.w_cabin:
            header = "Warning: Cargo Bay Exceeds Limits!"
            msg = "The width of the cargo bay exceeds preset boundary.\
                Width will be reset to 85% \that of the Cabin.\
                    Consider entering... ."
            error_message(header,msg)
            w_cargo = 0.85*self.w_cabin
        return w_cargo
    
    
    @Attribute
    def cargoSizeToPayload(self):
        V_cargo = self.cargo_length * self.cargo_width * self.h_cargo
        m_payload = V_cargo /self.rho_cargo
        return m_payload
    
    #@Attribute
    def payloadToCargoSize(self): ## Looking for relations...
        pass

    @Part
    def cargo_hold(self):
        return Box(length=self.cargo_length,width=self.cargo_width,height=self.cargo_height,
                   centered=False, 
                   position=self.loc.translate("y",self.cargo_width/2,"z",-self.cargo_height).rotate("z",-90,deg=True))
    
    @Attribute
    def corners(self):
        vertices = []
        corner_mark = []
        vrtx_coor = []
        component = self.cargo_hold.vertices
        #for component in components: ## In case of multiple objects
        for i in range(len(component)):
            vertex_cabin = component[i].point
            vertices.append(vertex_cabin)
            vrtx_coor.append([vertex_cabin.x,vertex_cabin.y,vertex_cabin.z])

        vrtx_coor = np.array(vrtx_coor)
        #print(np.array(vrtx_coor))
        vrtx_coor_s = vrtx_coor[np.argsort(vrtx_coor[:,0])]

        return [vertices, vrtx_coor_s]
    
    @Part(parse=False)
    def Cabin_vertices(self):
        #return self.corners[0]
        vertices = self.corners[0]
        corner_mark = []
        for corner in range(len(vertices)):
            marker = Sphere(radius=0.2,position=vertices[corner], color="red")
            corner_mark.append(marker)
        return corner_mark

if __name__ == '__main__':
    from parapy.gui import display

    obj = Cargo()
    display(obj)