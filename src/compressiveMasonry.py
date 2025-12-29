from structuralMasonry import Masonry
import math

class Wall(Masonry):
    """ compressive load capacity of a wall """
    def __init__(self, fb:float, fm:float, k:float, effective_height:float,
                 loaded_leaf_thickness:float, other_leaf_thickness:float, top_eccentricity:float = 5):
        super().__init__(fb, fm, k)
        self.effective_height = effective_height # [mm]
        self.loaded_leaf_thickness = loaded_leaf_thickness # [mm]
        self.other_leaf_thickness = other_leaf_thickness # [mm]
        self.top_eccentricity = top_eccentricity # [mm]

        # calculate the effective thickness of the wall
        self.effective_thickness = (self.loaded_leaf_thickness**3 + self.other_leaf_thickness**3)**(1/3)
        # calculate the eccentricity at the top of the wall
        self.einit = self.effective_height/450
        self.et = max(self.einit + self.top_eccentricity, 0.05 * self.loaded_leaf_thickness)

        # calculate the capacity reduction factor at the top of the wall
        self.top_slenderness_reduction = 1 - 2 * (self.et / self.loaded_leaf_thickness)

        # calculate the eccentricity at mid-height of the wall
        self.em = self.einit + self.top_eccentricity / 2 # assuming a linear reduction
        self.ek = (0.002 * self.creep_coefficient * (self.effective_height/self.effective_thickness) *
                   (self.loaded_leaf_thickness * self.em)**0.50)
        self.emk = max(self.em + self.ek, 0.05 * self.loaded_leaf_thickness)

        # calculate the slenderness corrective factor at mid-height of the wall
        A1 = 1 - 2 * (self.emk / self.loaded_leaf_thickness)
        lam = (self.effective_height / self.effective_thickness) * (self.fk / self.Em) ** 0.50
        u = (lam - 0.063) / (0.73 - 1.17 * (self.emk / self.loaded_leaf_thickness))
        self.mid_slenderness_reduction = A1 * math.e ** (-u ** 2 / 2)

        # calculate the design vertical capacity of the wall [kN/m]

        self.design_vertical_capacity = (min(self.top_slenderness_reduction, self.mid_slenderness_reduction) *
                                         self.fk * self.loaded_leaf_thickness / self.gm)

    # Return the design vertical load capacity of the wall
    def print_design_vertical_capacity(self)->None:
        print(f"hef = {self.effective_height:.2f}mm, twall = {self.loaded_leaf_thickness:.2f}mm")
        print(f"fk = {self.fk:.2f}MPa : Phi_top = {self.top_slenderness_reduction:.2f} : "
              f"Phi_mid = {self.mid_slenderness_reduction:.2f}")
        print(f"Design vertical load capacity of the wall = "
              f"{self.design_vertical_capacity:.2f} kN/m")

class Padstone(Masonry):
    def __init__(self, fb:float, fm:float, k:float, near_edge_distance:float, far_edge_distance:float,
                    wall_height:float, padstone_length:float, padstone_width:float):
        super().__init__(fb, fm, k)
        self.near_edge_distance = near_edge_distance
        self.far_edge_distance = far_edge_distance
        self.wall_height = wall_height
        self.padstone_length = padstone_length
        self.padstone_width = padstone_width
        self.padstone_area = self.padstone_length * self.padstone_width

        # calculate the effective spread at mid-height of the wall
        spread_angle = 60 * math.pi / 180
        x_max = (self.wall_height / 2) /math.tan(spread_angle)
        x = min(x_max, self.near_edge_distance)
        y = min(x_max, self.far_edge_distance)
        self.effective_length = x + y + self.padstone_length
        self.effective_area = self.effective_length * self.padstone_width

        area_ratio = min(self.padstone_area / self.effective_area, 0.45)

        # calculate the strength enhancement factor
        beta_min = 1
        beta_max = min(1.25 + (self.near_edge_distance/(2*self.wall_height)), 1.5)
        beta_calculated = min((1 + 0.30 * (self.near_edge_distance / self.wall_height))*(1.5 - 1.1 * area_ratio),beta_max)
        self.beta = max(beta_min, beta_calculated)

        # calculate the capacity of the padstone
        self.padstone_capacity = self.beta * self.padstone_area * self.fk / self.gm * 10**-3 # in kN

    def padstone_effective_udl(self, point_load:float)->float:
        """Calculate the effective UDL at mid-height of the wall from the point load"""
        return point_load / (self.effective_length * 10**-3)

    # --- Print functions --- #
    def print_padstone_capacity(self)->None:
        print(f"The beta factor = {self.beta:.2f}")
        print(f"Wall effective length at mid-height = {self.effective_length:.2f}mm")
        print(f"Padstone design capacity = {self.padstone_capacity:.2f}kN")

