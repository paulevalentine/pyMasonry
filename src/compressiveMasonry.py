from structuralMasonry import Masonry
import math

class Wall(Masonry):
    """ compressive load capacity of a wall """
    def __init__(self, fb:float, fm:float, k:float, effective_height:float,
                 loaded_leaf_thickness:float, other_leaf_thickness:float, mid_height_eccentricity:float = 5):
        super().__init__(fb, fm, k)
        self.effective_height = effective_height # [mm]
        self.loaded_leaf_thickness = loaded_leaf_thickness # [mm]
        self.other_leaf_thickness = other_leaf_thickness # [mm]
        self.mid_height_eccentricity = mid_height_eccentricity # [mm]

        # calculate the effective thickness of the wall
        self.effective_thickness = (self.loaded_leaf_thickness**3 + self.other_leaf_thickness**3)**(1/3)
        # calculate the mid-height eccentricity
        self.einit = self.effective_height/450
        self.emk = self.einit + self.mid_height_eccentricity

        # calculate the slenderness corrective factor at mid-height of the wall
        A1 = 1 - 2 * (self.emk / self.loaded_leaf_thickness)
        lam = (self.effective_height / self.effective_thickness) * (self.fk / self.Em) ** 0.50
        u = (lam - 0.063) / (0.73 - 1.17 * (self.emk / self.loaded_leaf_thickness))
        self.mid_slenderness_reduction = A1 * math.e ** (-u ** 2 / 2)

        # todo add a means of assessing top and bottom eccentricity calculations
        # calculate the design vertical capacity of the wall [kN/m]
        self.design_vertical_capacity = (self.mid_slenderness_reduction *
                                         self.fk * self.loaded_leaf_thickness / self.gm)

    # Return the design vertical load capacity of the wall
    def print_design_vertical_capacity(self)->None:
        print(f"hef = {self.effective_height:.2f}mm, twall = {self.loaded_leaf_thickness:.2f}mm")
        print(f"fk = {self.fk:.2f}MPa : Phi = {self.mid_slenderness_reduction:.2f}")
        print(f"Design vertical load capacity of the wall = "
              f"{self.design_vertical_capacity:.2f} kN/m")
