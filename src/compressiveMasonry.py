from structuralMasonry import Masonry
import math

class Wall(Masonry):
    """ compressive load capacity of a wall """
    def __init__(self, fb:float, fm:float, k:float, effective_height:float,
                 loaded_leaf_thickness:float, other_leaf_thickness:float, top_eccentricity:float = 5):
        super().__init__(fb, fm, k)

        self.effective_height: float = effective_height # [mm]
        self.loaded_leaf_thickness: float = loaded_leaf_thickness # [mm]
        self.other_leaf_thickness: float = other_leaf_thickness # [mm]
        self.top_eccentricity: float = top_eccentricity # [mm]

        # calculate the effective thickness of the wall
        self.effective_thickness: float = self.calc_effective_thickness(self.loaded_leaf_thickness,
                                                                 self.other_leaf_thickness)

        # calculate the eccentricity at the top of the wall
        self.einit: float = self.effective_height/450
        self.et: float = max(self.einit + self.top_eccentricity, 0.05 * self.loaded_leaf_thickness)

        # calculate the capacity reduction factor at the top of the wall
        self.top_slenderness_reduction: float = self.calc_strength_reduction(self.loaded_leaf_thickness,
                                                                             self.et)

        # calculate the eccentricity at mid-height of the wall
        self.em: float = self.einit + self.top_eccentricity / 2 # assuming a linear reduction
        self.ek: float = (0.002 * self.creep_coefficient * (self.effective_height/self.effective_thickness) *
                   (self.loaded_leaf_thickness * self.em)**0.50)
        self.emk: float = max(self.em + self.ek, 0.05 * self.loaded_leaf_thickness)

        # calculate the slenderness corrective factor at mid-height of the wall
        self.mid_slenderness_reduction: float = self.calc_slenderness_reduction(self.emk, self.loaded_leaf_thickness,
                                                                                self.effective_height,
                                                                                self.effective_thickness,
                                                                                self.fk, self.Em)

        # calculate the design vertical capacity of the wall [kN/m]

        self.design_vertical_capacity: float = self.calc_vertical_capacity(min(self.top_slenderness_reduction,
                                                                               self.mid_slenderness_reduction),
                                                                           self.fk, self.loaded_leaf_thickness, self.gm)

    # --- standard calculations ---
    def calc_effective_thickness(self,t1: float, t2: float)->float:
        """
        Calculates the effective stiffness of a cavity wall
        Parameters:
        -----------
        t1: float
            Thickness of the first leaf
        t2: float
            Thickness of the second leaf
        Returns:
        --------
            Effective thickness of the wall
        """
        return (t1**3 + t2**3)**(1/3)

    def calc_strength_reduction(self, wall_thickness: float, eccentricity: float)->float:
        """
        Strength reduction factor at the top or bottom of a wall
        Parameters:
        -----------
        wall_thickness: float
            Thickness of the wall leaf
        eccentricity: float
            Eccentricity of the load
        Returns:
        --------
            Strength reduction factor
        """
        return 1 - 2 * (eccentricity / wall_thickness)

    def calc_slenderness_reduction(self, eccentricity: float, wall_thickness:float,
                                   effective_height: float, effective_thickness: float,
                                   fk: float, Em: float)->float:
        """
        Calculates the reduction in strength at mid-hieght of a wall
        Parameters:
        -----------
        eccentricity: float
            The eccentricity at mid-height of the wall
        wall_thickness: float
            The thickness of the loaded leaf of the wall
        effective_height: float
            The effective height of the wall
        effective_thickness: float
            The effective thickness of the wall
        fk: float
            The characteristic strenth of the wall
        Em: float
            The modulus of elasticity of the wall
        Returns:
        --------
            Strength reduction factor at mid-height of the wall
        """
        A1 = 1 - 2 * (eccentricity / wall_thickness)
        lam = (effective_height / effective_thickness) * (fk / Em) ** 0.50
        u = (lam - 0.063) / (0.73 - 1.17 * (eccentricity / wall_thickness))
        return A1 * math.e ** (-u ** 2 / 2)

    def calc_vertical_capacity(self, reduction_factor: float, fk: float,
                               wall_thickness:float, partial_factor: float)->float:
        """
        Calculate the design vertical load capacity of a wall
        Parameters:
            reduction_factor: float
                The reduction factor at mid-height of the wall
            fk: float
                The characteristic strength of the wall
            wall_thickness: float
                The thickness of the loaded leaf of the wall
            partial_factor: float
                The partial factor for the design load
        Returns:
            Design vertical load capacity of the wall
        """

        return reduction_factor * fk * wall_thickness / partial_factor


    # --- standard print functions ---
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
        self.near_edge_distance: float = near_edge_distance
        self.far_edge_distance: float = far_edge_distance
        self.wall_height: float = wall_height
        self.padstone_length: float = padstone_length
        self.padstone_width: float = padstone_width
        self.padstone_area: float = self.padstone_length * self.padstone_width

        # calculate the effective spread at mid-height of the wa
        SPREAD_ANGLE: float = 60 * math.pi / 180
        x_max: float = (self.wall_height / 2) /math.tan(SPREAD_ANGLE)
        x: float = min(x_max, self.near_edge_distance)
        y: float = min(x_max, self.far_edge_distance)
        self.effective_length: float = x + y + self.padstone_length
        self.effective_area: float = self.effective_length * self.padstone_width

        area_ratio: float = min(self.padstone_area / self.effective_area, 0.45)

        # calculate the strength enhancement factor
        beta_min: float = 1.0
        beta_max: float  = min(1.25 + (self.near_edge_distance/(2*self.wall_height)), 1.5)
        beta_calculated: float = min((1 + 0.30 * (self.near_edge_distance / self.wall_height))*(1.5 - 1.1 * area_ratio),beta_max)
        self.beta: float  = max(beta_min, beta_calculated)

        # calculate the capacity of the padstone
        self.padstone_capacity: float = self.beta * self.padstone_area * self.fk / self.gm * 10**-3 # in kN

    # --- standard calculations ---
    def padstone_effective_udl(self, point_load:float)->float:
        """Calculate the effective UDL at mid-height of the wall from the point load"""
        return point_load / (self.effective_length * 10**-3)

    # --- Print functions --- #
    def print_padstone_capacity(self)->None:
        print(f"The beta factor = {self.beta:.2f}")
        print(f"Wall effective length at mid-height = {self.effective_length:.2f}mm")
        print(f"Padstone design capacity = {self.padstone_capacity:.2f}kN")

