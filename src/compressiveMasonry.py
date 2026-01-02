from IPython.core.display_functions import display

from structuralMasonry import Masonry
import math
import sympy as sp
sp.init_printing(order='none')

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
        self.effective_thickness: float = self.calc_effective_thickness()

        # calculate the eccentricity at the top of the wall
        print("Calculate the eccentricity at the top of the wall (mm):")
        heff, ec, etop = sp.symbols('h_eff, e_c, e_top')
        ec_eq = sp.Eq(ec, etop + heff/450)
        ec_val = ec_eq.subs({heff:self.effective_height, etop:self.top_eccentricity}).evalf(4)
        display(ec_eq, ec_val)
        self.ec_top: float = ec_val.rhs
        self.einit: float = self.effective_height / 450

        emin, t = sp.symbols('e_min, t')
        emin_eq = sp.Eq(emin, 0.05 * t)
        emin_val = emin_eq.subs({t:self.loaded_leaf_thickness}).evalf(4)
        display(emin_eq, emin_val)

        self.et: float = max(self.ec_top, 0.05 * self.loaded_leaf_thickness)

        et = sp.symbols('e_t')
        et_eq = sp.Eq(et, self.et)
        display(et_eq)

        # calculate the capacity reduction factor at the top of the wall
        self.top_slenderness_reduction: float = self.calc_strength_reduction(self.loaded_leaf_thickness,
                                                                             self.et)

        # calculate the eccentricity at mid-height of the wall
        print("Calculate the eccentricity at mid-height of the wall (mm):")

        em, etop, heff = sp.symbols('e_m, e_top, h_ef')
        em_eq = sp.Eq(em,heff/450 + etop/2)
        em_val = em_eq.subs({heff:self.effective_height, etop:self.top_eccentricity}).evalf(4)
        display(em_eq, em_val)
        self.em: float = em_val.rhs

        ek, teff, phi_creep, t = sp.symbols('e_k, t_ef, phi_creep, t')
        ek_eq = sp.Eq(ek, 0.002 * phi_creep * (heff/teff) * (t*em)**sp.Rational(1,2))
        ek_val = ek_eq.subs({teff:self.effective_thickness, phi_creep:self.creep_coefficient,
                             t:self.loaded_leaf_thickness, heff:self.effective_height,
                             em:self.em}).evalf(4)
        display(ek_eq, ek_val)

        self.ek: float = ek_val.rhs
        emk = sp.symbols('e_mk')
        emk_eq = sp.Eq(emk, ek + em)
        emk_val = emk_eq.subs({ek:ek_val.rhs, em:em_val.rhs}).evalf(4)
        display(emk_eq, emk_val)

        display(emin_val)

        self.emk: float = max(self.em + self.ek, 0.05 * self.loaded_leaf_thickness)

        # calculate the slenderness corrective factor at mid-height of the wall
        self.mid_slenderness_reduction: float = self.calc_slenderness_reduction()

        # calculate the design vertical capacity of the wall [kN/m]

        self.design_vertical_capacity: float = self.calc_vertical_capacity(min(self.top_slenderness_reduction,
                                                                               self.mid_slenderness_reduction),
                                                                           self.fk, self.loaded_leaf_thickness, self.gm)



    # --- standard calculations ---
    def calc_effective_thickness(self)->float:

        print("Calculate the effective thickness of the wall (mm):")
        teff, t1, t2 = sp.symbols("t_e, t_1, t_2")
        effective_thickness = sp.Eq(teff, (t1**3 + t2**3)**sp.Rational(1,3))
        val = effective_thickness.subs({t1: self.loaded_leaf_thickness, t2: self.other_leaf_thickness}).evalf(4)
        display(effective_thickness, val)
        return val.evalf().rhs

    def calc_strength_reduction(self, wall_thickness: float, eccentricity: float)->float:
        print(f"Calculate the strength reduction factor for the top of the wall:")

        phi, et, t = sp.symbols("phi, e_t, t")
        strength_reduction = sp.Eq(phi, 1 - 2 * (et / t))
        val = strength_reduction.subs({et: self.et, t: self.loaded_leaf_thickness}).evalf(2)
        display(strength_reduction, val)
        return val.evalf().rhs

    def calc_slenderness_reduction(self)->float:
        print("Calculate the strength reduction factor at mid-height of the wall:")
        A1, emk, t = sp.symbols("A_1, e_mk, t")
        eqA1 = sp.Eq(A1, 1 - 2 * (emk/t))
        display(eqA1)
        lam, hef, teff, fk, Em = sp.symbols("lambda, h_ef, t_eff, f_k, E_m")
        eqlam = sp.Eq(lam, (hef / teff) * (fk / Em) ** sp.Rational(1, 2))
        display(eqlam)
        u = sp.symbols("u")
        equ = sp.Eq(u, (lam - 0.063) / (0.73 - 1.17 * (emk /t)))
        display(equ)
        phi = sp.symbols("phi")
        eqphi = sp.Eq(phi, A1 * sp.exp(-u **2 / 2))
        display(eqphi)

        val_A1 = eqA1.subs({emk:self.emk, t:self.loaded_leaf_thickness}).evalf(3)
        display(val_A1)
        val_lam = eqlam.subs({hef:self.effective_height, teff:self.effective_thickness,
                             fk:self.fk, Em:self.Em}).evalf(3)
        display(val_lam)
        val_u = equ.subs({lam:val_lam.rhs, emk:self.emk, t:self.loaded_leaf_thickness}).evalf(3)
        display(val_u)

        val_phi = eqphi.subs({A1:val_A1.rhs, u:val_u.rhs}).evalf(2)
        display(val_phi)

        return val_phi.rhs

    def calc_vertical_capacity(self, reduction_factor: float, fk: float,
                               wall_thickness:float, partial_factor: float)->float:
        print("Calculate the design vertical load capacity of the wall (kN/m):")
        DVLR, phi_min, fk, gm, t = sp.symbols('DVLR, phi_min, f_k, gamma_m. t')
        eq_DVLR = sp.Eq(DVLR,phi_min * fk * t / gm)

        phi_min_eq = sp.Eq(phi_min, min(self.top_slenderness_reduction, self.mid_slenderness_reduction,
                                        self.mid_slenderness_reduction)).evalf(2)
        display(phi_min_eq)
        val_DVLR = eq_DVLR.subs({phi_min:phi_min_eq.evalf().rhs, fk:self.fk, t:self.loaded_leaf_thickness, gm:self.gm}).evalf(4)
        display(eq_DVLR, val_DVLR)
        return val_DVLR.rhs


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

