from IPython.core.display import Markdown
from IPython.core.display_functions import display
import sympy as sp

from structuralMasonry import Masonry

class Flexural(Masonry):
    def __init__(self, fbi: float, fmi: float, masonry_type_i:str, fbo: float,
                 fmo: float, masonry_type_o: str, masonry_ref_i: str ="",masonry_ref_o: str ="" ):
        super().__init__(fbi, fmi, masonry_type_i, masonry_ref_i)

        self.outer_leaf = Masonry(fbo, fmo, masonry_type_o, masonry_ref_o)

        self.masonry_ref_i: str = masonry_ref_i
        self.masonry_ref_o: str = masonry_ref_o

        self.inner_wall_thickness = 100
        self.outer_wall_thickness = 102.5

    def total_cracked_section(self, inner_load:float, outer_load:float)->float:
        inner_leaf_cracked = self.calc_cracked_section(self.inner_wall_thickness,
                                                       self.fk, inner_load, self.masonry_ref_i)
        outer_leaf_cracked = self.calc_cracked_section(self.outer_wall_thickness,
                                                       self.outer_leaf.fk, outer_load, self.masonry_ref_o)

        display(Markdown("Total cracked section:"))
        Mcrt, Mcri, Mcro = sp.symbols("M_ct M_ci M_co")
        Mct_eq = sp.Eq(Mcrt, Mcri + Mcro)
        Mct_val = Mct_eq.subs({Mcri:inner_leaf_cracked, Mcro:outer_leaf_cracked}).evalf(3)
        display(Mct_eq, Mct_val)
        return Mct_val.rhs

    def calc_cracked_section(self, thickness:float, char_strength:float, load:float, ref: str)->float:
        display(Markdown(f"{ref}"))
        fd, fk, gm, W, ws = sp.symbols('f_d f_k gamma_m, W, w_s')
        fd_eq = sp.Eq(fd, 1.1*fk/gm)
        fd_val = fd_eq.subs({fk:char_strength, gm:self.gm}).evalf(4)
        #display(Markdown("Calculate the design strength at base / head of wall (MPa)"))
        display(fd_eq, fd_val)

        #display(Markdown("Calculate stress block width (mm)"))
        ws_eq = sp.Eq(ws, W/fd)
        ws_val = ws_eq.subs({fd:fd_val.rhs, W: load}).evalf(3)
        display(ws_eq, ws_val)

        #display(Markdown("Calculate the lever arm (mm)"))
        la, t = sp.symbols('l_a t')
        la_eq = sp.Eq(la, (t-ws)/2)
        la_val = la_eq.subs({t:thickness, ws:ws_val.rhs}).evalf(3)
        display(la_eq, la_val)

        #display(Markdown("Calculate the cracked moment of resistance (kNm):"))
        Mcr =sp.symbols('M_cr')
        Mcr_eq = sp.Eq(Mcr, W * la * 10**-3)
        Mcr_val = Mcr_eq.subs({W:load, la:la_val.rhs}).evalf(3)
        display(Mcr_eq, Mcr_val)

        return Mcr_val.rhs



