import sympy as sp
from IPython.core.display import Markdown
from IPython.core.display_functions import display
globals()['Markdown'] = Markdown

class Masonry():
    def __init__(self,fb:float,fm:float,masonry_type:str)->None:

        """ define the main properties of the masonry
        Parameters:
            fb unit strength of the brick (MPa)
            fm compressive strength of masonry (MPa)
            masonry_type is either 'brick' or 'block' and used to generate the K value

            Cat I units assumed in execution class 2 is assumed for partial factors
        """

        self.fb: float = fb # normalised compressive strength of the unit
        self.fm: float = fm # compressive strength of the mortar

        if masonry_type=="brick":
            self.k: float = 0.50
        elif masonry_type=="block":
            self.k: float = 0.75
        else:
            self.k: float = 0.40 # lowest value in the table

        self.a: float = 0.70  # alpha for standard mortar beds
        self.b: float = 0.30  # beta factor
        self.gm: float = 2.7  # partial factor on material strength
        self.creep_coefficient: float = 1.5

        # calculate the characteristic compressive strength of the masonry
        self.fk: float = self.calc_characteristic_strength()
        # calculate the instant elastic modulus
        self.Em: float = self.calc_elastic_modulus()

    def calc_characteristic_strength(self) -> float:

        fk, k,fb, a, fm, b = sp.symbols('f_k, K,  f_b, alpha, f_m, beta')
        char_strength = sp.Eq(fk, k * fb**a * fm**b)
        display(Markdown("Calculate characteristic compressive strength of masonry (MPa):"))
        val = char_strength.subs({k:self.k, fb:self.fb, a:self.a, fm:self.fm, b:self.b}).evalf(3)
        display(char_strength)
        display(val)
        return val.evalf().rhs

    def calc_elastic_modulus(self) -> float:

        display(Markdown("Calculate the elastic modulus for the masonry (MPa):"))
        Em, fk = sp.symbols("E_m, f_k")
        elastic_modulus = sp.Eq(Em, 1000 * fk)
        val = elastic_modulus.subs({fk:self.fk}).evalf(3)
        display(elastic_modulus)
        display(val)

        return val.evalf().rhs

        #todo add shear strength calculation in basic form
        #todo add flexural strength calculation in basic form