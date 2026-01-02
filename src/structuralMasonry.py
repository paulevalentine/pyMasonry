import sympy as sp
from IPython.core.display_functions import display


class Masonry():
    def __init__(self,fb:float,fm:float,k:float)->None:
        """ define the main properties of the masonry """
        self.fb: float = fb # normalised compressive strenght of the unit
        self.fm: float = fm # comressive strength of the mortar
        self.k: float = k # factor for compressive strength
        '''
        K = 0.75 typically for block Group 1 units in general purpose mortar
        K = 0.50 typically for clay brick group 2 units in general purpose mortar
        '''
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
        print("Calculate characteristic compressive strength (MPa):")
        val = char_strength.subs({k:self.k, fb:self.fb, a:self.a, fm:self.fm, b:self.b}).evalf(3)
        display(char_strength)
        display(val)
        return val.evalf().rhs

    def calc_elastic_modulus(self) -> float:

        print("Calculate the elastic modulus for the masonry (MPa):")
        Em, fk = sp.symbols("E_m, f_k")
        elastic_modulus = sp.Eq(Em, 1000 * fk)
        val = elastic_modulus.subs({fk:self.fk}).evalf(3)
        display(elastic_modulus)
        display(val)

        return val.evalf().rhs

        #todo add shear strength calculation in basic form
        #todo add flexural strength calculation in basic form