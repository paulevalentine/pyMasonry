import sympy as sp

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

        self.fk = self.k * self.fb ** self.a * self.fm ** self.b

        # calculate the instant elastic modulus
        self.Em = 1000 * self.fk

        #todo add shear strength calculation in basic form
        #todo add flexural strength calculation in basic form