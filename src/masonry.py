# module for the design of masonry
# 02/06/243 - need to supply K factor added.
# modules to import
import math

class masonry():
    """ functions for the design of masonry elements """

    def __init__(self, fm, fb, K):
        """ define the main properties of the masonry """
        '''
        K = 0.75 typically for block Group 1 units in general purpose mortar
        K = 0.50 typically for clay brick group 2 units in general purpose mortar
        '''
        self.fm = fm # unit compressive strength of mortar
        self.fb = fb # unit compressive strength of masonry
        self.K = K
        self.a = 0.70 # alpha for standard mortar beds
        self.b = 0.30 # beta factor
        self.gm = 2.7 # partial factor on material strength
        # calculate the characteristic compressive strength of the masonry
        self.fk = self.K * self.fb**self.a * self.fm**self.b

class masonry_wall(masonry):
    """" load capacity of a wall """

    def __init__(self, fm,fb, hef, twall, touter, ec, K):
        super().__init__(fm,fb, K)
        self.hef = hef # effective height of the wall
        self.twall = twall # width of the loaded leaf
        self.touter = touter # width of any outer leaf

        """ calculate the mid-height slenderness reduction factor """
        teff = (self.touter**3 + self.twall**3)**(1/3)
        einit = self.hef/450 # mid height eccentricity 
        emk = ec + einit # total eccentricity at mid height
        A1 = 1 - 2*(emk/self.twall)
        lam = (self.hef/teff)*(1/(1000))**0.50
        u = (lam - 0.063)/(0.73 - 1.17*(emk/self.twall))
        self.phi = A1 * math.e**(-u**2/2)

        """ calculate the design vertical load capacity of a wall """
        self.wall_capacity = self.phi * self.fk * self.twall / self.gm

    def padstone_capacity(self, a1, a2, hc, lb, tb): # all in mm
        """ calculate the local capacity of an axially loaded padstone """
        '''
        a1 is the disstance from the edge of the padstone to the end of the masonry on side 1
        a2 is the disstance from the edge of the padstone to the end of the masonry on side 2
        '''
        Ab = lb * tb # area of the bearing
        # calculate half the spread to the midheightof wall
        spread = hc / (2 * math.tan(60 * math.pi / 180))
        if a2 < spread and a2 < spread:
            lef = a1 + lb + a2
        elif a1 < spread and a2 >= spread:
            lef = a1 + spread + lb
        elif a1 >= spread and a2 < spread:
            lef = a2 + spread + lb
        else:
            lef = lb + 2 * spread
        print(f"Effective length at mid height = {lef:.2f} mm")
        print(f"The possible spread through the wall = {spread:.2f}")
        Aef = lef * lb
        b1 = max((1 + 0.30 * (a1/hc)) * (1.5 - 1.1 * Ab / Aef),1)
        b2 = min(1.25 + a1 / (2 * hc), 1.5)
        bf = min(b1, b2)
        padCap =  lb * tb * bf * self.fk * 10**-3 / self.gm
        capacity_line_load = padCap / (lef*10**-3)
        return bf, capacity_line_load, padCap

    def flexural_capacity(self, fxk, gd):
        """ Calculation for the flexural capacity of the wall """
        z_wall = 1000 * self.twall**2 / 6
        fxk_mod = fxk + self.gm * gd
        mCap = fxk_mod * z_wall *10**-6 / self.gm # capacity in kNm.
        return mCap

    def calcEffectiveUDL(self, wallLength, windowOpenings, udl):
        ''' function to increase the udl for openigns'''
        netLength = wallLength - windowOpenings
        factor = wallLength / netLength
        self.effectiveUDL = factor * udl


# ---- printer functions ----

    def getVertCapacity(self):
        ''' function to print the results of the vertical load capacity '''
        print(f'Characteristic compressive strength = {self.fk:.2f}MPa')
        print(f'The slenderness correctof factor = {self.phi:.2f}')
        print(f"The design vertical load capacity of the wall = {self.wall_capacity:.2f}kN/m")

    def getPadstoneCapacity(self, a1, a2, hc, lb, tb):
        ''' function to print a padstone capacity '''
        bf, capacity_line_load, padCap = self.padstone_capacity(a1, a2, hc, lb, tb)
        print(f"Capacity enhancment factor = {bf:.2f}")
        print(f"Padstone capacity = {padCap:.2f} kN")
        #print(f"Quasi line load at mid-height from full capacity = {capacity_line_load:.2f}kN/m")

    def getEffectiveUDL (self, wallLength, windowOpenings, udl):
        ''' function to print the walls effective udl '''
        self.calcEffectiveUDL(wallLength, windowOpenings, udl)
        print(f"The effective udl on the wall = {self.effectiveUDL:.2f} kN/m")

class masonry_pier(masonry):

    def __init__(self, fm,fb, hef, lwall, twall, ec, K):
        super().__init__(fm,fb,K)

        
        self.hef = hef # effective height of the wall
        self.twall = twall # width of the pier
        self.lwall = lwall # length of the pier

        """ calculate the mid-height slenderness reduction factor """
        einit = self.hef/450 # mid height eccentricity
        emk = ec + einit  # total eccentricity at mid height
        A1 = 1 - 2*(emk/self.twall)
        lam = (self.hef/twall)*(1/(1000))**0.50
        u = (lam - 0.063)/(0.73 - 1.17*(emk/self.twall))
        self.phi = A1 * math.e**(-u**2/2)

        """ calculate the design vertical load capacity of a wall """
        self.pier_capacity = self.phi * self.fk * self.twall * self.lwall * 10**-3 / self.gm

    def getPierCapacity(self):
        print(f"The characteristic compressive strength of the masonry = {self.fk:.2f} MPa")
        print(f"The slenderness reduction factor for the pier = {self.phi:.2f}")
        print(f"The axial capacity of the {self.twall} x {self.lwall} mm pier is {self.pier_capacity:.2f} kN")




