import numpy as np
''' Amendments '''
'''
19/2/24 Top shear and top moments included within the analysis of retaining walls
'''
class earth():
    ''' class to module soil properties '''
    def __init__(self,psi, phi, delta, beta, gs):
	    # convert the angles to radians
        self.psi = np.radians(psi)
        self.phi = np.radians(phi)
        self.delta = np.radians(delta)
        self.beta = np.radians(beta)
        self.gs = gs
        a = (1/np.sin(self.psi)) * np.sin(self.psi - self.phi)
        b = np.sqrt(np.sin(self.psi + self.delta))
        c = np.sqrt((np.sin(self.phi+self.delta)*np.sin(self.phi-self.beta))/np.sin(self.psi - self.beta))
        d = np.sqrt(np.sin(self.psi - self.delta))
        e = np.sqrt((np.sin(self.phi+self.delta)*np.sin(self.phi+self.beta))/np.sin(self.psi - self.beta))
        self.kp = (a/(d-e))**2
        self.ka = (a / (b+c))**2

class retainingWall(earth):
    def __init__(self,psi,phi,delta,beta,gs, q, hs, ts, tb, lb, toe, ws, baseAngle, topShear, topMoment, wAdd):
        super().__init__(psi,phi,delta, beta,gs)
        
        self.mu = np.tan(np.radians(baseAngle))

        # calculate the earth effects on the wall
        self.fe = self.ka * (gs *  hs**2/2 + q*hs)*np.cos(self.delta) + topShear
        self.me = self.ka * (gs*hs**3/6 + q*hs**2/2)*np.cos(self.delta) + topShear * hs + topMoment
        self.te = self.fe * np.tan(self.delta)

        # calculate the gravity effects for the wall
        self.stem_weight = ws * ts * (hs - tb)
        self.base_weight = 24 * lb * tb
        self.earth_weight = (lb - toe - ts) * (hs - tb) * self.gs
        self.wAdd = wAdd
        self.total_weight = self.stem_weight + self.base_weight + self.earth_weight + self.te + self.wAdd
        # calculate the sliding resistance for the wall
        self.slidingResistance = self.total_weight * self.mu

        # calculate the lever arms from the toe of the wall for 
        # overturning calculations
        self.la_stem = toe + ts / 2
        self.la_base = lb / 2
        self.la_soil = toe + ts + (lb-toe-ts)/2
        self.la_fric = lb
        self.la_wAdd = toe  + 0.215 / 2 # lever arm to additional load added to the top of the wall
        self.gs = gs
        self.q = q
        self.topShear = topShear
        self.topMoment = topMoment

        # calculate the overturning resistance
        self.overturningResistance = self.la_stem * self.stem_weight + self.la_base * self.base_weight + self.la_soil * self.earth_weight + self.la_fric * self.te + self.wAdd * self.la_wAdd

        # calculate the factors of safety for overturning and sliding
        self.fosSliding = self.slidingResistance / self.fe
        self.fosOverturning = self.overturningResistance / self.me

        # calculate the bearing pressures under the base

        # calculate the eccentricity from the toe xBar
        netMoment = self.overturningResistance - self.me
        xBar = netMoment / self.total_weight

        # calculate the eccentricity from the centre line of the base ecc
        self.ecc = lb / 2 - xBar
        eccAbs = np.absolute(self.ecc)

        # check to see if ecc is outside of the kern of the base and
        # calculate the bearing pressure as applicable
        self.kern = lb / 6
        
        self.pmax = 0 # set a conservative default value

        if eccAbs > self.kern: # triangular pressure assumed
            if self.ecc > 0:
                self.pmax = 2 * self.total_weight / (3 * xBar)
                self.pmin = 0
            else:
                self.pmax = 2 * self.total_weight / (3 * (lb - xBar))
                self.pmin = 0
        else:
           self.pmax = (self.total_weight / lb) + (self.total_weight * eccAbs)/(lb**2 / 6)
           self.pmin = (self.total_weight / lb) - (self.total_weight * eccAbs)/(lb**2 / 6)

        # calculate the ULS bending moments at the base of the stem
        self.stemHeight = hs - tb

        # calculate the earth effects on the wall
        self.ulsfe = self.ka * (self.gs*1.35 *  self.stemHeight**2/2 + self.q*1.5*self.stemHeight)*np.cos(self.delta) + topShear * 1.5
        self.ulsme = self.ka * (self.gs*1.35 * self.stemHeight**3 / 6 + self.q*1.5*self.stemHeight**2/2)*np.cos(self.delta) + topShear * 1.5 * self.stemHeight + topMoment * 1.5
        self.ulste = self.fe * np.tan(self.delta)

    def forceAtDepth(self, z):
        self.ulsfe = self.ka * (self.gs*1.35 *  z**2/2 + self.q*1.5*z)*np.cos(self.delta) + self.topShear * 1.5
        self.ulsme = self.ka * (self.gs*1.35*z**3 / 6 + self.q*1.5*z**2/2)*np.cos(self.delta) + self.topShear * 1.5 * z + self.topMoment * 1.5
        self.ulste = self.fe * np.tan(self.delta)
        print(f"ULS Moment at depth {z:.3f}m = {self.ulsme:.2f} kNm")
        print(f"ULS Shear at depth {z:.3f}m = {self.ulsfe:.2f} kN")
        return self.ulsfe, self.ulsme
        

    def printStability(self):

        # Print the analysis for the retaining wall

        print(f"Coefficient of active earth pressure = {self.ka:.3f}")
        print(f"The horizontal shear = {self.fe:.2f} kN")
        print(f"The overturning moment = {self.me:.2f} kNm")
        print(f"The vertical shear on the wall back = {self.te:.2f} kN")
        print(f"The total weight of the wall = {self.total_weight:.2f} kN")
        print(f"The coefficient of base friction = {self.mu:.2f}")
        print(f"The sliding resistance = {self.slidingResistance:.2f} kN")
        print(f"The overturning resistance = {self.overturningResistance:.2f} kNm")
        print("")
        print(f"Eccentricity of the resultant load from \nthe center line of the base towards the toe = {self.ecc:.2f}m")
        print(f"The kern of the section = {self.kern:.2f}m")
        print("")
        print(f"The factor of safety against sliding = {self.fosSliding:.2f}")
        print(f"The factor of safety against overturning = {self.fosOverturning:.2f}")
        print(f"The contact pressure under the base at the toe = {self.pmax:.2f} kPa")
        print(f"The contact pressure under the base at the heel = {self.pmin:.2f} kPa")
        print("")
        print(f"The ULS bending moment at the base of the stem = {self.ulsme:.2f} kNm")
        print(f"The ULS shear force at the base of the stem = {self.ulsfe:.2f} kN")
    
    def masonryCheck(self,z,twall, fkx, gWall):
        ''' check the capacity of a masonry stem'''
        w = z * twall * gWall
        gd = w * 10**3 / (twall * 1000)
        fxd = fkx / 2.7
        zMod = 1000 * (twall*1000)**2 / 6
        mCap = fxd * zMod * 10**-6
        print(f"The moment capacity of the {twall:.2f}m thick stem = {mCap:.2f}kNm")
        self.forceAtDepth(z)
        if self.ulsme > mCap:
            print(f'{twall:.2f}m thick stem fails at this depth')
        else:
            print(f'{twall:.2f}m thick stem is adequate to {z:.2f}m below retained ground level')

