import math
import fluids
import numpy as np

import config_earth

"""
radiation3.py solves for radiation due to the enviorment for a particular datetime and altitude.
"""

class Radiation:
    # Constants
    I0 = 1358               # Direct Solar Radiation Level
    e = 0.016708            # Eccentricity of Earth's Orbit
    P0 = 101325             # Standard Atmospheric Pressure at Sea Level
    cloudElev = 3000        # (m)
    cloudFrac = 0.0         # Percent cloud coverage [0,1]
    cloudAlbedo = .65       # [0,1]
    albedoGround = .2       # Ground albedo [0,1]
    tGround = 293           # (K) Temperature of Ground
    emissGround = .95       # [0,1]
    SB = 5.670373E-8        # Stefan Boltzman Constant
    RE = 6371000            # (m) Radius of Earth
    radRef= .1              # [0,1] Balloon Reflectivity
    radTrans = .1           # [0,1] Balloon Transmitivity


    start_coord = config_earth.GNC['start_coord']
    t = config_earth.GNC['start_time']
    lat = math.radians(start_coord['lat'])
    Ls = t.timetuple().tm_yday
    d = config_earth.balloon_properties['d']
    #emissEnv = config_earth.balloon_properties['emissEnv']
    absEnvIR = config_earth.balloon_properties['absEnvIR']
    absEnv = config_earth.balloon_properties['absEnv']

    projArea = 0.25*math.pi*d*d
    surfArea = math.pi*d*d

    def setAbs(self,abs):
        absEnv = abs

    def setAbsIR(self,absIR):
        absEnv = absIR

    def getTemp(self, el):
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        return atm.T #+ 20

    def getPressure(self, el):
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        return atm.P

    def getDensity(self, el):
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        return atm.rho

    def getGravity(self, el):
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        return atm.g

    def get_SI0(self):
        """ Incident solar radiation above Earth's atmosphere (W/m^2)

        :returns: The incident solar radiation above Earths atm (W/m^2)
        :rtype: float
        """

        f = 2*math.pi*Radiation.Ls/365 #true anomaly
        e2 = pow(((1.+Radiation.e)/(1.-Radiation.e)),2) -1.
        return Radiation.I0*(1.+0.5*e2*math.cos(f))

    def get_declination(self):
        """Expression from http://en.wikipedia.org/wiki/Position_of_the_Sun

        :returns: Approximate solar declination (rad)
        :rtype: float
        """

        return -.4091*math.cos(2*math.pi*(Radiation.Ls+10)/365)

    def get_zenith(self, t, coord):
        """ Calculates solar zenith angle

        :param t: Lattitude (rad)
        :type t: Datetime
        :param coord: Solar Hour Angle (rad)
        :type coord: dict
        :returns: The approximate solar zenith angle (rad)
        :rtype: float
        """

        solpos = fluids.solar_position(t, coord["lat"], coord["lon"], Z = coord["alt"])
        zen = math.radians(solpos[0]) #get apparent zenith

        # For determining adjusted zenith at elevation
        # https://github.com/KosherJava/zmanim/blob/master/src/main/java/com/kosherjava/zmanim/util/AstronomicalCalculator.java#L176
        refraction = 4.478885263888294 / 60.
        solarRadius = 16 / 60.
        earthRadius = 6356.9; # in KM
        elevationAdjustment = math.acos(earthRadius / (earthRadius + (coord["alt"]/1000.)));


        adjusted_zen =zen - elevationAdjustment + math.radians(solarRadius + refraction)
        return adjusted_zen

    def get_air_mass(self,zen, el):
        """Air Mass at elevation

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The approximate air mass (unitless)
        :rtype: float
        """

        p = self.getPressure(el) #pressure at current elevation
        am = (p/Radiation.P0)*(math.sqrt(1229 + pow((614*math.cos(zen)),2))-614*math.cos(zen))
        return am

    def get_trans_atm(self,zen,el):
        """get zenith angle

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The atmospheric trasmittance (unitless)
        :rtype: float
        """

        am = self.get_air_mass(zen, el)
        trans = 0.5*(math.exp(-0.65*am) + math.exp(-0.095*am))
        return trans

    def get_direct_SI(self,zen,el):
        """Calculates Direct Solar Radiation

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: Tntensity of the direct solar radiation (W/m^2)
        :rtype: float
        """

        SI0 = self.get_SI0()
        trans = self.get_trans_atm(zen, el)
        if zen > math.pi/2 :
            direct_SI = 0
        else:
            direct_SI = trans*SI0

        return direct_SI

    def get_diffuse_SI(self,zen,el):
        """Calculates Diffuse Solar Radiation from sky

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The intensity of the diffuse solar radiation from the sky (W/m^2)
        :rtype: float
        """

        if(zen > math.pi/2.):
            return 0.0
        SI0 = self.get_SI0()
        trans = self.get_trans_atm(zen, el)
        if el < Radiation.cloudElev:
            return (1-Radiation.cloudFrac)*0.5*SI0*math.sin(math.pi/2.-zen)*(1.-trans)/(1-1.4*math.log(trans))
        else:
            return 0.5*SI0*math.sin(math.pi/2.-zen)*(1.-trans)/(1-1.4*math.log(trans))

    def get_reflected_SI(self,zen,el):
        """Calculates Reflected Solar Radiation from from the Earth's Surface

        :param zen: Solar Angle (rad)
        :type zen: float
        :param el: Elevation (m)
        :type el: float
        :returns: The intensity solar radiation reflected by the Earth (W/m^2)
        :rtype: float
        """

        incident_SI = self.get_SI0()
        tau_atm = self.get_trans_atm(zen,el)
        if el < Radiation.cloudElev:
            albedo = (1.-Radiation.cloudFrac)*Radiation.albedoGround;
        else:
            albedo = (1.-Radiation.cloudFrac)*(1-Radiation.cloudFrac)*Radiation.albedoGround + Radiation.cloudAlbedo*Radiation.cloudFrac

        return albedo*tau_atm*incident_SI*math.sin(math.pi/2.-zen)

    def get_earth_IR(self,el):
        """Calculates Infared Radiation emitted from Earth's surface

        :param el: Elevation (m)
        :type el: float
        :returns: Intensity of IR radiation emitted from earth (W/m^2)
        :rtype: float
        """
        p = self.getPressure(el)#pressure at current elevation
        IR_trans = 1.716-0.5*(math.exp(-0.65*p/Radiation.P0) + math.exp(-0.095*p/Radiation.P0))
        if el < Radiation.cloudElev:
            tEarth = Radiation.tGround
        else:
            clouds = fluids.atmosphere.ATMOSPHERE_1976(Radiation.cloudElev)
            tEarth = Radiation.tGround*(1.-Radiation.cloudFrac) + clouds.T*Radiation.cloudFrac
        return IR_trans*Radiation.emissGround*Radiation.SB*pow(tEarth,4)

    def get_sky_IR(self,el):
        """Calculates Infared Radiation emitted the from Sky

        :param el: Elevation (m)
        :type el: float
        :returns: Intensity of IR radiation emitted from sky (W/m^2)
        :rtype: float
        """

        return np.fmax(-0.03*el+300.,50.0)

    def get_rad_total(self,datetime,coord):
        """Total Radiation as a function of elevation, time of day, and balloon surface area

        """

        zen = self.get_zenith(datetime, coord)
        el = coord["alt"]

        radRef = Radiation.radRef + Radiation.radRef*Radiation.radRef +  Radiation.radRef*Radiation.radRef*Radiation.radRef
        totAbs = Radiation.absEnv # + Radiation.absEnv*Radiation.radTrans + Radiation.absEnv*Radiation.radTrans*radRef

        hca = math.asin(Radiation.RE/(Radiation.RE+el)) #half cone angle
        vf = 0.5*(1. - math.cos(hca)) #viewfactor

        direct_I = self.get_direct_SI(zen, el)
        power_direct = direct_I*totAbs*Radiation.projArea

        diffuse_I = self.get_diffuse_SI(zen, el)
        power_diffuse = diffuse_I*totAbs*(1.-vf)*Radiation.surfArea

        reflected_I = self.get_reflected_SI(zen, el)
        power_reflected = reflected_I*totAbs*vf*Radiation.surfArea

        earth_IR = self.get_earth_IR(el)
        power_earth_IR = earth_IR*Radiation.absEnvIR*vf*Radiation.surfArea

        sky_IR = self.get_sky_IR(el)
        power_sky_IR = sky_IR*Radiation.absEnvIR*(1.-vf)*Radiation.surfArea

        rad_tot_bal = power_direct + power_diffuse + power_reflected + power_earth_IR + power_sky_IR

        return rad_tot_bal
