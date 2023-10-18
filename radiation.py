"""
radiation solves for radiation due to the enviorment for a particular datetime and altitude.

"""

import math
import fluids
import numpy as np

import config_earth

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

    start_coord = config_earth.simulation['start_coord']
    t = config_earth.simulation['start_time']
    lat = math.radians(start_coord['lat'])
    Ls = t.timetuple().tm_yday
    d = config_earth.balloon_properties['d']
    emissEnv = config_earth.balloon_properties['emissEnv']
    absEnv = config_earth.balloon_properties['absEnv']

    projArea = 0.25*math.pi*d*d
    surfArea = math.pi*d*d

    def getTemp(self, el):
        atm = fluids.atmosphere.ATMOSPHERE_1976(el)
        return atm.T

    def getTempForecast(self, coord):
        r""" Looks up the forecast temperature at the current coordinate and altitude

        .. important:: TODO. This function is not operational yet. 

        :param coord: current coordinate
        :type coord: dict
        :returns: atmospheric temperature (k)
        :rtype: float

        """
        return temp

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
        r""" Incident solar radiation above Earth's atmosphere (W/m^2)

        .. math:: I_{sun,0}= I_0 \cdot [1+0.5(\frac{1+e}{1-e})^2-1) \cdot cos(f)]


        :returns: The incident solar radiation above Earths atm (W/m^2)
        :rtype: float

        """

        f = 2*math.pi*Radiation.Ls/365 #true anomaly
        e2 = pow(((1.+Radiation.e)/(1.-Radiation.e)),2) -1.
        return Radiation.I0*(1.+0.5*e2*math.cos(f))

    def get_declination(self):
        #This function is unused

        return -.4091*math.cos(2*math.pi*(Radiation.Ls+10)/365)

    def get_zenith(self, t, coord):
        """ Calculates adjusted solar zenith angle at elevation

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
        r"""Air Mass at elevation

        .. math:: AM = 1229+(614cos(\zeta)^2)^{\frac{1}{2}}-614cos(\zeta)

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
        r"""The amount of solar radiation that permeates through the atmosphere at a
        certain altitude, I_{sun} is driven by the atmospheric transmittance.

        .. math:: \tau_{atm}= \frac{1}{2}(e^{-0.65AM}+e^{-0.095AM})

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
        """Calculates total radiation sources as a function of altitude, time, and balloon surface area.

        The figure below shows how different altitudes effects the radiation sources on
        a particular date and coordinate for Tucson Arizona (at sruface level and 25 km altitudes)

        .. image:: ../../img/Tucson_Radiation_Comparison.png

        """

        zen = self.get_zenith(datetime, coord)
        el = coord["alt"]

        #radRef = Radiation.radRef + Radiation.radRef*Radiation.radRef +  Radiation.radRef*Radiation.radRef*Radiation.radRef
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
        power_earth_IR = earth_IR*Radiation.emissEnv*vf*Radiation.surfArea

        sky_IR = self.get_sky_IR(el)
        power_sky_IR = sky_IR*Radiation.emissEnv*(1.-vf)*Radiation.surfArea

        rad_tot_bal = power_direct + power_diffuse + power_reflected + power_earth_IR + power_sky_IR

        return rad_tot_bal
