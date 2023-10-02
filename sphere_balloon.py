import math
import numpy as np
import radiation

import config_earth

"""
sphere_balloon.py solves for the total heat transfer on the solar balloon.
"""

class Sphere_Balloon:
    """Initializes atmospheric properties from the earth configuration file"""
    Cp_air0 = config_earth.earth_properties['Cp_air0']
    Rsp_air = config_earth.earth_properties['Rsp_air']
    cv_air0 = config_earth.earth_properties['Cv_air0']
    cf = config_earth.balloon_properties['cp']

    RE = 6371000.0      # (m) Radius of Earth
    SB = 5.670373E-8    # Stefan Boltzman Constant

    def __init__(self):
        """Initializes all of the solar balloon paramaters from the configuration file"""
        self.d = config_earth.balloon_properties['d']
        self.emissEnv = config_earth.balloon_properties['emissEnv']

        self.surfArea = math.pi*self.d*self.d
        self.vol = math.pi*4/3*pow((self.d/2),3)


    def setEmiss(self,e):
        self.emissEnv = e

    def get_viscocity(self,T):
        """Calculates Kinematic Viscocity of Air at Temperature,T

        :param T: Temperature (K)
        :type Ra: float
        :returns: mu, Kinematic Viscocity of Air
        :rtype: float
        """
        #print("viscocity",T)
        return 1.458E-6*(np.sign(T) * (np.abs(T)) ** (1.5))/(T+110.4) #numpy power does not allow fractional powers of negative numbers. This is the workaround

    def get_conduction(self,T):
        """Calculates Thermal Diffusivity of Air at Temperature, T using Sutherland's Law of Thermal Diffusivity

        :param T: Temperature (K)
        :type Ra: float
        :returns: Thermal Diffusivity of Air (W/(m*K)
        :rtype: float
        """

        #print("conduction",T)



        return 0.0241*(np.sign(T) * (np.abs(T/273.15)) ** (0.9))

    def get_Pr(self,T):
        """Calculates Prantl Number

        :param T: Temperature (K)
        :type Ra: float
        :returns: Prantl Number
        :rtype: float
        """
        k = self.get_conduction(T) #Thermal diffusivity
        Pr = self.get_viscocity(T)*Sphere_Balloon.Cp_air0/k
        return Pr

    #-------------------------------------------SOLVE FOR T_S------------------------------------------------------------------'''

    def get_Nu_ext(self,Ra, Re, Pr):
        """Calculates External Nusselt Number

        :param Ra: Raleigh's number
        :type Ra: float
        :param Re: Reynold's number
        :type Re: float
        :param Pr: Prandtl Number
        :type Pr: float
        :returns: External Nusselt Number
        :rtype: float
        """

        Nu_n = 0.0
        if Ra < 1.5E8:
            Nu_n = 2.0 + 0.6*pow(Ra,0.25)
        else:
            Nu_n = 0.1*pow(Ra, 0.34)
        Nu_f = 0.0
        if Re < 5E4:
            try:
                Nu_f = 2 + 0.47*math.sqrt(Re)*pow(Pr, (1./3.))
            except:
                Nu_f = 2
        else:
            Nu_f = (0.0262*pow(Re, 0.8) - 615.)*pow(Pr, (1./3.));
        return np.fmax(Nu_f, Nu_n);

    def get_q_ext(self, T_s, el, v):
        """Calculate External Heat Transfer to balloon envelope

        :param zen: Surface Temperature of Envelope (K)
        :type zen: float
        :param el: Elevation (m)print fluids.atmosphere.solar_position(datetime.datetime(2018, 4, 15, 6, 43, 5), 51.0486, -114.07)[0]
        :type el: float
        :param el: velocity (m/s)
        :type el: float
        :returns: Power transferred from sphere to surrounding atmosphere due to convection(W)
        :rtype: float
        """

        rad = radiation.Radiation()
        T_atm = rad.getTemp(el)
        p_atm = rad.getPressure(el)
        rho_atm = rad.getDensity(el)
        g = rad.getGravity(el)

        Pr_atm = self.get_Pr(T_atm)

        T_avg = 0.5*(T_atm + T_s)
        rho_avg = p_atm/(Sphere_Balloon.Rsp_air*T_avg)
        Pr_avg = self.get_Pr(T_avg)

        exp_coeff = 1./T_avg;
        kin_visc = self.get_viscocity(T_avg)/rho_avg

        #Not sure if Raleighs number is the right equation here:
        Ra = Pr_avg*g*math.fabs(T_s-T_atm)*np.power(self.d,3)*exp_coeff/(kin_visc*kin_visc)
        Re = rho_atm*v*self.d/self.get_viscocity(T_atm)
        Nu = self.get_Nu_ext(Ra, Re, Pr_atm)
        k = self.get_conduction(T_avg)
        h = Nu*k/self.d
        return h*self.surfArea*(T_s-T_atm)

    def get_sum_q_surf(self,q_rad, T_s, el, v):
        """External Heat Transfer

        :param q_rad: Power input from external radiation (W)
        :type q_rad: float
        :param T_s: Surface Temperature of Envelope (K)
        :type T_s: float
        :param el: Elevation (m)
        :type el: float
        :param v: velocity (m/s)
        :type v: float
        :returns: The sum of power input to the balloon surface (W)
        :rtype: float
        """

        q_conv_loss = -self.get_q_ext(T_s, el, v)
        q_rad_lost = -self.emissEnv*Sphere_Balloon.SB*np.power(T_s,4)*self.surfArea
        return q_rad + q_conv_loss + q_rad_lost

    #--------------------------------------------SOLVE FOR T INT-------------------------------------------------------------

    def get_Nu_int(sef,Ra):
        """Calculates Internal Nusselt Number

        :param Ra: Raleigh's number
        :type Ra: float
        :returns: Internal Nusselt Number
        :rtype: float
        """

        if Ra < 1.35E8:
            return 2.5*(2+0.6*pow(Ra,0.25))
        else:
            return 0.325*pow(Ra, 0.333)

    def get_q_int(self,T_s, T_i, el):
        """Calculates Internal Heat Transfer

        :param T_s: Surface Temperature of Envelope (K)
        :type T_s: float
        :param el: Elevation (m)
        :type el: float
        :param v: velocity (m/s)
        :type v: float
        :returns: Internal Heat Transfer (W)
        :rtype: float
        """

        rad = radiation.Radiation()
        T_atm = rad.getTemp(el)
        p_atm = rad.getPressure(el)
        rho_atm = rad.getDensity(el)
        g = rad.getGravity(el)

        '''
        T_avg = 0.5*(T_s+T_i)
        rho_avg = p_atm/(Sphere_Balloon.Rsp_air*T_avg)
        Pr = self.get_Pr(T_avg)
        exp_coeff = 1./T_avg
        kin_visc = self.get_viscocity(T_avg)/rho_avg
        Ra = self.get_Pr(T_atm)*g*math.fabs(T_i-T_s)*pow(self.d,3)*exp_coeff/(kin_visc*kin_visc)
        Nu = self.get_Nu_int(Ra)
        k = self.get_conduction(T_avg)
        h = (Nu*k)/self.d
        q_int = h*self.surfArea*(T_s-T_i)
        '''

        T_avg = 0.5*(T_s+T_i)
        Pr = self.get_Pr(T_avg)
        rho_avg = p_atm/(Sphere_Balloon.Rsp_air*T_avg)
        mu = self.get_viscocity(T_avg)/rho_avg
        k = self.get_conduction(T_avg)

        inside = (np.power(rho_atm,2)*g*math.fabs(T_s-T_i)*Pr)/(T_i*np.power(mu,2))
        h = 0.13*k * np.sign(inside) * np.abs(inside) ** (1/3.)
        q_int = h*self.surfArea*(T_s-T_i)

        return q_int

    def get_sum_q_int(self, T_s, T_i, el):
        """Calculates sum of Internal Heat Transfer.

        .. note::

            Currently there are no initial heat sources. So this function returns the negative of *get_q_int()*

        :param T_s: Surface Temperature of Envelope (K)
        :type T_s: float
        :param el: Elevation (m)
        :type el: float
        :param v: velocity (m/s)
        :type v: float
        :returns: SUm of Internal Heat Transfer (W)
        :rtype: float
        """
        q_conv_int = self.get_q_int(T_s, T_i, el)
        return q_conv_int
