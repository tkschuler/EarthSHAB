import math
import radiation
import sphere_balloon
import config_earth  #Import parameters from configuration file.

""" solve_states.py uses numerical integration to solve for the dynamic response of the balloon.
"""

class SolveStates:
    def __init__(self):
        """Initializes all of the solar balloon paramaters from the configuration file"""

        self.Cp_air0 = config_earth.earth_properties['Cp_air0']
        self.Rsp_air = config_earth.earth_properties['Rsp_air']

        self.d = config_earth.balloon_properties['d']
        self.vol = math.pi*4/3*pow((self.d/2),3) #volume m^3
        self.surfArea = math.pi*self.d*self.d #m^2
        self.cs_area = math.pi*self.d*self.d/4.0 #m^2

        #self.emissEnv = config_earth.balloon_properties['emissEnv']
        self.areaDensityEnv = config_earth.balloon_properties['areaDensityEnv']
        self.mp = config_earth.balloon_properties['mp']
        self.mdot = 0
        self.massEnv = config_earth.balloon_properties['mEnv']
        self.Upsilon = config_earth.balloon_properties['Upsilon']

        self.vent = config_earth.simulation['vent']
        self.coord = config_earth.simulation['start_coord']
        self.t = config_earth.simulation['start_time']
        self.lat = math.radians(self.coord['lat'])
        self.Ls = self.t.timetuple().tm_yday
        self.min_alt = config_earth.simulation['min_alt']

        self.vm_coeff = .1 #virtual mass coefficient
        self.k = self.massEnv*config_earth.balloon_properties['cp'] #thermal mass coefficient

        self.dt = config_earth.dt

    def get_acceleration(self,v,el,T_s,T_i):
        """Solves for the acceleration of the solar balloon after one timestep (dt).

        :param T_s: Surface Temperature (K)
        :type T_s: float
        :param T_i: Internal Temperature (K)
        :type T_i: float
        :param el: Elevation (m)
        :type el: float
        :param v: Velocity (m)
        :type v: float

        :returns: acceleration of balloon (m/s^2)
        :rtype: float
        """

        rad = radiation.Radiation()
        T_atm = rad.getTemp(el)
        p_atm = rad.getPressure(el)
        rho_atm = rad.getDensity(el)
        g = rad.getGravity(el)


        rho_int = p_atm/(self.Rsp_air*T_i) # Internal air density

        Cd = .5 # Drag Coefficient
        F_b = (rho_atm - rho_int)*self.vol*g # Force due to buyoancy
        F_d =  Cd*(0.5*rho_atm*math.fabs(v)*v)*self.cs_area# Force due to Drag

        if F_d > 0:
            F_d = F_d * self.Upsilon
        vm = (self.massEnv + self.mp) + rho_atm*self.vol + self.vm_coeff*rho_atm*self.vol #Virtual Mass
        accel = ((F_b  - F_d - (self.massEnv + self.mp)*g)/vm)

        return accel

    def get_convection_vent(self,T_i,el):
        """Calculates the heat lost to the atmosphere due to venting

        :param T_i: Internal Temperature (K)
        :type T_i: float
        :param el: Elevation (m)
        :type el: float

        :returns: Convection due to Venting (unit?)
        :rtype: float
        """

        rad = radiation.Radiation()
        T_atm = rad.getTemp(el)

        Q_vent =  self.mdot*self.Cp_air0*(T_i-T_atm) # Convection due to released air
        return Q_vent


    def solveVerticalTrajectory(self,t,T_s,T_i,el,v,coord,alt_sp,v_sp):
        """This function numerically integrates and solves for the change in Surface Temperature, Internal Temperature, and accelleration
        after a timestep, dt.

        :param t: Datetime
        :type t: datetime
        :param T_s: Surface Temperature (K)
        :type T_s: float
        :param T_i: Internal Temperature (K)
        :type T_i: float
        :param el: Elevation (m)
        :type el: float
        :param v: Velocity (m)
        :type v: float
        :param alt_sp: Altitude Setpoint (m)
        :type alt_sp: float
        :param v_sp: Velocity Setpoint (m/s)
        :type v_sp: float

        :returns: Updated parameters after dt (seconds)
        :rtype: float [T_s,T_i,el,v]
        """

        bal = sphere_balloon.Sphere_Balloon()
        rad = radiation.Radiation()

        T_atm = rad.getTemp(el)
        p_atm = rad.getPressure(el)
        rho_atm = rad.getDensity(el)

        rho_int = p_atm/(self.Rsp_air*T_i)
        tm_air = rho_int*self.vol*self.Cp_air0

        #Numerically integrate change in Surface Temperature
        coord["alt"] = el
        q_rad  = rad.get_rad_total(t,coord)
        q_surf = bal.get_sum_q_surf(q_rad, T_s, el, v)
        q_int  = bal.get_sum_q_int(T_s, T_i, el)
        dT_sdt = (q_surf-q_int)/self.k

        #Numerically integrate change in Surface Temperature
        tm_air = rho_atm*self.vol*self.Cp_air0
        dT_idt = (q_int-self.get_convection_vent(T_i,el))/tm_air

        #Add the new surface and internal Temperatures
        T_s_new = T_s+dT_sdt*self.dt
        T_i_new = T_i+dT_idt*self.dt

        #solve for accellration, position, and velocity
        dzdotdt = self.get_acceleration(v,el,T_s,T_i)
        zdot = v + dzdotdt*self.dt
        z = el+zdot*self.dt

        #Add the new velocity and position
        if z < self.min_alt:
            v_new = 0
            el_new = self.min_alt
        else:
            v_new = zdot
            el_new = z

        # Venting commands for an altitude setpoint. Vent is either on or off.
        if el_new > alt_sp:
            self.mdot = self.vent

        if el_new < alt_sp:
            self.mdot = 0

        return [T_s_new,T_i_new,T_atm,el_new,v_new, q_rad, q_surf, q_int]
