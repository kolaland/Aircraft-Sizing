
def sizing_aircraft(AR):
    # IMPORT
    import math
    import numpy as np
    from scipy import interpolate

    # [-] Storing pi for ease
    pi = math.pi


    # INPUTS------------------------------------------------------------------------------------------------------------

    # [kg] Payload weight
    W_pay = 2
    # [-] Drag Coefficient
    CD0 = 0.025
    # [m] cruise altitude ASL
    alt = 1000
    # [kg/m^3] density of air at 1000 m [from engineeringtoolbox.com]
    rho = 1.225 * 0.9075
    # [deg] leading edge sweep:
    sweep_LE = 15
    # [deg] change in stall angle of attack
    del_alpha_stall = 1.179
    # [-] efficiency ratio
    mu_p = 0.70
    # [m^2] surface area
    S_area = 1.8
    # [W] Max Power
    P_max = 150
    # [-] Reynolds Number
    Re = 200000


    # WEIGHT ESTIMATION-------------------------------------------------------------------------------------------------

    # [-] Weight equation 1
    def eq1(W_TO):
        return 0.916 * (W_TO ** (-0.0795))
    # [-] Weight equation 2
    def eq2(W_TO):
        return 1 - ((0.11 * AR * W_pay) / (W_TO))
    # [-] Difference between eq1 and eq2
    def equations(W_TO):
        return eq1(W_TO) - eq2(W_TO)

    # [-] Using Bisection method
    def bisection_method(func, a, b, tol=1e-6, max_iter=1000):
        # [-] Ensuring that func(a) and func(b) have opposite signs
        if func(a) * func(b) > 0:
            return None
        # [-] Initial iteration
        iter_count = 0
        # [-] Bisection condition
        while (b - a) / 2.0 > tol and iter_count < max_iter:
            # [-] Midpoint
            c = (a + b) / 2.0
            # [-] For exact root
            if func(c) == 0:
                return c
            # [-] Root lies between a and c
            elif func(a) * func(c) < 0:
                b = c
            # [-] Root lies between c and b
            else:
                a = c
            # [-] Update iteration
            iter_count += 1

        # [-] Returns the midpoint as the best approximation of the root
        return (a + b) / 2.0

    # [-] Lower bound of the interval guess
    W_TO_lower = 1
    # [-] Upper bound of the interval guess
    W_TO_upper = 100
    # [-] Calling the bisection function to find the root
    W_TO_new = bisection_method(equations, W_TO_lower, W_TO_upper)
    if W_TO_new is not None:
        print(f"Calculated W_TO using Bisection method: {W_TO_new:.4g}")
    else:
        print("No solution found within the given interval.")

    # AEROD SIZING------------------------------------------------------------------------------------------------------

    # [-] Statement takes average of two e equations
    if 4 < AR < 10:
        ows_e1 = (4.61 * (1 - (0.045 * (AR**0.68))) * (math.cos(math.radians(sweep_LE))**0.15)) - 3.1
        ows_e2 = 2 / (2 - AR + math.sqrt(4 + (AR**2 * (1 + (math.tan(math.radians(sweep_LE))**2) ))))
        ows_e = (ows_e1 + ows_e2) / 2

    # [-] Statement for Equation 2
    elif 4 < AR < 15:
        ows_e = 2 / (2 - AR + math.sqrt(4 + (AR**2 * (1 + (math.tan(math.radians(sweep_LE))**2)))))

    # [-] Statement for Equation 1
    elif AR < 10:
        ows_e = (4.61 * (1 - (0.045 * (AR**0.68))) * (math.cos(math.radians(sweep_LE))**0.15)) - 3.1

    # [-] Error statement
    else:
        print('Enter an AR less than 15')

    # [-] Induced drag correction factor
    k = 1 / (pi * ows_e * AR)
    # [-] Maximum aerodynamic efficiency
    LD_max = 1 / (math.sqrt(4 * k * CD0))

    # [-] Load airfoil_data from same folder
    File_Data = np.loadtxt("airfoil_data.txt")
    # [deg] Angle of attack
    aoa = File_Data[:, 0]
    # [-] Lift coefficient
    CL = File_Data[:, 1]
    # [-] Interpolating aoa and CL using scipy.interpolate.interp1d
    CL_interp = interpolate.interp1d(aoa, CL, kind='linear', fill_value='extrapolate')
    # [deg] To find zero-lift aoa
    alpha_range = np.linspace(min(aoa), max(aoa), 1000)
    # [-] Gets CL range
    CL_range = CL_interp(alpha_range)
    # [-] Interpolates for alpha_zl
    alpha_zl = np.interp(0, CL_range, alpha_range)
    # [-] Extracts Clmax from maximum of CL data
    Clmax = np.max(CL)
    # [deg] Extracts corresponding angle of attack
    alpha_Clmax = aoa[np.argmax(CL)]

    # [-] Max lift coefficient - clean
    CLmax0 = 0.9 * Clmax
    # [-] 2D max lift coefficient
    CLmax = CLmax0 * math.cos(math.radians(sweep_LE))

    # [-] Finds the slope of aoa and CL data
    coefficients = np.polyfit(aoa, CL, 1)
    # [1/deg] Extracts the slope
    Cl_alpha = coefficients[0]
    # [-] Extracts the y-intercept
    Cl_0 = coefficients[1]

    # [1/deg] 2D lift curve slope
    CL_alpha = (Cl_alpha * AR) / (2 + math.sqrt((AR**2) + 4))
    # [deg] stall angle of attack
    alpha_stall = (CLmax / CL_alpha) + alpha_zl + del_alpha_stall

    W_TO_new = W_TO_new * 9.81
    W_TO = W_TO_new / 9.81

    # PERFORMANCE ANALYSIS----------------------------------------------------------------------------------------------

    # [m/s] Velocity for steepest climb angle
    V_gamma_max = (4 * (W_TO_new**2) * k) / (P_max * rho * mu_p * S_area)

    # [-] redundant, please ignore
    term1 = (P_max * mu_p) / (W_TO_new * V_gamma_max)
    term2 = (rho * (V_gamma_max**2) * S_area * CD0) / (2 * W_TO_new)
    term3 = (2 * k * W_TO_new) / (rho * (V_gamma_max**2) * S_area)
    result = term1 - term2 - term3
    # [deg] steepest climb angle, should be within bounds of -1 to 1 for asin
    # result = max(-1, min(1, result))

    gamma_max = math.asin(result)
    gamma_max = gamma_max * (180/pi)

    # [kg, -. deg. deg] OUTPUTS: Weight at takeoff W_TO, maximum aerodynamic efficiency, stall angle of attack, steepest climb angle
    return W_TO, LD_max, alpha_stall, gamma_max

# Calling sizing function with AR input
[W_TO, LDmax, alpha_stall, gamma_max] = sizing_aircraft(8)

# Prints outputs needed
print(f"W_TO [kg]: {W_TO:.4f}, LDmax [-]: {LDmax:.4f}, alpha_stall [deg]: {alpha_stall:.4f}, gamma_max [deg]: {gamma_max:.4f}")
