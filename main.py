# The Aeroinformatics Engineering 

#calculating lift 
'''
def lift(v, s, rho, cl):
    return 0.5 * v**2 * s * rho * cl

#calculating drag

def drag(v, s, rho, cd):
    return 0.5 * v**2 * s * rho * cd


#calculating lift to drag ratio 

v =  float(input("Enter the velocity of the aircraft in m/s: "))
s = float(input("Enter the wing area in m^2: "))

rho = float(input("Enter the air density in kg/m^3: "))

cl = float(input("Enter the lift coefficient: "))

cd = float(input("Enter the drag coefficient: "))

plane_lift = lift(v,s, rho, cl)
plane_drag = drag(v, s, rho, cd)

lift_to_drag_ratio = plane_lift/ plane_drag

print("The lift of the aircraft is: ", plane_lift, "N")
print("The drag of the aircraft is: ", plane_drag, "N")

print("The lift to drag rati
'''

import math
import csv

# ---------- Atmospheric Model ----------
def air_density(altitude):
    T0, P0, L, R, g = 288.15, 101325, 0.0065, 287.05, 9.80665
    T = T0 - L * altitude
    P = P0 * (T / T0) ** (g / (R * L))
    rho = P / (R * T)
    return rho

# ---------- Aerodynamic Model ----------
def coefficients(alpha_deg, CL0=0.7, a=0.9, CD0=0.044, k=0.024):
    CL = CL0 + a * alpha_deg
    CD = CD0 + k * (CL ** 2)
    return CL, CD

# ---------- Force Calculations ----------
def aerodynamic_forces(V, rho, S, CL, CD):
    L = 0.5 * rho * (V ** 2) * S * CL
    D = 0.5 * rho * (V ** 2) * S * CD
    return L, D

# ---------- Aeron Factor ----------
def aeron_factor(CL, CD):
    return CL / CD

# ---------- Autonomous Simulation ----------
def autonomous_flight_optimizer(altitude, V, S, alpha_range):
    rho = air_density(altitude)
    best_condition = {"AeronFactor": 0}
    data_log = []

    for alpha in alpha_range:
        CL, CD = coefficients(alpha)
        L, D = aerodynamic_forces(V, rho, S, CL, CD)
        AF = aeron_factor(CL, CD)
        data_log.append([alpha, CL, CD, L, D, AF])

        if AF > best_condition["AeronFactor"]:
            best_condition = {
                "Alpha": alpha,
                "CL": CL,
                "CD": CD,
                "Lift": L,
                "Drag": D,
                "AeronFactor": AF,
                "Rho": rho,
                "V": V,
                "S": S
            }

    # Log data to CSV
    with open("aeron_data_log.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Alpha", "CL", "CD", "Lift(N)", "Drag(N)", "AeronFactor"])
        writer.writerows(data_log)

    # Print the optimal result
    print("\n--- Optimal Flight Conditions ---")
    print(f"Altitude     : {altitude} m")
    print(f"Airspeed     : {V} m/s")
    print(f"Best AoA     : {best_condition['Alpha']}°")
    print(f"CL           : {best_condition['CL']:.3f}")
    print(f"CD           : {best_condition['CD']:.3f}")
    print(f"Lift         : {best_condition['Lift']:.2f} N")
    print(f"Drag         : {best_condition['Drag']:.2f} N")
    print(f"Aeron Factor : {best_condition['AeronFactor']:.2f}")
    print(f"Air Density  : {best_condition['Rho']:.3f} kg/m³")

    return best_condition

# ---------- Run Example ----------
if __name__ == "__main__":
    altitude = 3000   # m
    V = 75            # m/s
    S = 16            # m²
    alpha_range = range(0, 15)  # degrees
    autonomous_flight_optimizer(altitude, V, S, alpha_range)