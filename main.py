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


import math
import matplotlib.pyplot as plt
import numpy as np

# ---------- Atmospheric Model ----------
def air_density(altitude):
    T0, P0, L, R, g = 288.15, 101325, 0.0065, 287.05, 9.80665
    T = T0 - L * altitude
    P = P0 * (T / T0) ** (g / (R * L))
    rho = P / (R * T)
    return rho

# ---------- Aerodynamic Model ----------
def coefficients(alpha_deg, CL0=0.3, a=0.1, CD0=0.02, k=0.05):
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

# ---------- Dynamic Simulation & Visualization ----------
def dynamic_simulation(V=75, S=16, altitudes=[0, 2000, 5000]):
    alpha_range = np.linspace(0, 14, 30)  # AoA range (0°–14°)

    plt.figure(figsize=(10, 8))

    for altitude in altitudes:
        rho = air_density(altitude)
        AF_list, Lift_list, Drag_list = [], [], []

        for alpha in alpha_range:
            CL, CD = coefficients(alpha)
            L, D = aerodynamic_forces(V, rho, S, CL, CD)
            AF = aeron_factor(CL, CD)
            AF_list.append(AF)
            Lift_list.append(L)
            Drag_list.append(D)

        # Plot Aeron Factor Curve
        plt.plot(alpha_range, AF_list, label=f'Altitude: {altitude} m')

    plt.title("Aeron Factor vs Angle of Attack")
    plt.xlabel("Angle of Attack (°)")
    plt.ylabel("Aeron Factor (Cl/Cd)")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Optional: Secondary plots for Lift and Drag
    plt.figure(figsize=(10, 8))
    for altitude in altitudes:
        rho = air_density(altitude)
        L_list, D_list = [], []
        for alpha in alpha_range:
            CL, CD = coefficients(alpha)
            L, D = aerodynamic_forces(V, rho, S, CL, CD)
            L_list.append(L)
            D_list.append(D)
        plt.plot(alpha_range, L_list, label=f'Lift @ {altitude} m')
        plt.plot(alpha_range, D_list, linestyle='--', label=f'Drag @ {altitude} m')

    plt.title("Lift & Drag vs Angle of Attack")
    plt.xlabel("Angle of Attack (°)")
    plt.ylabel("Force (N)")
    plt.legend()
    plt.grid(True)
    plt.show()

# ---------- Run Dynamic Visualization ----------
if __name__ == "__main__":
    dynamic_simulation(
        V=75, 
        S=16, 
        altitudes=[0, 2000, 5000]
    )
    '''
import math
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ---------- Atmospheric Model ----------
def air_density(altitude):
    T0, P0, L, R, g = 288.15, 101325, 0.0065, 287.05, 9.80665
    T = T0 - L * altitude
    P = P0 * (T / T0) ** (g / (R * L))
    rho = P / (R * T)
    return rho

# ---------- Aerodynamic Model ----------
def coefficients(alpha_deg, CL0=0.3, a=0.1, CD0=0.02, k=0.05):
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
    return CL / CD if CD != 0 else 0

# ---------- Decision Engine ----------
def aeronix_decision(L, W, AF, prev_AF):
    if L < W:
        return "↑ Increase AoA (Low Lift)"
    elif AF < prev_AF * 0.9:
        return "↘ Reduce AoA (Efficiency Drop)"
    elif AF > prev_AF * 1.05:
        return "↗ Stable - Maintain"
    else:
        return "✓ Optimal Cruise"

# ---------- Real-Time Simulation ----------
def run_aeronix_simulation():
    # Aircraft Parameters
    S = 16          # Wing area (m²)
    W = 6000        # Weight (N)
    V = 75          # Initial airspeed (m/s)
    altitude = 0    # Initial altitude (m)
    alpha = 5       # Initial angle of attack (°)
    
    # Data lists
    time_data, lift_data, drag_data, af_data, alpha_data, decisions = [], [], [], [], [], []

    # Visualization setup
    plt.style.use("dark_background")
    fig, ax = plt.subplots(3, 1, figsize=(10, 8))
    plt.subplots_adjust(hspace=0.5)
    
    ax[0].set_title("Agent AeroniX - Lift & Drag vs Time")
    ax[1].set_title("Aeron Factor Evolution")
    ax[2].set_title("Angle of Attack Adjustment")

    line1, = ax[0].plot([], [], 'cyan', label='Lift (N)')
    line2, = ax[0].plot([], [], 'magenta', label='Drag (N)')
    line3, = ax[1].plot([], [], 'lime', label='Aeron Factor')
    line4, = ax[2].plot([], [], 'orange', label='AoA (°)')

    for a in ax:
        a.legend()
        a.grid(True, alpha=0.3)
    
    prev_AF = 0
    t = 0

    def update(frame):
        nonlocal altitude, V, alpha, prev_AF, t

        # Dynamic environment changes
        altitude += random.uniform(-10, 50)  # m change
        V += random.uniform(-2, 2)           # airspeed fluctuation
        rho = air_density(max(0, altitude))

        # Aerodynamic calculations
        CL, CD = coefficients(alpha)
        L, D = aerodynamic_forces(V, rho, S, CL, CD)
        AF = aeron_factor(CL, CD)
        decision = aeronix_decision(L, W, AF, prev_AF)
        
        # Store data
        time_data.append(t)
        lift_data.append(L)
        drag_data.append(D)
        af_data.append(AF)
        alpha_data.append(alpha)
        decisions.append(decision)
        t += 1
        
        # Decision updates
        if "Increase" in decision:
            alpha += 0.2
        elif "Reduce" in decision:
            alpha -= 0.2
        elif "Maintain" in decision:
            alpha += random.uniform(-0.05, 0.05)
        
        prev_AF = AF
        
        # Update plots
        line1.set_data(time_data, lift_data)
        line2.set_data(time_data, drag_data)
        line3.set_data(time_data, af_data)
        line4.set_data(time_data, alpha_data)
        
        for a in ax:
            a.relim()
            a.autoscale_view()

        # Live console output
        print(f"[t={t}s] Alt={altitude:.1f}m | V={V:.1f}m/s | α={alpha:.2f}° | AF={AF:.2f} | Decision: {decision}")

        return line1, line2, line3, line4

    ani = FuncAnimation(fig, update, interval=500)
    plt.show()

# ---------- Run Agent AeroniX ----------
if __name__ == "__main__":
    run_aeronix_simulation()