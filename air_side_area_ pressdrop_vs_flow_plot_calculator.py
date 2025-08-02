
import math
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
try:
    from CoolProp.CoolProp import PropsSI
    coolprop_available = True
except ImportError:
    coolprop_available = False

def air_properties_lookup(T_C):
    T_table = [0, 10, 20, 30, 40, 50, 60]
    mu_table = [1.71e-5, 1.75e-5, 1.81e-5, 1.87e-5, 1.92e-5, 1.98e-5, 2.03e-5]
    rho_table = [1.293, 1.247, 1.204, 1.165, 1.127, 1.093, 1.060]
    T_C = max(0, min(60, T_C))
    for i in range(len(T_table)-1):
        if T_table[i] <= T_C <= T_table[i+1]:
            frac = (T_C - T_table[i]) / (T_table[i+1] - T_table[i])
            mu = mu_table[i] + frac * (mu_table[i+1] - mu_table[i])
            rho = rho_table[i] + frac * (rho_table[i+1] - rho_table[i])
            return rho, mu
    return rho_table[-1], mu_table[-1]

def calculate_pressure_drop(
    tube_od_mm, tube_pitch_mm, row_pitch_mm,
    fin_thickness_mm, fpi, num_rows, face_width_m, face_height_m,
    air_flow_cmh, air_temp_C
):
    tube_od_m = tube_od_mm / 1000
    tube_pitch_m = tube_pitch_mm / 1000
    row_pitch_m = row_pitch_mm / 1000
    fin_thickness_m = fin_thickness_mm / 1000
    fins_per_m = fpi * 39.3701
    fin_spacing_m = 1 / fins_per_m
    frontal_area_m2 = face_width_m * face_height_m
    tubes_per_row = math.floor(face_width_m / tube_pitch_m)

    open_area_per_gap = face_width_m * (fin_spacing_m - fin_thickness_m)
    total_open_area = open_area_per_gap * fins_per_m * face_height_m
    frontal_tube_blockage = tubes_per_row * (math.pi / 4) * tube_od_m**2 * face_height_m
    net_free_flow_area = total_open_area - frontal_tube_blockage

    air_flow_m3s = air_flow_cmh / 3600
    air_velocity_ms = air_flow_m3s / net_free_flow_area if net_free_flow_area > 0 else 0

    if coolprop_available:
        T_K = air_temp_C + 273.15
        rho = PropsSI('D', 'T', T_K, 'P', 101325, 'Air')
        mu = PropsSI('V', 'T', T_K, 'P', 101325, 'Air')
        if mu <= 0 or math.isnan(mu):
            rho, mu = air_properties_lookup(air_temp_C)
    else:
        rho, mu = air_properties_lookup(air_temp_C)

    m_dot = rho * air_flow_m3s
    G = m_dot / net_free_flow_area if net_free_flow_area > 0 else 0
    perimeter_per_gap = 2 * face_width_m + tubes_per_row * (math.pi * tube_od_m)
    D_h = (4 * net_free_flow_area) / perimeter_per_gap if perimeter_per_gap > 0 else 0
    Re = G * D_h / mu if mu > 0 else 0
    f = 0.35 * Re ** -0.2 if Re > 0 else 0
    flow_depth = num_rows * row_pitch_m
    dP = f * (flow_depth / D_h) * (G**2) / (2 * rho) if D_h > 0 else 0

    return air_velocity_ms, Re, f, D_h, dP

st.title("Air-Side Pressure Drop vs Air Flow")

tube_od_mm = st.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_pitch_mm = st.number_input("Tube Pitch (mm)", value=25.4)
row_pitch_mm = st.number_input("Row Pitch (mm)", value=25.4)
fin_thickness_mm = st.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.number_input("Fins Per Inch (FPI)", value=12)
num_rows = st.number_input("Number of Rows", value=4)
face_width_m = st.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.number_input("Coil Face Height (m)", value=1.0)
air_temp_C = st.number_input("Air Temperature (°C)", value=35.0)

air_flow_range = st.slider("Air Flow Rate (m³/h)", 5000, 25000, (5000, 20000), step=1000)

flow_rates = []
pressure_drops = []

for flow in range(air_flow_range[0], air_flow_range[1] + 1, 1000):
    _, _, _, _, dp = calculate_pressure_drop(
        tube_od_mm, tube_pitch_mm, row_pitch_mm,
        fin_thickness_mm, fpi, num_rows,
        face_width_m, face_height_m,
        flow, air_temp_C
    )
    flow_rates.append(flow)
    pressure_drops.append(dp)

fig, ax = plt.subplots()
ax.plot(flow_rates, pressure_drops, marker='o')
ax.set_xlabel("Air Flow Rate (m³/h)")
ax.set_ylabel("Pressure Drop (Pa)")
ax.set_title("Pressure Drop vs Air Flow")
ax.grid(True)

st.pyplot(fig)
