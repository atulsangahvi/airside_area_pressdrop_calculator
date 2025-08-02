
import math
import streamlit as st
import pandas as pd
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

def calculate_air_side_results(
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
    fin_depth_m = num_rows * row_pitch_m
    tubes_per_row = math.floor(face_width_m / tube_pitch_m)
    total_tubes = tubes_per_row * num_rows
    tube_ext_area = total_tubes * (math.pi * tube_od_m)

    fin_area_per_fin = 2 * face_width_m * fin_depth_m
    total_gross_fin_area = fin_area_per_fin * fins_per_m
    hole_area_per_tube = (math.pi / 4) * tube_od_m**2
    total_hole_area = hole_area_per_tube * total_tubes * fins_per_m
    net_fin_area = total_gross_fin_area - total_hole_area
    total_air_side_area = (tube_ext_area + net_fin_area) * face_height_m

    open_area_per_gap = face_width_m * (fin_spacing_m - fin_thickness_m)
    total_open_area = open_area_per_gap * fins_per_m * face_height_m
    frontal_tube_blockage = tubes_per_row * (math.pi / 4) * tube_od_m**2 * face_height_m
    net_free_flow_area = total_open_area - frontal_tube_blockage
    percent_free_area = 100 * net_free_flow_area / frontal_area_m2
    air_flow_m3s = air_flow_cmh / 3600
    face_velocity = air_flow_m3s / frontal_area_m2
    u_max = air_flow_m3s / net_free_flow_area if net_free_flow_area > 0 else 0

    if coolprop_available:
        T_K = air_temp_C + 273.15
        rho = PropsSI('D', 'T', T_K, 'P', 101325, 'Air')
        mu = PropsSI('V', 'T', T_K, 'P', 101325, 'Air')
        if mu <= 0 or math.isnan(mu):
            rho, mu = air_properties_lookup(air_temp_C)
    else:
        rho, mu = air_properties_lookup(air_temp_C)

    m_dot = rho * air_flow_m3s

    passage_height = tube_pitch_m - tube_od_m
    passage_width = fin_spacing_m - fin_thickness_m
    A_min_cell = passage_height * passage_width
    P_wet_cell = 2 * (passage_height + passage_width)
    D_h = (4 * A_min_cell) / P_wet_cell if P_wet_cell > 0 else 0

    Re = (rho * u_max * D_h) / mu if mu > 0 else 0
    f = 0.25 * Re**-0.25 if Re > 0 else 0
    flow_depth = num_rows * row_pitch_m
    dP = f * (flow_depth / D_h) * (m_dot / net_free_flow_area)**2 / (2 * rho) if D_h > 0 else 0

    return {
        "Tubes per row": tubes_per_row,
        "Total tubes": total_tubes,
        "Fin depth (m)": fin_depth_m,
        "Tube external area (m²)": tube_ext_area,
        "Net fin area (m²)": net_fin_area,
        "Total air side area (m²)": total_air_side_area,
        "Free flow area (m²)": net_free_flow_area,
        "Free flow area (%)": percent_free_area,
        "Face velocity (m/s)": face_velocity,
        "Max fin passage velocity (m/s)": u_max,
        "Hydraulic diameter (m)": D_h,
        "Air density (kg/m³)": rho,
        "Air viscosity (Pa·s)": mu,
        "Mass flow rate (kg/s)": m_dot,
        "Reynolds number (corrected)": Re,
        "Friction factor": f,
        "Air-side Pressure Drop (Pa)": dP
    }

st.title("Air-Side Area and Pressure Drop Calculator (Corrected Reynolds Number)")

tube_od_mm = st.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_pitch_mm = st.number_input("Tube Pitch (mm)", value=25.4)
row_pitch_mm = st.number_input("Row Pitch (mm)", value=25.4)
fin_thickness_mm = st.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.number_input("Fins Per Inch (FPI)", value=12)
num_rows = st.number_input("Number of Rows", value=4)
face_width_m = st.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.number_input("Coil Face Height (m)", value=1.0)
air_flow_cmh = st.number_input("Air Flow Rate (m³/h)", value=10000)
air_temp_C = st.number_input("Air Temperature (°C)", value=35.0)

if st.button("Calculate"):
    results = calculate_air_side_results(
        tube_od_mm, tube_pitch_mm, row_pitch_mm,
        fin_thickness_mm, fpi, num_rows,
        face_width_m, face_height_m,
        air_flow_cmh, air_temp_C
    )
    df = pd.DataFrame(list(results.items()), columns=["Parameter", "Value"])
    df["Value"] = df["Value"].apply(lambda x: f"{x:.6f}" if isinstance(x, float) else x)
    st.table(df)
