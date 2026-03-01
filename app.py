import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# ======================================
# Page Config
# ======================================
st.set_page_config(page_title="Power System Professional Tool", layout="wide")

# ======================================
# Dark Theme
# ======================================
dark_mode = st.sidebar.checkbox("Dark Mode")
if dark_mode:
    plt.style.use("dark_background")
else:
    plt.style.use("default")

st.title("⚡ Power Transmission Analysis")

# ======================================
# Helper Function (Polar Display)
# ======================================
def polar_str(value, unit=""):
    mag = np.abs(value)
    ang = np.degrees(np.angle(value))
    return f"{mag:.4f} ∠ {ang:.2f}° {unit}"

# ======================================
# 1) Operating Condition
# ======================================
st.sidebar.header("Operating Condition")
st.sidebar.header("Stability Study")

# Reactive Limits
st.sidebar.header("Reactive Limits")
Qmax = st.sidebar.number_input("Qmax (MVAR)", value=300.0)
Qmin = st.sidebar.number_input("Qmin (MVAR)", value=-300.0)

run_anim = st.sidebar.checkbox("Run Stability Animation")

mode = st.sidebar.selectbox("Mode", ["Load", "No Load"])

if mode == "Load":
    Pr = st.sidebar.number_input("Receiving Power Pr (MW)", value=100.0)
    PF = st.sidebar.slider("Power Factor", 0.0, 1.0, 0.9)
    pf_type = st.sidebar.selectbox("PF Type", ["Lagging", "Leading"])

    phi = np.arccos(PF) if PF != 0 else 0
    if pf_type == "Leading":
        phi = -phi

    Pr_phase = Pr / 3
    Qr = Pr * np.tan(phi) if PF != 0 else 0
    Qr_phase = Qr / 3
    Sr = Pr_phase + 1j * Qr_phase
    load_mode = "Loaded Operation"

else:
    Pr = 0.0
    Qr = 0.0
    Sr = 0 + 0j
    load_mode = "No-Load (Line charging only)"

st.write("### Receiving End")
st.write("Mode:", load_mode)
st.write("Sr =", polar_str(Sr, "MVA"))

# ======================================
# 2) Line Parameters
# ======================================
st.sidebar.header("Line Parameters")
method = st.sidebar.selectbox(
    "ABCD Input Method",
    ["Enter A,B,C directly", "Enter Z,Y,length", "Enter r,L,C,length"]
)

if method == "Enter A,B,C directly":
    A_mag = st.sidebar.number_input("|A|", value=0.986)
    A_ang = st.sidebar.number_input("Angle A", value=0.1)
    B_mag = st.sidebar.number_input("|B|", value=60.0)
    B_ang = st.sidebar.number_input("Angle B", value=83.0)
    C_mag = st.sidebar.number_input("|C|", value=0.0004)
    C_ang = st.sidebar.number_input("Angle C", value=90.0)

    A = A_mag*np.exp(1j*np.radians(A_ang))
    B = B_mag*np.exp(1j*np.radians(B_ang))
    C = C_mag*np.exp(1j*np.radians(C_ang))
    D = A

elif method == "Enter Z,Y,length":
    R = st.sidebar.number_input("R (ohm/km)", value=0.1)
    X = st.sidebar.number_input("X (ohm/km)", value=0.4)
    G = st.sidebar.number_input("G (S/km)", value=0.0)
    B_sh = st.sidebar.number_input(
    "B (S/km)",
    min_value=0.0,
    value=0.000001,
    format="%.6f"
)
    length = st.sidebar.number_input("Length (km)", value=100.0)

    Z = (R + 1j*X)*length
    Y = (G + 1j*B_sh)*length

    theta = np.sqrt(Z*Y)
    A = np.cosh(theta)
    B = Z*np.sinh(theta)/theta
    C = Y*np.sinh(theta)/theta
    D = A

else:
    r = st.sidebar.number_input("r (ohm/km)", value=0.1)
    L = st.sidebar.number_input("L (H/km)", value=1e-3)
    C_val = st.sidebar.number_input("C (F/km)", value=1e-8)
    length = st.sidebar.number_input("Length (km)", value=100.0)

    w = 2*np.pi*50
    Z = (r + 1j*w*L)*length
    Y = (1j*w*C_val)*length

    theta = np.sqrt(Z*Y)
    A = np.cosh(theta)
    B = Z*np.sinh(theta)/theta
    C = Y*np.sinh(theta)/theta
    D = A

st.write("### ABCD Parameters")
st.write("A =", polar_str(A))
st.write("B =", polar_str(B, "ohm"))
st.write("C =", polar_str(C, "S"))
st.write("D =", polar_str(D))

# ======================================
# 3) Voltage & Current
# ======================================
Vr_mag = st.sidebar.number_input("Receiving Voltage Vr (kV)", value=220.0)

st.sidebar.header("Fault Simulation")
fault = st.sidebar.selectbox("Fault Type", ["No Fault", "3-Phase Fault at Receiving End"])
fault_impedance = st.sidebar.number_input("Fault Impedance (ohm)", value=0.01)

Vr = Vr_mag + 0j
Ir = np.conj(Sr) / Vr if Vr != 0 else 0

if fault == "3-Phase Fault at Receiving End":
    Vr = fault_impedance * Ir

Vs = A*Vr + B*Ir
Is = C*Vr + D*Ir

st.write("### Voltages & Currents")
st.write("Vr =", polar_str(Vr, "kV"))
st.write("Ir =", polar_str(Ir, "kA"))
st.write("Vs =", polar_str(Vs, "kV"))
st.write("Is =", polar_str(Is, "kA"))

# ======================================
# 4) Performance Analysis
# ======================================
Ss = 3*Vs*np.conj(Is)
Ps = np.real(Ss)
Qs = np.imag(Ss)

if Qs > Qmax:
    Q_status = "Overexcited (Above Qmax)"
elif Qs < Qmin:
    Q_status = "Underexcited (Below Qmin)"
else:
    Q_status = "Within Limits"

st.write("Reactive Status:", Q_status)

Loss = Ps - Pr
Efficiency = (Pr/Ps)*100 if Ps != 0 else 0
Voltage_reg = ((np.abs(Vs)-Vr_mag)/Vr_mag)*100

Pmax_calc = (np.abs(Vs)*Vr_mag)/np.abs(B) if B != 0 else 0
Loading = (Pr/Pmax_calc)*100 if Pmax_calc !=0 else 0
Margin = Pmax_calc - Pr

stability = "Stable" if Pr < Pmax_calc else "Unstable"
if Loading < 50:
    condition = "Light Load"
elif Loading < 80:
    condition = "Normal Load"
else:
    condition = "Heavy Load"

st.write("### Performance")
st.write("Ps =", f"{Ps:.2f} MW")
st.write("Qs =", f"{Qs:.2f} MVAR")
st.write("Loss =", f"{Loss:.2f} MW")
st.write("Efficiency =", f"{Efficiency:.2f} %")
st.write("Voltage Regulation =", f"{Voltage_reg:.2f} %")
st.write("Pmax =", f"{Pmax_calc:.2f} MW")
st.write("Loading =", f"{Loading:.2f} %")
st.write("Stability Margin =", f"{Margin:.2f} MW")

# ======================================
# PV Curve
# ======================================
st.subheader("PV Curve (Voltage Stability)")
P_range = np.linspace(0, Pmax_calc*1.2, 50)
V_list = []
for P_test in P_range:
    Q_test = P_test * np.tan(phi) if PF != 0 else 0
    Sr_test = (P_test + 1j*Q_test)/3
    Ir_test = np.conj(Sr_test)/Vr if Vr !=0 else 0
    Vs_test = A*Vr + B*Ir_test
    V_list.append(np.abs(Vs_test))

fig_pv, ax_pv = plt.subplots()
ax_pv.plot(P_range, V_list)
ax_pv.axvline(Pmax_calc, linestyle='--', color='red', label="Pmax")
ax_pv.set_xlabel("Load P (MW)")
ax_pv.set_ylabel("Sending Voltage |Vs| (kV)")
ax_pv.set_title("PV Curve")
ax_pv.grid()
st.pyplot(fig_pv)

# ======================================
# KPIs
# ======================================
st.subheader("System Status")
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Ps (MW)", f"{Ps:.2f}")
c2.metric("Loss (MW)", f"{Loss:.2f}")
c3.metric("Efficiency %", f"{Efficiency:.2f}")
c4.metric("Loading %", f"{Loading:.2f}")
c5.metric("Stability", stability)
st.write("Operating Condition:", condition)
st.write("Stability Margin:", f"{Margin:.2f} MW")

# ======================================
# Power Circle with vectors
# ======================================

alpha = np.angle(A)
beta = np.angle(B)
ratio = np.abs(A/B)

Rr = ratio*(Vr_mag**2)
theta_r = beta - alpha
Crx = Rr*np.cos(theta_r)
Cry = Rr*np.sin(theta_r)

Vs_mag = np.abs(Vs)
Rs = ratio*(Vs_mag**2)
Csx = Rs*np.cos(beta)
Csy = Rs*np.sin(beta)

theta_plot = np.linspace(0, 2*np.pi, 400)
Pr_circle = Crx + Rr*np.cos(theta_plot)
Qr_circle = Cry + Rr*np.sin(theta_plot)
Ps_circle = Csx + Rs*np.cos(theta_plot)
Qs_circle = Csy + Rs*np.sin(theta_plot)

focus = max(abs(Pr), abs(Ps), abs(Qr), abs(Qs))
window = max(focus*3, 0.3*max(Rr, Rs))

# ======================================
# Stability Animation (Fixed + Lecture Style)
# ======================================
if run_anim:
    st.subheader("Operating Point Movement (Stability Animation)")

    load_values = np.linspace(0, Pmax_calc*1.2, 60)

    fig_anim, ax_anim = plt.subplots(figsize=(10,10))

    # Power circles
    ax_anim.plot(Pr_circle, Qr_circle, 'b', linewidth=2.5, label='Receiving Circle')
    ax_anim.plot(Ps_circle, Qs_circle, 'r', linewidth=2.5, label='Sending Circle')

    # ======================================
    # Stability Path (Dashed like lecture)
    # ======================================
    P_path = []
    Q_path = []

    for P_test in load_values:
        Q_test = P_test * np.tan(phi) if PF != 0 else 0
        P_path.append(P_test)
        Q_path.append(Q_test)

    ax_anim.plot(P_path, Q_path,
                 linestyle='--',
                 linewidth=2,
                 color='green',
                 label='Stability Path')

    # Mark operating point
    ax_anim.scatter(Pr, Qr, color='blue', s=120, zorder=6)
    ax_anim.scatter(Ps, Qs, color='red', s=120, zorder=6)

    # Pmax point
    ax_anim.scatter(Pmax_calc, 0, color='purple', s=120, zorder=7)
    ax_anim.text(Pmax_calc, 0,
                 f'  Pmax\n{Pmax_calc:.0f} MW',
                 color='purple', fontsize=12)

    # Operating vectors
    scale = 0.06 * max(Rr, Rs)
    ax_anim.arrow(0, 0, Pr, Qr,
                  head_width=scale,
                  length_includes_head=True,
                  color='blue')
    ax_anim.arrow(0, 0, Ps, Qs,
                  head_width=scale,
                  length_includes_head=True,
                  color='red')

    # Transfer line
    ax_anim.plot([Pr, Ps], [Qr, Qs],
                 'k', linewidth=2,
                 label="Power Transfer Line")

    # ======================================
    # Smart Zoom (Fixed indentation error)
    # ======================================
    max_P_anim = max(
        np.max(np.abs(Pr_circle)),
        np.max(np.abs(Ps_circle)),
        abs(Pmax_calc),
        abs(Pr),
        abs(Ps)
    )

    max_Q_anim = max(
        np.max(np.abs(Qr_circle)),
        np.max(np.abs(Qs_circle)),
        abs(Qr),
        abs(Qs)
    )

    window_anim = max(max_P_anim, max_Q_anim) * 1.2

    ax_anim.set_xlim(-window_anim, window_anim)
    ax_anim.set_ylim(-window_anim, window_anim)

    ax_anim.set_aspect('equal')
    ax_anim.grid(True, linestyle='--', alpha=0.4)
    ax_anim.set_xlabel("Active Power P (MW)")
    ax_anim.set_ylabel("Reactive Power Q (MVAR)")
    ax_anim.set_title("Stability Path (Lecture Style)")
    ax_anim.legend()

    st.pyplot(fig_anim)

# ======================================
# Professional Academic Power Circle (Lecture Level)
# ======================================
fig, ax = plt.subplots(figsize=(8,8))

# Circles (thinner for clarity)
ax.plot(Pr_circle, Qr_circle, 'b', linewidth=1.8, label="Receiving Circle")
ax.plot(Ps_circle, Qs_circle, 'r', linewidth=1.8, label="Sending Circle")

# Centers
ax.scatter(Crx, Cry, color='blue', s=40, zorder=6)
ax.scatter(Csx, Csy, color='red', s=40, zorder=6)
ax.text(Crx, Cry, " Cr", fontsize=12, color='blue')
ax.text(Csx, Csy, " Cs", fontsize=12, color='red')

# Line of centers (important in lecture)
ax.plot([Crx, Csx], [Cry, Csy],
        linestyle='--', color='black', linewidth=1.2)

# Origin
ax.scatter(0, 0, color='black', s=50)
ax.text(0, 0, " O", fontsize=12)

# ======================================
# Operating Points
# ======================================
ax.scatter(Pr, Qr, color='blue', s=120, edgecolors='black', zorder=7)
ax.scatter(Ps, Qs, color='red', s=120, edgecolors='black', zorder=7)

ax.text(Pr, Qr, "  k", fontsize=13, color='blue')
ax.text(Ps, Qs, "  k'", fontsize=13, color='red')

# ======================================
# Power Triangles (inside circles)
# ======================================
# Receiving triangle
ax.plot([0, Pr], [0, 0], 'b')
ax.plot([Pr, Pr], [0, Qr], 'b')
ax.plot([0, Pr], [0, Qr], 'b--')

# Sending triangle
ax.plot([0, Ps], [0, 0], 'r')
ax.plot([Ps, Ps], [0, Qs], 'r')
ax.plot([0, Ps], [0, Qs], 'r--')

# ======================================
# Power Factor Angles
# ======================================
phi_r = np.arctan2(Qr, Pr) if Pr != 0 else 0
phi_s = np.arctan2(Qs, Ps) if Ps != 0 else 0

arc_radius = 0.15 * window

theta_r = np.linspace(0, phi_r, 40)
ax.plot(arc_radius*np.cos(theta_r),
        arc_radius*np.sin(theta_r),
        color='blue', linewidth=2)

ax.text(arc_radius*1.2*np.cos(phi_r/2),
        arc_radius*1.2*np.sin(phi_r/2),
        "φr", fontsize=13, color='blue')

theta_s = np.linspace(0, phi_s, 40)
ax.plot(arc_radius*np.cos(theta_s),
        arc_radius*np.sin(theta_s),
        color='red', linewidth=2)

ax.text(arc_radius*1.2*np.cos(phi_s/2),
        arc_radius*1.2*np.sin(phi_s/2),
        "φs", fontsize=13, color='red')

# ======================================
# Stability Limit & Tangent at Pmax
# ======================================
ax.scatter(Pmax_calc, 0, color='purple', s=100, zorder=8)
ax.text(Pmax_calc, 0,
        f"  Pmax\n{Pmax_calc:.0f}",
        fontsize=12, color='purple')

# Tangent line (vertical approximation from lecture)
ax.plot([Pmax_calc, Pmax_calc],
        [-window, window],
        linestyle='--',
        color='purple',
        linewidth=1.5,
        label="Stability Limit")

# ======================================
# Transfer Line
# ======================================
ax.plot([Pr, Ps], [Qr, Qs], 'g', linewidth=2)

# Reactive limits
ax.axhline(Qmax, linestyle='--', color='orange', linewidth=1)
ax.axhline(Qmin, linestyle='--', color='orange', linewidth=1)

# ======================================
# Smart Scaling (important for clarity)
# ======================================
max_P = max(abs(Pr_circle).max(),
            abs(Ps_circle).max(),
            abs(Pmax_calc),
            abs(Pr),
            abs(Ps))

max_Q = max(abs(Qr_circle).max(),
            abs(Qs_circle).max(),
            abs(Qr),
            abs(Qs))

window = max(max_P, max_Q) * 1.1

ax.set_xlim(-window, window)
ax.set_ylim(-window, window)

ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.4)
ax.set_xlabel("P (MW)", fontsize=12)
ax.set_ylabel("Q (MVAR)", fontsize=12)
ax.set_title("Power Circle Diagram (Lecture Style)", fontsize=15)
ax.legend()

st.pyplot(fig)


# ======================================
# Voltage Along Line
# ======================================
st.subheader("Voltage Profile")

line_length = st.sidebar.number_input("Line Length for Profile (km)", value=100.0)
distance = np.linspace(0, line_length, 100)
V_profile = np.abs(Vr) + (np.abs(Vs)-np.abs(Vr))*(distance/line_length)

fig2, ax2 = plt.subplots(figsize=(10,4))
ax2.plot(distance, V_profile, linewidth=2)
ax2.set_xlabel("Distance (km)")
ax2.set_ylabel("Voltage (kV)")
ax2.set_title("Voltage Along the Line")
ax2.grid(True, linestyle='--', alpha=0.4)

st.pyplot(fig2)








