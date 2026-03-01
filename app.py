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
# الرسمة الاحترافية (مطابقة للمحاضرة)
# ======================================

# حل مشكلة الـ NameError: بنعرف D بتساوي A لأن الخط متماثل
D = A 

# 1. الحسابات الهندسية للمراكز (Centers)
# زوايا ثوابت الخط
alpha_rad = np.angle(A)
beta_rad = np.angle(B)
beta_deg = np.degrees(beta_rad)
alpha_deg = np.degrees(alpha_rad)

# مركز دائرة الاستقبال (نقطة n) -> الإحداثيات: -(A/B)*Vr^2
n_mag = (np.abs(A) / np.abs(B)) * (Vr_mag**2)
nx, ny = -n_mag * np.cos(beta_rad - alpha_rad), n_mag * np.sin(beta_rad - alpha_rad)

# مركز دائرة الإرسال (نقطة n') -> الإحداثيات: (D/B)*Vs^2
ns_mag = (np.abs(D) / np.abs(B)) * (Vs_mag**2)
nsx, nsy = ns_mag * np.cos(beta_rad - alpha_rad), -ns_mag * np.sin(beta_rad - alpha_rad)

fig, ax = plt.subplots(figsize=(12, 12))

# 2. رسم الدوائر (نفس الحسابات الأصلية عندك بس بشكل منقط)
ax.plot(Pr_circle, Qr_circle, color='purple', linestyle='--', alpha=0.3)
ax.plot(Ps_circle, Qs_circle, color='purple', linestyle='--', alpha=0.3)

# 3. رسم المحاور (Active & Reactive Power) باللون الأحمر
ax.axhline(0, color='red', linewidth=1.2)
ax.axvline(0, color='red', linewidth=1.2)
ax.text(window*0.9, 5, 'Active power', color='black', fontsize=12, fontweight='bold')
ax.text(5, window*0.9, 'Reactive power', color='black', fontsize=12, fontweight='bold', rotation=90)

# 4. المتجهات الرئيسية (Vectors)
# المتجه من الأصل لـ n والمتجه من n لـ k
ax.annotate('', xy=(0, 0), xytext=(nx, ny), arrowprops=dict(arrowstyle='->', color='blue', lw=2))
ax.annotate('', xy=(Pr, Qr), xytext=(nx, ny), arrowprops=dict(arrowstyle='->', color='blue', lw=2))

# المتجه من الأصل لـ n' والمتجه من n' لـ k'
ax.annotate('', xy=(0, 0), xytext=(nsx, nsy), arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
ax.annotate('', xy=(Ps, Qs), xytext=(nsx, nsy), arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))

# 5. تسمية النقاط (Labels) كما في المحاضرة
ax.text(nx, ny, ' n', fontsize=16, color='blue', fontweight='bold', va='bottom')
ax.text(Pr, Qr, ' k', fontsize=16, color='blue', fontweight='bold', va='top')
ax.text(nsx, nsy, " n'", fontsize=16, color='darkgreen', fontweight='bold', va='top')
ax.text(Ps, Qs, " k'", fontsize=16, color='darkgreen', fontweight='bold', va='bottom')

# تسمية المتجهات (القيم)
ax.text(nx/2, ny/2, r'$\frac{A}{B}V_r^2$', fontsize=14, color='blue', ha='right')
ax.text((nx+Pr)/2, (ny+Qr)/2, r'$\frac{V_s V_r}{B}$', fontsize=14, color='blue')

# 6. رسم الزوايا (Arcs)
# زاوية beta - alpha
from matplotlib.patches import Arc
angle_ba = Arc((0, 0), n_mag*0.4, n_mag*0.4, theta1=180-np.degrees(beta_rad-alpha_rad), theta2=180, color='red', lw=1.5)
ax.add_patch(angle_ba)
ax.text(nx*0.2, ny*0.1, r'$\beta-\alpha$', color='red', fontsize=13)

# زاوية delta (بين k و الخط العمودي لـ n)
delta_val = np.degrees(np.angle(Vs) - np.angle(Vr))
angle_delta = Arc((nx, ny), n_mag*0.5, n_mag*0.5, theta1=-np.degrees(beta_rad-alpha_rad), 
                  theta2=-np.degrees(beta_rad-alpha_rad)+delta_val, color='magenta', lw=2)
ax.add_patch(angle_delta)
ax.text(nx + (Pr-nx)*0.3, ny + (Qr-ny)*0.3, r'$\delta$', color='magenta', fontsize=15, fontweight='bold')

# 7. خط الربط (Power Transfer Line)
ax.plot([Pr, Ps], [Qr, Qs], 'g', linewidth=2, label="Transfer Line")

# 8. ضبط الحدود والشكل النهائي
ax.set_xlim(-window*1.5, window*1.5)
ax.set_ylim(-window*1.5, window*1.5)
ax.set_aspect('equal')
ax.grid(True, linestyle=':', alpha=0.3)
ax.set_title("Combined Sending and Receiving-end Power Circle Diagram", fontsize=16, pad=20)

st.pyplot(fig)
# ======================================
# Voltage Along Line
# ======================================
st.subheader("Voltage Profile")
line_length = st.sidebar.number_input("Line Length for Profile (km)", value=100.0)
distance = np.linspace(0, line_length, 100)
V_profile = np.abs(Vr) + (np.abs(Vs)-np.abs(Vr))*(distance/line_length)

fig2, ax2 = plt.subplots()
ax2.plot(distance, V_profile)
ax2.set_xlabel("Distance (km)")
ax2.set_ylabel("Voltage (kV)")
ax2.set_title("Voltage Along the Line")
ax2.grid()
st.pyplot(fig2)

