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

# --- الجزء الخاص بالرسم المتوافق مع المحاضرة ---

fig, ax = plt.subplots(figsize=(12, 12))

# 1. إعدادات الزوايا والمراكز (بناءً على الكود الأصلي)
alpha_rad = np.angle(A)
beta_rad = np.angle(B)
beta_deg = np.degrees(beta_rad)
alpha_deg = np.degrees(alpha_rad)

# مركز دائرة الاستقبال (نقطة n)
# n = (A/B) * Vr^2 بزاوية (beta - alpha) في الربع الثاني/الثالث هندسياً
# ملاحظة: في الرسمة التقليدية يتم عكس الإشارات لتناسب اتجاه القدرة
n_mag = (np.abs(A) / np.abs(B)) * (Vr_mag**2)
nx, ny = -n_mag * np.cos(beta_rad - alpha_rad), n_mag * np.sin(beta_rad - alpha_rad)

# مركز دائرة الإرسال (نقطة n')
ns_mag = (np.abs(D) / np.abs(B)) * (Vs_mag**2) # D غالباً تساوي A
nsx, nsy = n_mag * np.cos(beta_rad - alpha_rad), -n_mag * np.sin(beta_rad - alpha_rad)

# 2. رسم المحاور الأساسية (Active & Reactive Power)
ax.axhline(0, color='red', linewidth=1.5) # محور P
ax.axvline(0, color='red', linewidth=1.5) # محور Q
ax.text(window*0.8, 10, 'Active power', fontsize=12, fontweight='bold')
ax.text(10, window*0.8, 'Reactive power', fontsize=12, fontweight='bold', rotation=90)

# 3. رسم الدوائر المنقطة (Circles)
ax.plot(Pr_circle, Qr_circle, color='purple', linestyle='--', alpha=0.5, label='Receiving Circle')
ax.plot(Ps_circle, Qs_circle, color='purple', linestyle='--', alpha=0.5, label='Sending Circle')

# 4. رسم المتجهات الرئيسية (Vectors as lines with arrows)
# متجه n (المركز)
ax.annotate('', xy=(0, 0), xytext=(nx, ny), arrowprops=dict(arrowstyle='->', color='blue', lw=2))
ax.text(nx-20, ny+20, 'n', fontsize=15, color='blue', fontweight='bold')

# متجه n' (مركز الإرسال)
ax.annotate('', xy=(0, 0), xytext=(nsx, nsy), arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
ax.text(nsx+10, nsy-20, "n'", fontsize=15, color='darkgreen', fontweight='bold')

# المتجه الواصل لنقطة التشغيل (Radius Vector)
ax.annotate('', xy=(Pr, Qr), xytext=(nx, ny), arrowprops=dict(arrowstyle='->', color='blue', lw=2))
ax.text(Pr+5, Qr+5, 'k', fontsize=15, color='blue', fontweight='bold')

ax.annotate('', xy=(Ps, Qs), xytext=(nsx, nsy), arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))
ax.text(Ps+5, Qs+5, "k'", fontsize=15, color='darkgreen', fontweight='bold')

# 5. إضافة الزوايا (Arcs) - لمحاكاة شكل المحاضرة
# زاوية beta - alpha
arc_angle = patches.Arc((0, 0), window*0.2, window*0.2, theta1=180-(beta_deg-alpha_deg), theta2=180, 
                        color='black', linewidth=1.5)
ax.add_patch(arc_angle)
ax.text(nx*0.3, ny*0.1, r'$\beta-\alpha$', fontsize=12)

# زاوية delta (بين نصف القطر والمحور الشاقولي لنقطة n)
delta_deg = np.degrees(np.angle(Vs) - np.angle(Vr))
arc_delta = patches.Arc((nx, ny), window*0.3, window*0.3, theta1=- (beta_deg-alpha_deg), theta2=- (beta_deg-alpha_deg) + delta_deg, 
                        color='red', linewidth=2)
ax.add_patch(arc_delta)
ax.text(nx+20, ny-40, r'$\delta$', color='red', fontsize=14, fontweight='bold')

# 6. خط نقل القدرة (الواصل بين k و k')
ax.plot([Pr, Ps], [Qr, Qs], 'orange', linewidth=3, label='Power Transfer')

# 7. تحسين المظهر العام (Scaling)
limit = max(abs(nx), abs(ny), abs(Pr), abs(Ps)) * 1.4
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_aspect('equal')
ax.grid(True, which='both', linestyle=':', alpha=0.4)

# إخفاء الإطارات الخارجية لجعلها تبدو كرسمة "سبورة"
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

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






