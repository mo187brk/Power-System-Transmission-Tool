import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch

# ======================================
# 1. Page Config & Dark Dashboard Styling
# ======================================
st.set_page_config(page_title="Power Transmission Analysis", layout="wide", page_icon="⚡")

st.markdown("""
    <style>
    .stApp { background: linear-gradient(135deg, #151928 0%, #0b0f19 100%); color: #ffffff; }
    h1, h2, h3, h4 { color: #ffffff !important; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
    div[data-testid="stWidgetLabel"] p, label p, .stSlider label p, .stNumberInput label p, .stSelectbox label p {
        color: #ffffff !important; font-weight: 600 !important; font-size: 14px !important;
    }
    div[data-testid="metric-container"] {
        background-color: #1e2233; padding: 15px; border-radius: 12px;
        box-shadow: 0 8px 24px #000000; border: 1px solid #2d3246; transition: transform 0.3s ease;
    }
    div[data-testid="metric-container"]:hover { transform: translateY(-5px); box-shadow: 0 12px 32px #00d2ff33; }
    div[data-testid="stMetricLabel"] p { color: #a0a5b5 !important; font-size: 15px !important; font-weight: 600 !important;}
    div[data-testid="stMetricValue"] { color: #ffffff !important; font-weight: 800 !important; font-size: 1.8rem !important;}
    div[data-testid="stButton"] button {
        background: linear-gradient(90deg, #d83cb8 0%, #f87462 100%) !important;
        color: #ffffff !important; border-radius: 8px !important; border: none !important;
        font-weight: 700 !important; font-size: 15px !important; height: 44px !important;
        width: 100% !important; margin-top: 27px !important;
        box-shadow: 0 4px 15px rgba(216, 60, 184, 0.4) !important; transition: all 0.3s ease !important;
    }
    div[data-testid="stButton"] button:hover {
        box-shadow: 0 6px 20px rgba(216, 60, 184, 0.6) !important;
        transform: translateY(-2px) !important;
        background: linear-gradient(90deg, #f87462 0%, #d83cb8 100%) !important;
    }
    hr { border-color: #2d3246; margin-top: 15px; margin-bottom: 15px; }
    </style>
""", unsafe_allow_html=True)

# ======================================
# Helper Functions
# ======================================
def polar_str(value, unit=""):
    mag = np.abs(value)
    ang = np.degrees(np.angle(value))
    return f"{mag:.4f} ∠ {ang:.2f}° {unit}"

def style_dark_plot(ax, fig, title, xlabel, ylabel):
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)
    ax.tick_params(colors='#ffffff')
    for spine in ax.spines.values():
        spine.set_color('#444a5e')
    ax.grid(True, linestyle=':', color='#444a5e', alpha=0.6)
    ax.set_title(title, fontsize=14, fontweight='bold', color='#ffffff', pad=15)
    ax.set_xlabel(xlabel, color='#ffffff', fontweight='bold')
    ax.set_ylabel(ylabel, color='#ffffff', fontweight='bold')

def fix_legend_color(legend):
    if legend:
        for text in legend.get_texts():
            text.set_color('white')

# ======================================
# 2. Session State Management
# ======================================
if 'history' not in st.session_state:
    st.session_state.history = []

# ======================================
# 3. Header & Control Panel
# ======================================
_, col_title = st.columns([1, 11])
with col_title:
    st.title("⚡ Power Transmission Analysis")

st.markdown("<hr>", unsafe_allow_html=True)

st.markdown("<span style='color:#00d2ff; font-weight:bold;'>1. Operating Conditions</span>", unsafe_allow_html=True)
c1, c2, c3, c4 = st.columns(4)
with c1:
    mode = st.selectbox("Operating Mode", ["Load", "No Load"])
with c2:
    Pr = st.number_input("Receiving Power Pr (MW)", value=100.0) if mode == "Load" else 0.0
with c3:
    PF = st.slider("Power Factor", 0.0, 1.0, 0.9) if mode == "Load" else 0.0
with c4:
    pf_type = st.selectbox("PF Type", ["Lagging", "Leading"]) if mode == "Load" else "Lagging"

st.markdown("<br><span style='color:#00d2ff; font-weight:bold;'>2. Line Parameters</span>", unsafe_allow_html=True)
c5, c6, c7, c8 = st.columns(4)
with c5:
    method = st.selectbox("ABCD Input Method", ["Enter A,B,C directly", "Enter Z,Y,length", "Enter r,L,C,length"])
    line_length = st.number_input("Line Length (km) [For Profile]", value=100.0)

with c6:
    if method == "Enter A,B,C directly":
        A_mag = st.number_input("|A|", value=0.986)
        A_ang = st.number_input("Angle A (°)", value=0.1)
    elif method == "Enter Z,Y,length":
        R = st.number_input("R (ohm/km)", value=0.1)
        X = st.number_input("X (ohm/km)", value=0.4)
    else:
        r = st.number_input("r (ohm/km)", value=0.1)
        L = st.number_input("L (H/km)", value=1e-3, format="%.5f")

with c7:
    if method == "Enter A,B,C directly":
        B_mag = st.number_input("|B|", value=60.0)
        B_ang = st.number_input("Angle B (°)", value=83.0)
    elif method == "Enter Z,Y,length":
        G = st.number_input("G (S/km)", value=0.0)
        B_sh = st.number_input("B (S/km)", min_value=0.0, value=0.000001, format="%.6f")
    else:
        C_val = st.number_input("C (F/km)", value=1e-8, format="%.9f")

with c8:
    if method == "Enter A,B,C directly":
        C_mag = st.number_input("|C|", value=0.0004, format="%.5f")
        C_ang = st.number_input("Angle C (°)", value=90.0)

st.markdown("<br><span style='color:#00d2ff; font-weight:bold;'>3. Constraints, Faults & Actions</span>", unsafe_allow_html=True)
c9, c10, c11, c12 = st.columns(4)
with c9:
    Vr_mag = st.number_input("Receiving Voltage Vr (kV) [Line-to-Line]", value=220.0)
    run_anim = st.checkbox("Run Stability Animation", value=False)
with c10:
    fault = st.selectbox("Fault Simulation", ["No Fault", "3-Phase Fault at Receiving End"])
    Qmax = st.number_input("Qmax (MVAR)", value=300.0)
with c11:
    fault_impedance = st.number_input("Fault Impedance (Ω)", value=0.01)
    Qmin = st.number_input("Qmin (MVAR)", value=-300.0)
with c12:
    save_flag = st.button("💾 Save Condition")

# ======================================
# 4. Core Mathematical Calculations (CORRECTED)
# ======================================

# --- Power angle ---
phi = np.arccos(PF) if PF > 0 else 0.0
if pf_type == "Leading":
    phi = -phi

# --- Total 3-phase reactive power ---
Qr = Pr * np.tan(phi) if PF > 0 else 0.0

# --- ABCD Parameters ---
Z, Y = None, None
if method == "Enter A,B,C directly":
    A = A_mag * np.exp(1j * np.radians(A_ang))
    B = B_mag * np.exp(1j * np.radians(B_ang))
    C = C_mag * np.exp(1j * np.radians(C_ang))
    D = A
elif method == "Enter Z,Y,length":
    Z = (R + 1j * X) * line_length
    Y = (G + 1j * B_sh) * line_length
    theta = np.sqrt(Z * Y)
    A = np.cosh(theta)
    B = Z * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Z
    C = Y * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Y
    D = A
else:
    w = 2 * np.pi * 50
    Z = (r + 1j * w * L) * line_length
    Y = (1j * w * C_val) * line_length
    theta = np.sqrt(Z * Y)
    A = np.cosh(theta)
    B = Z * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Z
    C = Y * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Y
    D = A

# =====================================================================
# FIX: Use per-phase voltage (line-to-line / sqrt(3)) for current calc
# All powers are 3-phase totals (MW, MVAR); voltages are LINE-TO-LINE (kV)
# Ir calculation uses per-phase voltage to get correct per-phase current
# =====================================================================
Vr_phase = Vr_mag / np.sqrt(3)   # Phase voltage (kV)
Vr = Vr_phase + 0j               # Reference phasor (real axis)

# --- Per-phase apparent power (MVA) ---
Sr_phase = (Pr / 3) + 1j * (Qr / 3)   # MVA per phase

# --- Per-phase receiving current (kA) ---
if np.abs(Vr) > 1e-10:
    Ir = np.conj(Sr_phase) / np.conj(Vr)   # kA per phase
else:
    Ir = 0.0 + 0j

# --- Fault simulation ---
if fault == "3-Phase Fault at Receiving End":
    Vr = fault_impedance * Ir   # Vr drops to fault voltage

# --- Sending end (phase quantities) ---
Vs_phase = A * Vr + B * Ir       # kV phase
Is_phase = C * Vr + D * Ir       # kA per phase

# --- Line-to-line voltages ---
Vs_mag = np.abs(Vs_phase) * np.sqrt(3)   # kV line-to-line

# --- 3-phase powers ---
Ss_3ph = 3 * Vs_phase * np.conj(Is_phase)   # MVA
Ps = np.real(Ss_3ph)   # MW
Qs = np.imag(Ss_3ph)   # MVAR

Loss = Ps - Pr
Efficiency = (Pr / Ps) * 100 if Ps > 1e-6 else 0.0
Voltage_reg = (Vs_mag - Vr_mag) / Vr_mag * 100

# =====================================================================
# FIX: Power Circle Diagram calculations using LINE voltages
# Standard formula:
#   Centre of RPCD: n_r = (A*Vr²/B) at angle -(β-α) from origin
#   Radius Rr = Vr * Vs_line / |B|
#   Centre of SPCD: n_s = (A*Vs²/B) at angle +(β-α) from origin  
#   Radius Rs = Vr * Vs_line / |B|
# =====================================================================
A_magnitude = np.abs(A)
B_magnitude = np.abs(B)
beta  = np.angle(B, deg=True)    # β in degrees
alpha = np.angle(A, deg=True)    # α in degrees
delta = np.angle(Vs_phase, deg=True)   # power angle δ (Vs leads Vr)

# Pmax (receiving): maximum receivable power
Pmax_calc = (Vs_mag * Vr_mag) / B_magnitude  # MW (3-phase)
# Shift by centre offset:
Pmax_calc = Pmax_calc - (A_magnitude * Vr_mag**2 / B_magnitude) * np.cos(np.radians(beta - alpha))

Loading  = (Pr / Pmax_calc) * 100 if Pmax_calc > 1e-6 else 0.0
Margin   = Pmax_calc - Pr
stability = "Stable" if Pr < Pmax_calc else "Unstable"

# ======================================
# 5. Save Logic
# ======================================
if save_flag:
    st.session_state.history.append({
        "Mode": mode, "Pr (MW)": round(Pr, 2), "PF": round(PF, 3),
        "Ps (MW)": round(Ps, 2), "Qs (MVAR)": round(Qs, 2),
        "Loading %": round(Loading, 1), "Stability": stability
    })

# ======================================
# 6. Top KPIs (Metrics)
# ======================================
st.markdown("<hr>", unsafe_allow_html=True)
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Sending Power (Ps)", f"{Ps:.2f} MW")
m2.metric("Power Loss", f"{Loss:.2f} MW")
m3.metric("Efficiency", f"{Efficiency:.2f} %")
m4.metric("Line Loading", f"{Loading:.2f} %")

stability_color = "#5cb85c" if stability == "Stable" else "#d9534f"
m5.markdown(f"""
    <div data-testid="metric-container" style="border-left: 5px solid {stability_color};">
        <label style="color: #a0a5b5; font-size: 15px; font-weight: 600;">System Stability</label>
        <div style="color: {stability_color}; font-weight: 800; font-size: 1.8rem; margin-top: -5px;">{stability}</div>
    </div>
""", unsafe_allow_html=True)

# ======================================
# 7. Tabs
# ======================================
st.markdown("<br>", unsafe_allow_html=True)
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Power Circle Diagram", "📈 PV Curve",
    "⚡ Voltage Profile", "📋 Detailed Status", "💾 Saved Cases"
])

# ============================================================
# TAB 1 — Power Circle Diagram (CORRECTED like lecture)
# ============================================================
with tab1:
    col_plot, col_info = st.columns([3, 1])
    with col_plot:

        # -------------------------------------------------------
        # LECTURE-CORRECT PCD GEOMETRY
        # Using line-to-line voltages (kV) throughout
        #
        # RPCD centre  n_r:
        #   n_rx = -(A|Vr|²/|B|) * cos(β-α)
        #   n_ry = -(A|Vr|²/|B|) * sin(β-α)
        #
        # SPCD centre  n_s:
        #   n_sx = +(A|Vs|²/|B|) * cos(β-α)    [sign convention from lecture]
        #   n_sy = -(A|Vs|²/|B|) * sin(β-α)    → actually +(β-δ) side
        #
        # Radius (both circles) = |Vr|*|Vs| / |B|
        #
        # Standard formulas (3-phase, line kV):
        #   Pr = -A|Vr|²/|B| * cos(β-α) + |Vr||Vs|/|B| * cos(β-δ)
        #   Qr = -A|Vr|²/|B| * sin(β-α) + |Vr||Vs|/|B| * sin(β-δ)
        #   Ps = +A|Vs|²/|B| * cos(β-α) - |Vr||Vs|/|B| * cos(β+δ)
        #   Qs = +A|Vs|²/|B| * sin(β-α) - |Vr||Vs|/|B| * sin(β+δ)
        # -------------------------------------------------------

        beta_r  = np.radians(beta)
        alpha_r = np.radians(alpha)
        delta_r = np.radians(delta)

        # RPCD centre (n_r)
        n_rx = -(A_magnitude * Vr_mag**2 / B_magnitude) * np.cos(beta_r - alpha_r)
        n_ry = -(A_magnitude * Vr_mag**2 / B_magnitude) * np.sin(beta_r - alpha_r)

        # SPCD centre (n_s)  — note positive sign from lecture
        n_sx = (A_magnitude * Vs_mag**2 / B_magnitude) * np.cos(beta_r - alpha_r)
        n_sy = (A_magnitude * Vs_mag**2 / B_magnitude) * np.sin(beta_r - alpha_r)

        # Radius (same for both in this network)
        R_circle = (Vr_mag * Vs_mag) / B_magnitude

        # Circle points
        theta_arr = np.linspace(0, 2 * np.pi, 400)
        Pr_circle = n_rx + R_circle * np.cos(theta_arr)
        Qr_circle = n_ry + R_circle * np.sin(theta_arr)
        Ps_circle = n_sx + R_circle * np.cos(theta_arr)
        Qs_circle = n_sy + R_circle * np.sin(theta_arr)

        # Pmax on RPCD (operating angle α_op = β-α for RPCD,
        # Pmax is the rightmost point = n_rx + R_circle)
        Pmax_rpcd = n_rx + R_circle

        if run_anim:
            # --- Stability animation ---
            st.markdown(
                "<h5 style='color:white; text-align:center;'>Operating Point Movement (Stability Animation)</h5>",
                unsafe_allow_html=True
            )
            load_values = np.linspace(0, Pmax_rpcd * 1.2, 40)
            fig_anim, ax_anim = plt.subplots(figsize=(8, 8))
            style_dark_plot(ax_anim, fig_anim, "", "Active Power P (MW)", "Reactive Power Q (MVAR)")

            ax_anim.plot(Pr_circle, Qr_circle, '#00d2ff', linewidth=2, label='Receiving (RPCD)', alpha=0.7)
            ax_anim.plot(Ps_circle, Qs_circle, '#ff8d72', linewidth=2, label='Sending (SPCD)', alpha=0.7)

            for P_test in load_values:
                Q_test = P_test * np.tan(phi) if PF > 0 else 0
                color = '#5cb85c' if P_test <= Pmax_rpcd else '#d9534f'
                ax_anim.scatter(P_test, Q_test, color=color, s=25, zorder=5)

            ax_anim.scatter(Pmax_rpcd, 0, color='#e14eca', s=120, zorder=6)
            ax_anim.text(Pmax_rpcd, 0, f' Pmax\n{Pmax_rpcd:.0f}', color='#e14eca', fontsize=12)

            focus = max(abs(Pmax_rpcd), R_circle) * 1.5
            ax_anim.set_xlim(-focus, focus * 1.3)
            ax_anim.set_ylim(-focus, focus)
            ax_anim.axhline(0, color='#626b82', linewidth=1.2)
            ax_anim.axvline(0, color='#626b82', linewidth=1.2)
            ax_anim.set_aspect('equal')
            leg = ax_anim.legend(facecolor='#1e2233', edgecolor='#444a5e')
            fix_legend_color(leg)
            st.pyplot(fig_anim)

        else:
            # -------------------------------------------------------
            # MAIN PCD PLOT — lecture style
            # -------------------------------------------------------
            fig, ax = plt.subplots(figsize=(9, 9))
            style_dark_plot(ax, fig, "Combined Power Circle Diagram", "Active Power P (MW)", "Reactive Power Q (MVAR)")

            # --- Draw the two circles ---
            ax.plot(Pr_circle, Qr_circle, color='#00d2ff', linewidth=2.0, label='RPCD (Receiving)', zorder=3)
            ax.plot(Ps_circle, Qs_circle, color='#ff8d72', linewidth=2.0, label='SPCD (Sending)', zorder=3)

            # --- Qmax / Qmin limits ---
            ax.axhline(Qmax, linestyle='--', color='#e14eca', alpha=0.6, label=f"Qmax = {Qmax:.0f}")
            ax.axhline(Qmin, linestyle='--', color='#e14eca', alpha=0.6, label=f"Qmin = {Qmin:.0f}")

            # --- Centre points ---
            ax.scatter(n_rx, n_ry, color='#00d2ff', s=60, marker='x', zorder=6, linewidths=2)
            ax.scatter(n_sx, n_sy, color='#ff8d72', s=60, marker='x', zorder=6, linewidths=2)
            ax.text(n_rx + 2, n_ry + 2, 'n_r', color='#00d2ff', fontsize=10)
            ax.text(n_sx + 2, n_sy + 2, 'n_s', color='#ff8d72', fontsize=10)

            # --- Radius vectors from centre to operating point ---
            # RPCD: n_r → Sr
            ax.annotate("", xy=(Pr, Qr), xytext=(n_rx, n_ry),
                        arrowprops=dict(arrowstyle='->', color='#00d2ff', lw=1.5))
            # SPCD: n_s → Ss
            ax.annotate("", xy=(Ps, Qs), xytext=(n_sx, n_sy),
                        arrowprops=dict(arrowstyle='->', color='#ff8d72', lw=1.5))

            # --- Vectors from ORIGIN (0,0) to Sr and Ss ---
            ax.annotate("", xy=(Pr, Qr), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->', color='#00d2ff', lw=2.5))
            ax.annotate("", xy=(Ps, Qs), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='->', color='#ff8d72', lw=2.5))

            # --- Sr, Ss labels ---
            ax.scatter(Pr, Qr, color='#00d2ff', s=120, edgecolors='white', zorder=7)
            ax.scatter(Ps, Qs, color='#ff8d72', s=120, edgecolors='white', zorder=7)
            ax.text(Pr + 2, Qr + 2, 'Sr', fontsize=13, color='#00d2ff', fontweight='bold')
            ax.text(Ps + 2, Qs - 8, 'Ss', fontsize=13, color='#ff8d72', fontweight='bold')

            # --- Power Transfer line Sr → Ss ---
            ax.plot([Pr, Ps], [Qr, Qs], color='#5cb85c', linewidth=1.8,
                    linestyle='--', label='Power Transfer Sr→Ss', zorder=4)

            # --- Pmax point ---
            ax.scatter(Pmax_rpcd, 0, color='#e14eca', s=140, zorder=8, marker='^')
            ax.text(Pmax_rpcd + 2, 3, f'Pmax\n{Pmax_rpcd:.0f} MW',
                    color='#e14eca', fontsize=11, fontweight='bold')

            # --- φr angle arc (from +x axis to Sr vector) ---
            phi_deg = np.degrees(np.arctan2(Qr, Pr))
            arc_r = R_circle * 0.15
            arc_angles = np.linspace(0, np.radians(phi_deg), 40)
            ax.plot(arc_r * np.cos(arc_angles), arc_r * np.sin(arc_angles),
                    color='#00d2ff', linewidth=1.0, linestyle=':')
            ax.text(arc_r * 1.1, arc_r * 0.3,
                    f'φr = {abs(phi_deg):.1f}°', color='#00d2ff', fontsize=10)

            # --- δ angle label ---
            ax.text(0.02 * R_circle, -0.08 * R_circle,
                    f'δ = {delta:.1f}°',
                    color='white', fontsize=11,
                    bbox=dict(facecolor='#1e2233', alpha=0.8, edgecolor='#444a5e', boxstyle='round,pad=0.3'))

            # --- β - α angle annotation ---
            ax.text(n_rx + R_circle * 0.05, n_ry - R_circle * 0.12,
                    f'β−α = {beta - alpha:.1f}°', color='#aaaaaa', fontsize=9)

            # --- Axes lines ---
            ax.axhline(0, color='#626b82', linewidth=1.2)
            ax.axvline(0, color='#626b82', linewidth=1.2)

            # --- Equal aspect & window ---
            all_vals = [abs(Pr), abs(Qr), abs(Ps), abs(Qs), R_circle, abs(n_rx), abs(n_ry), Pmax_rpcd]
            window = max(all_vals) * 1.4 if all_vals else 100
            ax.set_xlim(-window * 0.6, window * 1.2)
            ax.set_ylim(-window * 0.8, window * 0.8)
            ax.set_aspect('equal', adjustable='box')

            # --- Legend ---
            leg = ax.legend(facecolor='#1e2233', edgecolor='#444a5e', fontsize=9, loc='lower right')
            fix_legend_color(leg)

            st.pyplot(fig)

    with col_info:
        st.markdown("#### Point Details")
        st.info(f"**Pr:** {Pr:.2f} MW")
        st.info(f"**Qr:** {Qr:.2f} MVAR")
        st.info(f"**Ps:** {Ps:.2f} MW")
        st.info(f"**Qs:** {Qs:.2f} MVAR")
        st.info(f"**Pmax:** {Pmax_rpcd:.2f} MW")
        st.info(f"**Margin:** {Margin:.2f} MW")
        st.info(f"**δ:** {delta:.2f}°")
        st.info(f"**β−α:** {beta - alpha:.2f}°")

# ============================================================
# TAB 2 — PV Curve
# ============================================================
with tab2:
    st.markdown("#### PV Curve (Voltage Stability)")
    P_range = np.linspace(0, Pmax_rpcd * 1.2, 60)
    V_list = []
    for P_test in P_range:
        Q_test = P_test * np.tan(phi) if PF > 0 else 0
        Sr_test_ph = (P_test / 3) + 1j * (Q_test / 3)
        Ir_test = np.conj(Sr_test_ph) / np.conj(Vr) if np.abs(Vr) > 1e-10 else 0
        Vs_test = A * Vr + B * Ir_test
        V_list.append(np.abs(Vs_test) * np.sqrt(3))   # line voltage

    fig_pv, ax_pv = plt.subplots(figsize=(10, 4))
    style_dark_plot(ax_pv, fig_pv, "PV Curve", "Load P (MW)", "Sending Voltage |Vs| line-to-line (kV)")
    ax_pv.plot(P_range, V_list, color='#e14eca', linewidth=2)
    ax_pv.axvline(Pmax_rpcd, linestyle='--', color='#ff8d72', label=f"Pmax ({Pmax_rpcd:.1f} MW)")
    ax_pv.scatter(Pr, Vs_mag, color='#00d2ff', s=100, zorder=5, label=f"Operating point")
    leg_pv = ax_pv.legend(facecolor='#1e2233', edgecolor='#444a5e')
    fix_legend_color(leg_pv)
    st.pyplot(fig_pv)

# ============================================================
# TAB 3 — Voltage Profile
# ============================================================
with tab3:
    st.markdown("#### Voltage Profile (along the line)")
    Z_total = Z if Z is not None else B
    Y_total = Y if Y is not None else C

    if line_length > 0 and Z_total is not None and Y_total is not None:
        Z_per_km = Z_total / line_length
        Y_per_km = Y_total / line_length
        gamma = np.sqrt(Z_per_km * Y_per_km)
        Zc = np.sqrt(Z_per_km / Y_per_km) if np.abs(Y_per_km) > 1e-15 else 1.0

        distance = np.linspace(0, line_length, 200)
        V_profile = []
        for x in distance:
            Vx = Vr * np.cosh(gamma * x) + Ir * Zc * np.sinh(gamma * x)
            V_profile.append(np.abs(Vx) * np.sqrt(3))   # line voltage

        fig2, ax2 = plt.subplots(figsize=(10, 4))
        style_dark_plot(ax2, fig2,
                        "Voltage Profile along the Transmission Line",
                        "Distance from Receiving End (km)",
                        "Voltage (kV) line-to-line")
        ax2.plot(distance, V_profile, color='#00d2ff', linewidth=2)
        ax2.invert_xaxis()
        st.pyplot(fig2)
    else:
        st.info("Voltage profile requires Z,Y parameters. Please use 'Enter Z,Y,length' or 'Enter r,L,C,length'.")

# ============================================================
# TAB 4 — Detailed Status
# ============================================================
with tab4:
    st.markdown("#### 📋 System Detailed Status")
    c_s1, c_s2 = st.columns(2)
    with c_s1:
        st.markdown("**Voltages & Currents (Phase quantities)**")
        st.code(
            f"Vr (phase) = {polar_str(Vr, 'kV')}\n"
            f"Ir (phase) = {polar_str(Ir, 'kA')}\n"
            f"Vs (phase) = {polar_str(Vs_phase, 'kV')}\n"
            f"Is (phase) = {polar_str(Is_phase, 'kA')}\n"
            f"Vs (line)  = {Vs_mag:.4f} kV\n"
            f"Vr (line)  = {Vr_mag:.4f} kV"
        )
        st.markdown("**Performance Metrics**")
        reactive_status = (
            'Overexcited' if Qs > Qmax
            else 'Underexcited' if Qs < Qmin
            else 'Within Limits'
        )
        st.code(
            f"Power angle δ  = {delta:.4f}°\n"
            f"Voltage Reg    = {Voltage_reg:.4f} %\n"
            f"Efficiency     = {Efficiency:.4f} %\n"
            f"Reactive Status= {reactive_status}"
        )
    with c_s2:
        st.markdown("**ABCD Parameters**")
        st.code(
            f"A = {polar_str(A)}\n"
            f"B = {polar_str(B, 'Ω')}\n"
            f"C = {polar_str(C, 'S')}\n"
            f"D = {polar_str(D)}"
        )
        st.markdown("**PCD Key Values**")
        st.code(
            f"n_r  centre = ({n_rx:.2f}, {n_ry:.2f}) MW/MVAR\n"
            f"n_s  centre = ({n_sx:.2f}, {n_sy:.2f}) MW/MVAR\n"
            f"Radius Rr   = {R_circle:.2f} MVA\n"
            f"Pmax (RPCD) = {Pmax_rpcd:.2f} MW\n"
            f"β − α       = {beta - alpha:.2f}°\n"
            f"A|Vr|²/|B|  = {A_magnitude * Vr_mag**2 / B_magnitude:.2f} MVA"
        )

# ============================================================
# TAB 5 — Saved Cases
# ============================================================
with tab5:
    st.markdown("#### 💾 Saved Cases Management")
    if len(st.session_state.history) > 0:
        headers = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
        headers[0].markdown("**Mode**");    headers[1].markdown("**Pr (MW)**")
        headers[2].markdown("**PF**");      headers[3].markdown("**Ps (MW)**")
        headers[4].markdown("**Qs (MVAR)**"); headers[5].markdown("**Load %**")
        headers[6].markdown("**Stability**"); headers[7].markdown("**Action**")
        st.markdown("<hr style='margin: 5px 0;'>", unsafe_allow_html=True)

        for i, case in enumerate(st.session_state.history):
            cols = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
            cols[0].write(case["Mode"])
            cols[1].write(f"{case['Pr (MW)']}")
            cols[2].write(f"{case['PF']}")
            cols[3].write(f"{case['Ps (MW)']}")
            cols[4].write(f"{case['Qs (MVAR)']}")
            cols[5].write(f"{case['Loading %']}")
            s_color = "#5cb85c" if case['Stability'] == "Stable" else "#d9534f"
            cols[6].markdown(
                f"<span style='color:{s_color}; font-weight:bold;'>{case['Stability']}</span>",
                unsafe_allow_html=True
            )
            if cols[7].button("🗑️", key=f"del_{i}"):
                st.session_state.history.pop(i)
                st.rerun()
    else:
        st.info("No cases saved yet. Adjust parameters and click 'Save Condition'.")
