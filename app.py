import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ======================================
# 1. Page Config & Styling
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
        box-shadow: 0 4px 15px rgba(216, 60, 184, 0.4) !important;
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
# Helpers
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
    ax.grid(True, linestyle=':', color='#444a5e', alpha=0.5)
    ax.set_title(title, fontsize=14, fontweight='bold', color='#ffffff', pad=15)
    ax.set_xlabel(xlabel, color='#ffffff', fontweight='bold')
    ax.set_ylabel(ylabel, color='#ffffff', fontweight='bold')

def fix_legend_color(legend):
    if legend:
        for text in legend.get_texts():
            text.set_color('white')

# ======================================
# 2. Session State
# ======================================
if 'history' not in st.session_state:
    st.session_state.history = []

# ======================================
# 3. Header
# ======================================
_, col_title = st.columns([1, 11])
with col_title:
    st.title("⚡ Power Transmission Analysis")

st.markdown("<hr>", unsafe_allow_html=True)

# ======================================
# 4. Inputs
# ======================================
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
    line_length = st.number_input("Line Length (km)", value=100.0)

with c6:
    if method == "Enter A,B,C directly":
        A_mag = st.number_input("|A|", value=0.986)
        A_ang = st.number_input("Angle A (°)", value=0.1)
    elif method == "Enter Z,Y,length":
        R_ohm = st.number_input("R (ohm/km)", value=0.1)
        X_ohm = st.number_input("X (ohm/km)", value=0.4)
    else:
        r_ohm = st.number_input("r (ohm/km)", value=0.1)
        L_val = st.number_input("L (H/km)", value=1e-3, format="%.5f")

with c7:
    if method == "Enter A,B,C directly":
        B_mag = st.number_input("|B|", value=60.0)
        B_ang = st.number_input("Angle B (°)", value=83.0)
    elif method == "Enter Z,Y,length":
        G_val = st.number_input("G (S/km)", value=0.0)
        B_sh  = st.number_input("B (S/km)", min_value=0.0, value=1e-6, format="%.6f")
    else:
        C_val = st.number_input("C (F/km)", value=1e-8, format="%.9f")

with c8:
    if method == "Enter A,B,C directly":
        C_mag = st.number_input("|C|", value=0.0004, format="%.5f")
        C_ang = st.number_input("Angle C (°)", value=90.0)

st.markdown("<br><span style='color:#00d2ff; font-weight:bold;'>3. Constraints & Actions</span>", unsafe_allow_html=True)
c9, c10, c11, c12 = st.columns(4)
with c9:
    Vr_mag = st.number_input("Receiving Voltage Vr (kV) — Line-to-Line", value=220.0)
    run_anim = st.checkbox("Run Stability Animation", value=False)
with c10:
    fault = st.selectbox("Fault Simulation", ["No Fault", "3-Phase Fault at Receiving End"])
    Qmax  = st.number_input("Qmax (MVAR)", value=300.0)
with c11:
    fault_impedance = st.number_input("Fault Impedance (Ω)", value=0.01)
    Qmin  = st.number_input("Qmin (MVAR)", value=-300.0)
with c12:
    save_flag = st.button("💾 Save Condition")

# ======================================
# 5. Core Calculations  ← ALL FIXED
# ======================================

# Power angle
phi = np.arccos(np.clip(PF, 0, 1)) if PF > 0 else 0.0
if pf_type == "Leading":
    phi = -phi

# 3-phase reactive power (MW/MVAR)
Qr = Pr * np.tan(phi) if PF > 0 else 0.0

# ABCD parameters
Z_line, Y_line = None, None
if method == "Enter A,B,C directly":
    A = A_mag * np.exp(1j * np.radians(A_ang))
    B = B_mag * np.exp(1j * np.radians(B_ang))
    C = C_mag * np.exp(1j * np.radians(C_ang))
    D = A
elif method == "Enter Z,Y,length":
    Z_line = (R_ohm + 1j * X_ohm) * line_length
    Y_line = (G_val + 1j * B_sh)  * line_length
    theta  = np.sqrt(Z_line * Y_line)
    A = np.cosh(theta)
    B = Z_line * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Z_line
    C = Y_line * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Y_line
    D = A
else:
    w      = 2 * np.pi * 50
    Z_line = (r_ohm + 1j * w * L_val) * line_length
    Y_line = (1j * w * C_val)          * line_length
    theta  = np.sqrt(Z_line * Y_line)
    A = np.cosh(theta)
    B = Z_line * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Z_line
    C = Y_line * np.sinh(theta) / theta if np.abs(theta) > 1e-10 else Y_line
    D = A

# FIX 1 — Use PHASE voltage for current calculation
Vr_phase = Vr_mag / np.sqrt(3)          # kV (phase)
Vr_c     = Vr_phase + 0j                 # reference phasor

# FIX 2 — Per-phase power → per-phase current
Sr_ph = (Pr / 3) + 1j * (Qr / 3)        # MVA per phase
Ir    = (np.conj(Sr_ph) / np.conj(Vr_c)
         if np.abs(Vr_c) > 1e-10 else 0j)  # kA per phase

# Fault
if fault == "3-Phase Fault at Receiving End":
    Vr_c = fault_impedance * Ir

# Sending end (phase quantities)
Vs_ph = A * Vr_c + B * Ir               # kV phase
Is_ph = C * Vr_c + D * Ir               # kA phase

# Line voltages
Vs_mag = np.abs(Vs_ph) * np.sqrt(3)     # kV line-to-line

# 3-phase powers
Ss_3ph = 3 * Vs_ph * np.conj(Is_ph)
Ps     = np.real(Ss_3ph)                # MW
Qs     = np.imag(Ss_3ph)                # MVAR

Loss        = Ps - Pr
Efficiency  = (Pr / Ps * 100) if Ps > 1e-6 else 0.0
Voltage_reg = (Vs_mag - Vr_mag) / Vr_mag * 100

# Power angle δ (Vs leads Vr)
delta = np.degrees(np.angle(Vs_ph))

# Scalar ABCD magnitudes
Am = np.abs(A)
Bm = np.abs(B)
beta  = np.degrees(np.angle(B))    # β
alpha = np.degrees(np.angle(A))    # α

# ======================================================
# FIX 3 — CORRECT PCD GEOMETRY  (lecture formulas)
#
# Using LINE-TO-LINE voltages throughout.
#
# RPCD centre  n_r:
#   n_rx = -(|A|·Vr²/|B|)·cos(β−α)
#   n_ry = -(|A|·Vr²/|B|)·sin(β−α)
#
# SPCD centre  n_s:
#   n_sx = +(|A|·Vs²/|B|)·cos(β−α)
#   n_sy = -(|A|·Vs²/|B|)·sin(β−α)   ← NEGATIVE
#
# Radius (both circles):
#   R = |Vr|·|Vs| / |B|
#
# Pmax (rightmost point of RPCD):
#   Pmax = n_rx + R
# ======================================================
bma_rad = np.radians(beta - alpha)

AV2B_r = Am * Vr_mag**2 / Bm
n_rx   = -AV2B_r * np.cos(bma_rad)
n_ry   = -AV2B_r * np.sin(bma_rad)

AV2B_s = Am * Vs_mag**2 / Bm
n_sx   =  AV2B_s * np.cos(bma_rad)
n_sy   = -AV2B_s * np.sin(bma_rad)    # ← negative sign (lecture)

R_circle = Vr_mag * Vs_mag / Bm

Pmax_calc = n_rx + R_circle
Loading   = (Pr / Pmax_calc * 100) if Pmax_calc > 1e-6 else 0.0
Margin    = Pmax_calc - Pr
stability = "Stable" if Pr < Pmax_calc else "Unstable"

# ======================================
# 6. Save
# ======================================
if save_flag:
    st.session_state.history.append({
        "Mode": mode, "Pr (MW)": round(Pr, 2), "PF": round(PF, 3),
        "Ps (MW)": round(Ps, 2), "Qs (MVAR)": round(Qs, 2),
        "Loading %": round(Loading, 1), "Stability": stability
    })

# ======================================
# 7. KPIs
# ======================================
st.markdown("<hr>", unsafe_allow_html=True)
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Sending Power (Ps)", f"{Ps:.2f} MW")
m2.metric("Power Loss",         f"{Loss:.2f} MW")
m3.metric("Efficiency",         f"{Efficiency:.2f} %")
m4.metric("Line Loading",       f"{Loading:.2f} %")

stab_color = "#5cb85c" if stability == "Stable" else "#d9534f"
m5.markdown(f"""
    <div data-testid="metric-container" style="border-left:5px solid {stab_color};">
        <label style="color:#a0a5b5;font-size:15px;font-weight:600;">System Stability</label>
        <div style="color:{stab_color};font-weight:800;font-size:1.8rem;margin-top:-5px;">{stability}</div>
    </div>
""", unsafe_allow_html=True)

# ======================================
# 8. Tabs
# ======================================
st.markdown("<br>", unsafe_allow_html=True)

tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Power Circle Diagram",
    "📈 PV Curve",
    "⚡ Voltage Profile",
    "📋 Detailed Status",
    "💾 Saved Cases"
])

# ─────────────────────────────────────────────────────────
# TAB 1 — Power Circle Diagram
# ─────────────────────────────────────────────────────────
with tab1:

    col_plot, col_info = st.columns([3,1])

    with col_plot:

        # ==========================================
        # Circle Coordinates
        # ==========================================
        theta_arr = np.linspace(0, 2*np.pi, 800)

        Pr_circle = n_rx + R_circle*np.cos(theta_arr)
        Qr_circle = n_ry + R_circle*np.sin(theta_arr)

        Ps_circle = n_sx + R_circle*np.cos(theta_arr)
        Qs_circle = n_sy + R_circle*np.sin(theta_arr)

        # ==========================================
        # Figure
        # ==========================================
        fig, ax = plt.subplots(
            figsize=(10,10),
            facecolor='#0b1020'
        )

        ax.set_facecolor('#0b1020')

        # ==========================================
        # Style
        # ==========================================
        ax.tick_params(colors='white')

        for spine in ax.spines.values():
            spine.set_color('#555555')

        ax.grid(
            True,
            linestyle='--',
            linewidth=0.6,
            color='#7d7d7d',
            alpha=0.22
        )

        ax.set_title(
            "Combined Power Circle Diagram",
            fontsize=18,
            color='white',
            fontweight='bold',
            pad=20
        )

        ax.set_xlabel(
            "Active Power P (MW)",
            color='white',
            fontsize=13,
            fontweight='bold'
        )

        ax.set_ylabel(
            "Reactive Power Q (MVAR)",
            color='white',
            fontsize=13,
            fontweight='bold'
        )

        # ==========================================
        # Glow Effect
        # ==========================================
        ax.plot(
            Pr_circle,
            Qr_circle,
            color='#00d2ff',
            linewidth=8,
            alpha=0.05
        )

        ax.plot(
            Ps_circle,
            Qs_circle,
            color='#ff8d72',
            linewidth=8,
            alpha=0.05
        )

        # ==========================================
        # Main Circles
        # ==========================================
        ax.plot(
            Pr_circle,
            Qr_circle,
            color='#00d2ff',
            linewidth=2.8,
            label='RPCD (Receiving)',
            zorder=3
        )

        ax.plot(
            Ps_circle,
            Qs_circle,
            color='#ff8d72',
            linewidth=2.8,
            label='SPCD (Sending)',
            zorder=3
        )

        # ==========================================
        # Axes
        # ==========================================
        ax.axhline(
            0,
            color='white',
            linewidth=1.2,
            alpha=0.35
        )

        ax.axvline(
            0,
            color='white',
            linewidth=1.2,
            alpha=0.35
        )

        # ==========================================
        # Q Limits
        # ==========================================
        ax.axhline(
            Qmax,
            linestyle='--',
            linewidth=1.5,
            color='#d84fff',
            alpha=0.5
        )

        ax.axhline(
            Qmin,
            linestyle='--',
            linewidth=1.5,
            color='#d84fff',
            alpha=0.5
        )

        # ==========================================
        # Centers
        # ==========================================
        ax.scatter(
            n_rx,
            n_ry,
            color='#00d2ff',
            marker='x',
            s=130,
            linewidths=3,
            zorder=7
        )

        ax.scatter(
            n_sx,
            n_sy,
            color='#ff8d72',
            marker='x',
            s=130,
            linewidths=3,
            zorder=7
        )

        # ==========================================
        # Center Labels
        # ==========================================
        ax.text(
            n_rx + R_circle*0.03,
            n_ry + R_circle*0.05,
            "n_r",
            color='#00d2ff',
            fontsize=13,
            fontweight='bold'
        )

        ax.text(
            n_sx + R_circle*0.03,
            n_sy - R_circle*0.07,
            "n_s",
            color='#ff8d72',
            fontsize=13,
            fontweight='bold'
        )

        # ==========================================
        # Construction Lines
        # ==========================================
        ax.plot(
            [0, n_rx],
            [0, n_ry],
            linestyle='--',
            color='#00d2ff',
            linewidth=1.4,
            alpha=0.35
        )

        ax.plot(
            [0, n_sx],
            [0, n_sy],
            linestyle='--',
            color='#ff8d72',
            linewidth=1.4,
            alpha=0.35
        )

        # ==========================================
        # Operating Points
        # ==========================================
        ax.scatter(
            Pr,
            Qr,
            color='#00d2ff',
            edgecolors='white',
            s=160,
            linewidths=1.5,
            zorder=9
        )

        ax.scatter(
            Ps,
            Qs,
            color='#ff8d72',
            edgecolors='white',
            s=160,
            linewidths=1.5,
            zorder=9
        )

        # ==========================================
        # Labels
        # ==========================================
        ax.text(
            Pr + R_circle*0.025,
            Qr + R_circle*0.04,
            "Sr",
            color='#00d2ff',
            fontsize=15,
            fontweight='bold'
        )

        ax.text(
            Ps + R_circle*0.025,
            Qs - R_circle*0.05,
            "Ss",
            color='#ff8d72',
            fontsize=15,
            fontweight='bold'
        )

        # ==========================================
        # Radius Lines
        # ==========================================
        ax.plot(
            [n_rx, Pr],
            [n_ry, Qr],
            linestyle='--',
            linewidth=2,
            color='#00d2ff'
        )

        ax.plot(
            [n_sx, Ps],
            [n_sy, Qs],
            linestyle='--',
            linewidth=2,
            color='#ff8d72'
        )

        # ==========================================
        # Origin Vectors
        # ==========================================
        ax.annotate(
            "",
            xy=(Pr,Qr),
            xytext=(0,0),
            arrowprops=dict(
                arrowstyle='->',
                linewidth=2.8,
                color='#00d2ff'
            )
        )

        ax.annotate(
            "",
            xy=(Ps,Qs),
            xytext=(0,0),
            arrowprops=dict(
                arrowstyle='->',
                linewidth=2.8,
                color='#ff8d72'
            )
        )

        # ==========================================
        # Projection Line
        # ==========================================
        ax.plot(
            [Pr,Pr],
            [0,Qr],
            linestyle='--',
            linewidth=1,
            color='white',
            alpha=0.25
        )

        # ==========================================
        # Pmax
        # ==========================================
        ax.plot(
            [n_rx,Pmax_calc],
            [n_ry,0],
            linestyle=':',
            linewidth=2.2,
            color='#ff4fd8'
        )

        ax.scatter(
            Pmax_calc,
            0,
            color='#ff4fd8',
            marker='^',
            s=220,
            zorder=10
        )

        ax.text(
            Pmax_calc + R_circle*0.03,
            R_circle*0.04,
            f'Pmax = {Pmax_calc:.0f} MW',
            color='#ff4fd8',
            fontsize=14,
            fontweight='bold'
        )

        # ==========================================
        # PF Arc
        # ==========================================
        phi_deg_plot = np.degrees(phi)

        arc_r = R_circle * 0.18

        arc_theta = np.linspace(
            0,
            np.radians(phi_deg_plot),
            80
        )

        ax.plot(
            arc_r*np.cos(arc_theta),
            arc_r*np.sin(arc_theta),
            color='#00d2ff',
            linewidth=1.5,
            linestyle=':'
        )

        ax.text(
            arc_r*1.05,
            arc_r*0.25,
            f'φr = {abs(phi_deg_plot):.1f}°',
            color='#00d2ff',
            fontsize=12
        )

        # ==========================================
        # Delta
        # ==========================================
        delta = np.degrees(
            np.angle(Vs_ph) - np.angle(Vr_c)
        )

        ax.text(
            R_circle*0.02,
            -R_circle*0.09,
            f'δ = {delta:.1f}°',
            fontsize=12,
            color='white',
            bbox=dict(
                facecolor='#1a1f33',
                edgecolor='#666666',
                boxstyle='round,pad=0.3',
                alpha=0.9
            )
        )

        # ==========================================
        # Beta Alpha
        # ==========================================
        ax.text(
            n_rx + R_circle*0.04,
            n_ry - R_circle*0.10,
            f'β−α = {beta-alpha:.1f}°',
            color='#bbbbbb',
            fontsize=11
        )

        # ==========================================
        # Limits
        # ==========================================
        all_x = [
            Pr,
            Ps,
            n_rx + R_circle,
            n_rx - R_circle,
            n_sx + R_circle,
            n_sx - R_circle,
            Pmax_calc,
            0
        ]

        all_y = [
            Qr,
            Qs,
            n_ry + R_circle,
            n_ry - R_circle,
            n_sy + R_circle,
            n_sy - R_circle,
            0
        ]

        x_min = min(all_x)
        x_max = max(all_x)

        y_min = min(all_y)
        y_max = max(all_y)

        pad_x = (x_max - x_min) * 0.15
        pad_y = (y_max - y_min) * 0.15

        ax.set_xlim(
            x_min - pad_x,
            x_max + pad_x*1.6
        )

        ax.set_ylim(
            y_min - pad_y,
            y_max + pad_y
        )

        ax.set_aspect('equal')

        # ==========================================
        # Legend
        # ==========================================
        leg = ax.legend(
            facecolor='#1a1f33',
            edgecolor='#555555',
            fontsize=11,
            loc='lower right'
        )

        for txt in leg.get_texts():
            txt.set_color('white')

        st.pyplot(fig)

    # =====================================================
    # SIDE INFO
    # =====================================================
    with col_info:

        st.markdown("#### Point Details")

        st.info(f"**Pr:** {Pr:.2f} MW")
        st.info(f"**Qr:** {Qr:.2f} MVAR")

        st.info(f"**Ps:** {Ps:.2f} MW")
        st.info(f"**Qs:** {Qs:.2f} MVAR")

        st.info(f"**Pmax:** {Pmax_calc:.2f} MW")

        st.info(f"**Margin:** {Margin:.2f} MW")

        st.info(f"**δ:** {delta:.2f}°")

        st.info(f"**β−α:** {beta-alpha:.2f}°")

        st.info(f"**Radius:** {R_circle:.2f} MVA")
# ─────────────────────────────────────────────────────────
# TAB 2 — PV Curve
# ─────────────────────────────────────────────────────────
with tab2:
    st.markdown("#### PV Curve (Voltage Stability)")
    P_range = np.linspace(0, Pmax_calc * 1.2, 60)
    V_list  = []
    for P_t in P_range:
        Q_t    = P_t * np.tan(phi) if PF > 0 else 0
        Sr_t   = (P_t / 3) + 1j * (Q_t / 3)
        Ir_t   = (np.conj(Sr_t) / np.conj(Vr_c)
                  if np.abs(Vr_c) > 1e-10 else 0j)
        Vs_t   = A * Vr_c + B * Ir_t
        V_list.append(np.abs(Vs_t) * np.sqrt(3))

    fig_pv, ax_pv = plt.subplots(figsize=(10, 4))
    style_dark_plot(ax_pv, fig_pv, "PV Curve",
                    "Load P (MW)", "Sending Voltage |Vs| L-L (kV)")
    ax_pv.plot(P_range, V_list, color='#e14eca', lw=2)
    ax_pv.axvline(Pmax_calc, ls='--', color='#ff8d72',
                  label=f"Pmax ({Pmax_calc:.1f} MW)")
    ax_pv.scatter(Pr, Vs_mag, color='#00d2ff', s=100,
                  zorder=5, label="Operating point")
    leg_pv = ax_pv.legend(facecolor='#1e2233', edgecolor='#444a5e')
    fix_legend_color(leg_pv)
    st.pyplot(fig_pv)

# ─────────────────────────────────────────────────────────
# TAB 3 — Voltage Profile
# ─────────────────────────────────────────────────────────
with tab3:
    st.markdown("#### Voltage Profile along the Line")
    if line_length > 0 and Z_line is not None and Y_line is not None:
        Z_pk  = Z_line / line_length
        Y_pk  = Y_line / line_length
        gamma = np.sqrt(Z_pk * Y_pk)
        Zc    = (np.sqrt(Z_pk / Y_pk)
                 if np.abs(Y_pk) > 1e-15 else 1.0)
        dist  = np.linspace(0, line_length, 200)
        V_prof = [np.abs(Vr_c * np.cosh(gamma * x)
                         + Ir * Zc * np.sinh(gamma * x)) * np.sqrt(3)
                  for x in dist]
        fig2, ax2 = plt.subplots(figsize=(10, 4))
        style_dark_plot(ax2, fig2,
                        "Voltage Profile along the Transmission Line",
                        "Distance from Receiving End (km)",
                        "Voltage (kV) L-L")
        ax2.plot(dist, V_prof, color='#00d2ff', lw=2)
        ax2.invert_xaxis()
        st.pyplot(fig2)
    else:
        st.info("Voltage profile needs Z & Y. Use 'Enter Z,Y,length' or 'Enter r,L,C,length'.")

# ─────────────────────────────────────────────────────────
# TAB 4 — Detailed Status
# ─────────────────────────────────────────────────────────
with tab4:
    st.markdown("#### 📋 System Detailed Status")
    cs1, cs2 = st.columns(2)
    with cs1:
        st.markdown("**Voltages & Currents (per-phase)**")
        st.code(
            f"Vr (phase) = {polar_str(Vr_c,  'kV')}\n"
            f"Ir (phase) = {polar_str(Ir,     'kA')}\n"
            f"Vs (phase) = {polar_str(Vs_ph,  'kV')}\n"
            f"Is (phase) = {polar_str(Is_ph,  'kA')}\n"
            f"Vs (L-L)   = {Vs_mag:.4f} kV\n"
            f"Vr (L-L)   = {Vr_mag:.4f} kV"
        )
        st.markdown("**Performance**")
        r_stat = ('Overexcited'  if Qs > Qmax
                  else 'Underexcited' if Qs < Qmin
                  else 'Within Limits')
        st.code(
            f"Power angle δ  = {delta:.4f}°\n"
            f"Voltage Reg    = {Voltage_reg:.4f} %\n"
            f"Efficiency     = {Efficiency:.4f} %\n"
            f"Reactive Status= {r_stat}"
        )
    with cs2:
        st.markdown("**ABCD Parameters**")
        st.code(
            f"A = {polar_str(A)}\n"
            f"B = {polar_str(B, 'Ω')}\n"
            f"C = {polar_str(C, 'S')}\n"
            f"D = {polar_str(D)}"
        )
        st.markdown("**PCD Key Values**")
        st.code(
            f"n_r  = ({n_rx:.2f} , {n_ry:.2f})  MW / MVAR\n"
            f"n_s  = ({n_sx:.2f} , {n_sy:.2f})  MW / MVAR\n"
            f"R    = {R_circle:.2f}  MVA  (both circles)\n"
            f"Pmax = {Pmax_calc:.2f}  MW\n"
            f"β−α  = {beta - alpha:.2f}°\n"
            f"|A|Vr²/|B| = {AV2B_r:.2f}  MVA\n"
            f"|A|Vs²/|B| = {AV2B_s:.2f}  MVA"
        )

# ─────────────────────────────────────────────────────────
# TAB 5 — Saved Cases
# ─────────────────────────────────────────────────────────
with tab5:
    st.markdown("#### 💾 Saved Cases")
    if st.session_state.history:
        hdrs = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
        for h, lbl in zip(hdrs, ["Mode","Pr (MW)","PF","Ps (MW)",
                                   "Qs (MVAR)","Load %","Stability","Del"]):
            h.markdown(f"**{lbl}**")
        st.markdown("<hr style='margin:5px 0'>", unsafe_allow_html=True)
        for i, case in enumerate(st.session_state.history):
            cols = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
            cols[0].write(case["Mode"])
            cols[1].write(case["Pr (MW)"])
            cols[2].write(case["PF"])
            cols[3].write(case["Ps (MW)"])
            cols[4].write(case["Qs (MVAR)"])
            cols[5].write(case["Loading %"])
            sc = "#5cb85c" if case["Stability"] == "Stable" else "#d9534f"
            cols[6].markdown(
                f"<span style='color:{sc};font-weight:bold;'>{case['Stability']}</span>",
                unsafe_allow_html=True
            )
            if cols[7].button("🗑️", key=f"del_{i}"):
                st.session_state.history.pop(i)
                st.rerun()
    else:
        st.info("No cases saved yet.")
