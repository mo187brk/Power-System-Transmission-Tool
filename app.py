import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as pe

# ======================================
# 1. Page Config & Styling
# ======================================
st.set_page_config(page_title="Power Transmission Analysis", layout="wide", page_icon="⚡")

st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@400;700&family=Share+Tech+Mono&display=swap');

    .stApp { background: linear-gradient(135deg, #0d1117 0%, #0a0e1a 100%); color: #e6edf3; }

    h1, h2, h3, h4 {
        color: #58a6ff !important;
        font-family: 'Orbitron', monospace !important;
        letter-spacing: 1px;
    }

    div[data-testid="stWidgetLabel"] p, label p,
    .stSlider label p, .stNumberInput label p, .stSelectbox label p {
        color: #8b949e !important; font-weight: 600 !important; font-size: 13px !important;
        font-family: 'Share Tech Mono', monospace !important;
        text-transform: uppercase; letter-spacing: 0.5px;
    }

    div[data-testid="metric-container"] {
        background: linear-gradient(135deg, #161b22 0%, #0d1117 100%);
        padding: 15px; border-radius: 10px;
        box-shadow: 0 0 15px rgba(88,166,255,0.1);
        border: 1px solid #21262d;
        transition: all 0.3s ease;
    }
    div[data-testid="metric-container"]:hover {
        transform: translateY(-3px);
        box-shadow: 0 0 25px rgba(88,166,255,0.25);
        border-color: #58a6ff;
    }

    div[data-testid="stMetricLabel"] p { color: #8b949e !important; font-size: 12px !important; font-weight: 600 !important; text-transform: uppercase; letter-spacing: 1px; }
    div[data-testid="stMetricValue"] { color: #58a6ff !important; font-weight: 800 !important; font-size: 1.6rem !important; font-family: 'Share Tech Mono', monospace !important; }

    div[data-testid="stButton"] button {
        background: linear-gradient(90deg, #1f6feb 0%, #388bfd 100%) !important;
        color: #ffffff !important; border-radius: 8px !important; border: none !important;
        font-weight: 700 !important; font-size: 14px !important; height: 44px !important;
        width: 100% !important; margin-top: 27px !important;
        box-shadow: 0 0 15px rgba(56,139,253,0.3) !important;
        font-family: 'Share Tech Mono', monospace !important;
        letter-spacing: 1px !important; transition: all 0.3s ease !important;
    }
    div[data-testid="stButton"] button:hover {
        box-shadow: 0 0 30px rgba(56,139,253,0.6) !important;
        transform: translateY(-2px) !important;
    }

    .stTabs [data-baseweb="tab"] {
        font-family: 'Share Tech Mono', monospace !important;
        color: #8b949e !important; font-size: 13px !important;
        text-transform: uppercase; letter-spacing: 1px;
    }
    .stTabs [aria-selected="true"] { color: #58a6ff !important; }

    hr { border-color: #21262d; margin: 15px 0; }

    code { background: #161b22 !important; color: #79c0ff !important;
           font-family: 'Share Tech Mono', monospace !important; }
    .stCode { border: 1px solid #21262d !important; border-radius: 8px !important; }
    </style>
""", unsafe_allow_html=True)

# ======================================
# 2. Helper Functions
# ======================================
def polar_str(value, unit=""):
    mag = np.abs(value)
    ang = np.degrees(np.angle(value))
    return f"{mag:.4f} ∠ {ang:.2f}° {unit}"

def style_dark_plot(ax, fig, title="", xlabel="", ylabel=""):
    """Apply consistent dark theme to matplotlib plots."""
    fig.patch.set_facecolor('#0d1117')
    ax.set_facecolor('#0d1117')
    ax.tick_params(colors='#8b949e', labelsize=10)
    for spine in ax.spines.values():
        spine.set_color('#21262d')
        spine.set_linewidth(0.8)
    ax.grid(True, linestyle=':', color='#21262d', linewidth=0.6, alpha=0.8)
    if title:
        ax.set_title(title, fontsize=13, fontweight='bold', color='#58a6ff', pad=12,
                     fontfamily='monospace')
    if xlabel:
        ax.set_xlabel(xlabel, color='#8b949e', fontweight='bold', fontsize=11)
    if ylabel:
        ax.set_ylabel(ylabel, color='#8b949e', fontweight='bold', fontsize=11)

def fix_legend(legend):
    if legend:
        legend.get_frame().set_facecolor('#161b22')
        legend.get_frame().set_edgecolor('#21262d')
        for text in legend.get_texts():
            text.set_color('#c9d1d9')
            text.set_fontsize(10)

# ======================================
# 3. Session State
# ======================================
if 'history' not in st.session_state:
    st.session_state.history = []

# ======================================
# 4. Header
# ======================================
st.title("⚡ Power Transmission Line Analysis")
st.markdown("<hr>", unsafe_allow_html=True)

# ======================================
# 5. Input Panel
# ======================================
st.markdown("<span style='color:#58a6ff; font-weight:bold; font-family:monospace; font-size:13px; text-transform:uppercase; letter-spacing:1px;'>► 01 | OPERATING CONDITIONS</span>", unsafe_allow_html=True)

c1, c2, c3, c4 = st.columns(4)
with c1:
    mode = st.selectbox("Operating Mode", ["Load", "No Load"])
with c2:
    Pr = st.number_input("Receiving Power Pr (MW)", value=100.0, min_value=0.0) if mode == "Load" else 0.0
with c3:
    PF = st.slider("Power Factor", 0.0, 1.0, 0.9) if mode == "Load" else 0.0
with c4:
    pf_type = st.selectbox("PF Type", ["Lagging", "Leading"]) if mode == "Load" else "Lagging"

st.markdown("<br><span style='color:#58a6ff; font-weight:bold; font-family:monospace; font-size:13px; text-transform:uppercase; letter-spacing:1px;'>► 02 | LINE PARAMETERS</span>", unsafe_allow_html=True)

c5, c6, c7, c8 = st.columns(4)
with c5:
    method = st.selectbox("ABCD Input Method", ["Enter A,B,C directly", "Enter Z,Y,length", "Enter r,L,C,length"])
    line_length = st.number_input("Line Length (km) [For Profile]", value=100.0)
with c6:
    if method == "Enter A,B,C directly":
        A_mag = st.number_input("|A|", value=0.986, format="%.4f")
        A_ang = st.number_input("Angle A (°)", value=0.1)
    elif method == "Enter Z,Y,length":
        R_km = st.number_input("R (Ω/km)", value=0.1)
        X_km = st.number_input("X (Ω/km)", value=0.4)
    else:
        r_km = st.number_input("r (Ω/km)", value=0.1)
        L_km = st.number_input("L (H/km)", value=1e-3, format="%.5f")
with c7:
    if method == "Enter A,B,C directly":
        B_mag = st.number_input("|B| (Ω)", value=60.0)
        B_ang = st.number_input("Angle B (°)", value=83.0)
    elif method == "Enter Z,Y,length":
        G_km = st.number_input("G (S/km)", value=0.0)
        B_sh_km = st.number_input("B (S/km)", value=0.000001, format="%.6f")
    else:
        C_km = st.number_input("C (F/km)", value=1e-8, format="%.9f")
with c8:
    if method == "Enter A,B,C directly":
        C_mag = st.number_input("|C| (S)", value=0.0004, format="%.5f")
        C_ang = st.number_input("Angle C (°)", value=90.0)

st.markdown("<br><span style='color:#58a6ff; font-weight:bold; font-family:monospace; font-size:13px; text-transform:uppercase; letter-spacing:1px;'>► 03 | CONSTRAINTS & ACTIONS</span>", unsafe_allow_html=True)

c9, c10, c11, c12 = st.columns(4)
with c9:
    Vr_mag = st.number_input("Receiving Voltage Vr (kV)", value=220.0)
with c10:
    fault = st.selectbox("Fault Simulation", ["No Fault", "3-Phase Fault at Receiving End"])
    Qmax = st.number_input("Qmax (MVAR)", value=300.0)
with c11:
    fault_impedance = st.number_input("Fault Impedance (Ω)", value=0.01)
    Qmin = st.number_input("Qmin (MVAR)", value=-300.0)
with c12:
    save_flag = st.button("💾 Save Condition")

# ======================================
# 6. ABCD Parameter Calculation
# ======================================
Z_line = None
Y_line = None

if method == "Enter A,B,C directly":
    A = A_mag * np.exp(1j * np.radians(A_ang))
    B = B_mag * np.exp(1j * np.radians(B_ang))
    C = C_mag * np.exp(1j * np.radians(C_ang))
    D = A  # Symmetric network: D = A

elif method == "Enter Z,Y,length":
    Z_line = (R_km + 1j * X_km) * line_length
    Y_line = (G_km + 1j * B_sh_km) * line_length
    gamma = np.sqrt(Z_line * Y_line)
    # Avoid division by zero
    sinh_g = np.sinh(gamma)
    if abs(gamma) > 1e-10:
        A = np.cosh(gamma)
        B = Z_line * sinh_g / gamma
        C = Y_line * sinh_g / gamma
    else:
        A = 1.0 + 0j
        B = Z_line
        C = Y_line
    D = A

else:  # r, L, C per km
    w = 2 * np.pi * 50
    Z_line = (r_km + 1j * w * L_km) * line_length
    Y_line = (1j * w * C_km) * line_length
    gamma = np.sqrt(Z_line * Y_line)
    sinh_g = np.sinh(gamma)
    if abs(gamma) > 1e-10:
        A = np.cosh(gamma)
        B = Z_line * sinh_g / gamma
        C = Y_line * sinh_g / gamma
    else:
        A = 1.0 + 0j
        B = Z_line
        C = Y_line
    D = A

# ======================================
# 7. Power & Phasor Calculations
# ======================================
phi = np.arccos(np.clip(PF, 0, 1)) if PF != 0 else 0
if pf_type == "Leading":
    phi = -phi

# Per-phase quantities (Vr is reference phasor)
Vr = Vr_mag + 0j  # kV, reference

Pr_phase = Pr / 3      # MW per phase
Qr = Pr * np.tan(phi) if PF != 0 else 0.0
Qr_phase = Qr / 3     # MVAR per phase
Sr_phase = Pr_phase + 1j * Qr_phase  # MVA per phase

# Receiving end current
Ir = np.conj(Sr_phase) / np.conj(Vr) if abs(Vr) > 1e-10 else 0.0

# Apply fault if selected
if fault == "3-Phase Fault at Receiving End":
    Vr = fault_impedance * Ir  # Fault voltage

# Sending end phasors
Vs = A * Vr + B * Ir
Is = C * Vr + D * Ir

# Three-phase power
Ss = 3 * Vs * np.conj(Is)
Ps = np.real(Ss)
Qs = np.imag(Ss)

Loss_MW = Ps - Pr
Efficiency = (Pr / Ps) * 100 if abs(Ps) > 1e-10 else 0.0
Voltage_reg = ((np.abs(Vs) - Vr_mag) / Vr_mag) * 100 if Vr_mag > 0 else 0.0

# ======================================
# 8. CORRECTED Power Circle Parameters
# ======================================
# From textbook derivation:
# Receiving-end circle:
#   Center n:  Nrx = -(|A|/|B|)*Vr² * cos(β-α)
#              Nry = -(|A|/|B|)*Vr² * sin(β-α)
#   Radius Rr = Vr * |Vs| / |B|   (depends on Vs)
#
# Sending-end circle:
#   Center n': Nsx = -(|D|/|B|)*Vs² * cos(β-δ_D)  where D=A for symmetric
#              Nsy = -(|D|/|B|)*Vs² * sin(β-δ_D)
#   Radius Rs = Vr * |Vs| / |B|   (same radius)
#
# Pmax (receiving) = Rr - Nrx  = VrVs/|B| + (A/B)*Vr²*cos(β-α)

A_mag_val = np.abs(A)
A_ang_val = np.angle(A)   # α in radians
B_mag_val = np.abs(B)
B_ang_val = np.angle(B)   # β in radians
Vs_mag    = np.abs(Vs)

beta_alpha = B_ang_val - A_ang_val  # β - α

# Receiving-end circle center (CORRECTED — negative sign from derivation)
Nrx = -(A_mag_val / B_mag_val) * Vr_mag**2 * np.cos(beta_alpha)
Nry = -(A_mag_val / B_mag_val) * Vr_mag**2 * np.sin(beta_alpha)

# Sending-end circle center (D = A for symmetric network)
D_mag_val = np.abs(D)
D_ang_val = np.angle(D)   # δ_D = α
beta_delta = B_ang_val - D_ang_val  # β - δ_D = β - α (same for symmetric)
Nsx = -(D_mag_val / B_mag_val) * Vs_mag**2 * np.cos(beta_delta)
Nsy = -(D_mag_val / B_mag_val) * Vs_mag**2 * np.sin(beta_delta)

# Common radius (both circles share same radius VrVs/B)
R_circle = (Vr_mag * Vs_mag) / B_mag_val

# CORRECTED Pmax: including the center offset
Pmax_calc = R_circle - Nrx   # = VrVs/B + (A/B)*Vr²*cos(β-α)

# Stability
Loading  = (Pr / Pmax_calc) * 100 if Pmax_calc > 1e-10 else 0.0
Margin   = Pmax_calc - Pr
stability = "Stable" if Pr < Pmax_calc else "Unstable"

# ======================================
# 9. Save Logic
# ======================================
if save_flag:
    st.session_state.history.append({
        "Mode": mode, "Pr (MW)": round(Pr, 2), "PF": round(PF, 3),
        "Ps (MW)": round(Ps, 2), "Qs (MVAR)": round(Qs, 2),
        "Qr (MVAR)": round(Qr, 2),
        "Loading %": round(Loading, 1), "Stability": stability
    })

# ======================================
# 10. KPI Metrics
# ======================================
st.markdown("<hr>", unsafe_allow_html=True)
m1, m2, m3, m4, m5, m6 = st.columns(6)
m1.metric("Ps (MW)", f"{Ps:.2f}")
m2.metric("Qs (MVAR)", f"{Qs:.2f}")
m3.metric("Line Loss (MW)", f"{Loss_MW:.2f}")
m4.metric("Efficiency (%)", f"{Efficiency:.2f}")
m5.metric("Line Loading (%)", f"{Loading:.2f}")

stability_color = "#3fb950" if stability == "Stable" else "#f85149"
m6.markdown(f"""
    <div data-testid="metric-container" style="border-left: 4px solid {stability_color}; padding:15px; border-radius:10px; background:linear-gradient(135deg,#161b22,#0d1117); border-top:1px solid #21262d; border-right:1px solid #21262d; border-bottom:1px solid #21262d;">
        <p style="color:#8b949e; font-size:12px; font-weight:600; margin:0 0 4px 0; text-transform:uppercase; letter-spacing:1px; font-family:monospace;">STABILITY</p>
        <p style="color:{stability_color}; font-weight:800; font-size:1.4rem; margin:0; font-family:monospace;">{stability}</p>
    </div>
""", unsafe_allow_html=True)

# ======================================
# 11. Tabs
# ======================================
st.markdown("<br>", unsafe_allow_html=True)
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊  Power Circle Diagram",
    "📈  PV Curve",
    "⚡  Voltage Profile",
    "📋  System Status",
    "💾  Saved Cases"
])

# ============================================================
# TAB 1 — CORRECTED Power Circle Diagram
# ============================================================
with tab1:
    col_plot, col_info = st.columns([3, 1])

    with col_plot:
        fig, ax = plt.subplots(figsize=(9, 9))
        style_dark_plot(ax, fig,
                        "Combined Power Circle Diagram",
                        "Active Power   P  (MW)",
                        "Reactive Power   Q  (MVAR)")

        # --- Plot full circles ---
        theta = np.linspace(0, 2 * np.pi, 500)

        # Receiving-end circle (teal)
        Pr_circle = Nrx + R_circle * np.cos(theta)
        Qr_circle = Nry + R_circle * np.sin(theta)
        ax.plot(Pr_circle, Qr_circle, color='#1abc9c', linewidth=2.0,
                label='Receiving-end circle', zorder=3)

        # Sending-end circle (blue, dashed)
        Ps_circle = Nsx + R_circle * np.cos(theta)
        Qs_circle = Nsy + R_circle * np.sin(theta)
        ax.plot(Ps_circle, Qs_circle, color='#3b82f6', linewidth=2.0,
                linestyle='--', dashes=(6, 3), label='Sending-end circle', zorder=3)

        # Qmax / Qmin lines
        ax.axhline(Qmax, linestyle='--', color='#f59e0b', linewidth=1.2,
                   alpha=0.7, label=f'Qmax = {Qmax:.0f} MVAR')
        ax.axhline(Qmin, linestyle='--', color='#f59e0b', linewidth=1.2,
                   alpha=0.7, label=f'Qmin = {Qmin:.0f} MVAR')

        # --- Center points ---
        ax.plot(Nrx, Nry, 'o', color='#1abc9c', markersize=8, zorder=5)
        ax.annotate('n', (Nrx, Nry),
                    textcoords="offset points", xytext=(-18, 5),
                    fontsize=12, color='#1abc9c', fontweight='bold',
                    fontfamily='monospace')

        ax.plot(Nsx, Nsy, 's', color='#3b82f6', markersize=8, zorder=5)
        ax.annotate("n'", (Nsx, Nsy),
                    textcoords="offset points", xytext=(8, 5),
                    fontsize=12, color='#3b82f6', fontweight='bold',
                    fontfamily='monospace')

        # --- n → n' dashed connector ---
        ax.plot([Nrx, Nsx], [Nry, Nsy], '--', color='#6b7280',
                linewidth=0.8, alpha=0.6, zorder=2)

        # --- CORRECTED arrows: from center to operating point K ---
        # Sr vector: from n to K=(Pr, Qr)
        ax.annotate("", xy=(Pr, Qr), xytext=(Nrx, Nry),
                    arrowprops=dict(arrowstyle='->', color='#1abc9c',
                                   lw=2.2, mutation_scale=18))

        # Ss vector: from n' to K=(Pr, Qr)
        ax.annotate("", xy=(Pr, Qr), xytext=(Nsx, Nsy),
                    arrowprops=dict(arrowstyle='->', color='#3b82f6',
                                   lw=2.2, mutation_scale=18))

        # --- VrVs/B radius arm (n to Pmax on horizontal through n) ---
        ax.annotate("", xy=(Pmax_calc, Nry), xytext=(Nrx, Nry),
                    arrowprops=dict(arrowstyle='->', color='#9ca3af',
                                   lw=1.4, mutation_scale=14))
        ax.text((Nrx + Pmax_calc) / 2, Nry + 6,
                "VrVs/|B|", ha='center', va='bottom',
                fontsize=9, color='#9ca3af', fontfamily='monospace')

        # --- Operating point K ---
        ax.plot(Pr, Qr, 'o', color='#f87171', markersize=11,
                zorder=7, label=f'Operating point K')
        ax.annotate(f'  K\n  Sr = {Pr:.0f}+j{Qr:.0f}',
                    (Pr, Qr), textcoords="offset points",
                    xytext=(12, 6), fontsize=10, color='#f87171',
                    fontfamily='monospace')

        # --- Pmax marker on P-axis ---
        ax.plot(Pmax_calc, Nry, 'D', color='#a855f7',
                markersize=9, zorder=6, label=f'Pmax = {Pmax_calc:.0f} MW')
        ax.annotate(f'  Pm\n  {Pmax_calc:.0f} MW',
                    (Pmax_calc, Nry), textcoords="offset points",
                    xytext=(8, -22), fontsize=10, color='#a855f7',
                    fontfamily='monospace')

        # --- P₀ = projection of n onto P-axis ---
        ax.plot(Nrx, 0, 'v', color='#1abc9c', markersize=7, zorder=5)
        ax.annotate(f'P₀\n{Nrx:.0f}',
                    (Nrx, 0), textcoords="offset points",
                    xytext=(-10, -28), fontsize=9, color='#1abc9c',
                    fontfamily='monospace')

        # Dashed vertical from n to P-axis
        ax.plot([Nrx, Nrx], [Nry, 0], '--', color='#1abc9c',
                linewidth=0.8, alpha=0.5, zorder=2)

        # --- Pr marker on P-axis ---
        ax.plot(Pr, 0, '^', color='#f87171', markersize=8, zorder=6)
        ax.annotate(f'Pr\n{Pr:.0f}',
                    (Pr, 0), textcoords="offset points",
                    xytext=(4, -26), fontsize=9, color='#f87171',
                    fontfamily='monospace')
        # Dashed vertical from K to P-axis
        ax.plot([Pr, Pr], [Qr, 0], '--', color='#f87171',
                linewidth=0.8, alpha=0.5, zorder=2)

        # --- Qr marker on Q-axis ---
        ax.plot(0, Qr, '<', color='#f87171', markersize=7, zorder=6)
        ax.annotate(f'Qr\n{Qr:.0f}',
                    (0, Qr), textcoords="offset points",
                    xytext=(-52, 0), fontsize=9, color='#f87171',
                    fontfamily='monospace')
        # Dashed horizontal from K to Q-axis
        ax.plot([0, Pr], [Qr, Qr], '--', color='#f87171',
                linewidth=0.8, alpha=0.5, zorder=2)

        # --- Power angle δ annotation (at center n) ---
        delta_ang = np.degrees(np.angle(Vs) - np.angle(Vr))
        # φr angle (from P-axis to Sr vector, at n)
        phi_r_deg = np.degrees(np.arctan2(Qr - Nry, Pr - Nrx))
        ax.text(Nrx + 20, Nry + 8,
                f"φᵣ = {abs(phi_r_deg):.1f}°",
                fontsize=9, color='#e2e8f0', fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#161b22',
                          edgecolor='#374151', alpha=0.85))

        # β - α angle annotation
        bma_deg = np.degrees(beta_alpha)
        ax.text(Nrx - 15, Nry - 20,
                f"β−α = {bma_deg:.1f}°",
                fontsize=9, color='#60a5fa', fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#161b22',
                          edgecolor='#374151', alpha=0.85))

        # --- Current operating Ps point on sending circle ---
        ax.plot(Ps, Qs, 's', color='#3b82f6', markersize=9,
                zorder=7, label=f"Sending point K'")
        ax.annotate(f"  K'\n  Ss={Ps:.0f}+j{Qs:.0f}",
                    (Ps, Qs), textcoords="offset points",
                    xytext=(10, 6), fontsize=9, color='#3b82f6',
                    fontfamily='monospace')

        # Power loss line K → K'
        ax.plot([Pr, Ps], [Qr, Qs], '-', color='#f59e0b',
                linewidth=1.5, label='Power loss (K→K\')', zorder=4)

        # --- Axes through origin ---
        ax.axhline(0, color='#374151', linewidth=1.0, zorder=1)
        ax.axvline(0, color='#374151', linewidth=1.0, zorder=1)

        # --- Auto window ---
        all_x = [Pr_circle.min(), Pr_circle.max(),
                 Ps_circle.min(), Ps_circle.max(), Pr, Ps, Pmax_calc, 0, Nrx, Nsx]
        all_y = [Qr_circle.min(), Qr_circle.max(),
                 Qs_circle.min(), Qs_circle.max(), Qr, Qs, 0, Nry, Nsy]
        xpad = (max(all_x) - min(all_x)) * 0.18
        ypad = (max(all_y) - min(all_y)) * 0.18
        ax.set_xlim(min(all_x) - xpad, max(all_x) + xpad)
        ax.set_ylim(min(all_y) - ypad, max(all_y) + ypad)
        ax.set_aspect('equal')

        leg = ax.legend(loc='upper left', fontsize=9,
                        facecolor='#161b22', edgecolor='#21262d')
        fix_legend(leg)
        st.pyplot(fig)
        plt.close(fig)

    with col_info:
        st.markdown("#### 📌 Circle Data")
        st.markdown(f"""
        <div style='font-family:monospace; font-size:12px; line-height:2;
                    background:#161b22; padding:14px; border-radius:8px;
                    border:1px solid #21262d; color:#c9d1d9;'>
        <b style='color:#1abc9c'>Receiving Circle</b><br>
        Center n:<br>
        &nbsp;Nrx = {Nrx:.2f} MW<br>
        &nbsp;Nry = {Nry:.2f} MVAR<br>
        Radius = {R_circle:.2f} MVA<br><br>
        <b style='color:#3b82f6'>Sending Circle</b><br>
        Center n':<br>
        &nbsp;Nsx = {Nsx:.2f} MW<br>
        &nbsp;Nsy = {Nsy:.2f} MVAR<br>
        Radius = {R_circle:.2f} MVA<br><br>
        <b style='color:#a855f7'>Power Limits</b><br>
        Pmax = {Pmax_calc:.2f} MW<br>
        Margin = {Margin:.2f} MW<br>
        Loading = {Loading:.1f}%<br><br>
        <b style='color:#f87171'>Operating Point K</b><br>
        Pr = {Pr:.2f} MW<br>
        Qr = {Qr:.2f} MVAR<br>
        Ps = {Ps:.2f} MW<br>
        Qs = {Qs:.2f} MVAR
        </div>
        """, unsafe_allow_html=True)

# ============================================================
# TAB 2 — PV Curve (Corrected: fixed Vs, vary Pr, solve Vr)
# ============================================================
with tab2:
    st.markdown("#### 📈 PV Curve — Voltage Stability (Fixed Vs, Varying Load)")

    P_range = np.linspace(0, Pmax_calc * 1.1, 200)
    Vr_upper = []
    Vr_lower = []

    Vs_fixed = np.abs(Vs)

    for P_test in P_range:
        Q_test = P_test * np.tan(phi) if PF > 1e-6 else 0.0
        # From: Vr² - Vr*(Vs_fixed/A_mag_val)*cos(...) + ... = 0  (simplified)
        # Use per-phase, iterative: start from rated Vr
        # Or use circle geometry: Vr² = (A/B)² * Vr⁴ + ... 
        # Practical approach: solve |Vs|² = |A*Vr + B*Ir|² with Ir from load
        # Use two-point quadratic: Vr_high and Vr_low from quadratic
        # Sr_test / Vr = Ir* → |Vs|² = |A|²Vr² + 2Re(A*B*Sr_test*) + |B|²|Sr_test|²/Vr²
        # Let u = Vr², Sr_mag = |(P+jQ)/3|
        Sr_test_phase = (P_test + 1j * Q_test) / 3.0
        Sr_mag = np.abs(Sr_test_phase)
        A2 = A_mag_val**2
        B2 = B_mag_val**2
        # |Vs|² = A²u + 2Re(A*B*(P-jQ)/3) + B²(P²+Q²)/(9u)
        # → A²u² + [2Re(AB̄ S*) - |Vs|²]u + B²|S|²/9 = 0  (×u)
        AB_conj = A * np.conj(B)
        coeff_linear = 2 * np.real(AB_conj * np.conj(Sr_test_phase)) - Vs_fixed**2
        coeff_const  = B2 * (Sr_mag**2)
        discriminant = coeff_linear**2 - 4 * A2 * coeff_const
        if discriminant < 0:
            Vr_upper.append(np.nan)
            Vr_lower.append(np.nan)
        else:
            u1 = (-coeff_linear + np.sqrt(discriminant)) / (2 * A2)
            u2 = (-coeff_linear - np.sqrt(discriminant)) / (2 * A2)
            Vr_upper.append(np.sqrt(max(u1, 0)) if u1 >= 0 else np.nan)
            Vr_lower.append(np.sqrt(max(u2, 0)) if u2 >= 0 else np.nan)

    Vr_upper = np.array(Vr_upper)
    Vr_lower = np.array(Vr_lower)

    fig2, ax2 = plt.subplots(figsize=(10, 5))
    style_dark_plot(ax2, fig2, "PV Curve (Nose Curve)",
                    "Load Active Power P (MW)", "Receiving Voltage |Vr| (kV)")

    ax2.plot(P_range, Vr_upper, color='#1abc9c', linewidth=2.2,
             label='Upper (stable) branch')
    ax2.plot(P_range, Vr_lower, color='#f87171', linewidth=2.2,
             linestyle='--', label='Lower (unstable) branch')

    # Operating point
    ax2.axvline(Pr, color='#f59e0b', linestyle=':', linewidth=1.5,
                label=f'Current Pr = {Pr:.0f} MW')
    ax2.axhline(Vr_mag, color='#a855f7', linestyle=':', linewidth=1.5,
                label=f'Rated Vr = {Vr_mag:.0f} kV')
    ax2.plot(Pr, Vr_mag, 'o', color='#f87171', markersize=10, zorder=5,
             label='Operating point')

    # Nose point (Pmax at nose)
    valid_mask = ~np.isnan(Vr_upper) & ~np.isnan(Vr_lower)
    if valid_mask.any():
        nose_idx = np.nanargmin(np.abs(Vr_upper - Vr_lower) + 
                                np.where(valid_mask, 0, np.inf))
        ax2.plot(P_range[nose_idx], Vr_upper[nose_idx], '*',
                 color='#a855f7', markersize=14, zorder=6,
                 label=f'Nose point ≈ {P_range[nose_idx]:.0f} MW')

    leg2 = ax2.legend(fontsize=9)
    fix_legend(leg2)
    st.pyplot(fig2)
    plt.close(fig2)

# ============================================================
# TAB 3 — Voltage Profile (Long Line)
# ============================================================
with tab3:
    st.markdown("#### ⚡ Voltage Profile along the Transmission Line")

    # Use correct Z and Y per-unit-length
    if Z_line is not None and Y_line is not None:
        Z_total = Z_line
        Y_total = Y_line
    else:
        # Approximate from ABCD when entered directly
        # B ≈ Z*sinh(γl)/γ ≈ Z for short line; use B as Z_total approximation
        # This is an approximation — shown with a warning
        st.warning("⚠️ Voltage profile uses approximate Z/Y derived from B and C parameters. "
                   "For accurate profiles, use 'Enter Z,Y,length' or 'Enter r,L,C,length' mode.")
        Z_total = B  # Approximate: B ≈ Z for short line
        Y_total = C  # Approximate: C ≈ Y for short line

    if line_length > 0 and abs(Z_total) > 1e-10:
        Z_per_km = Z_total / line_length
        Y_per_km = Y_total / line_length
        gamma = np.sqrt(Z_per_km * Y_per_km)
        Zc = np.sqrt(Z_per_km / Y_per_km) if abs(Y_per_km) > 1e-15 else 1.0 + 0j

        distance = np.linspace(0, line_length, 300)
        V_mag = []
        for x in distance:
            # V(x) = Vr*cosh(γx) + Ir*Zc*sinh(γx)
            Vx = Vr * np.cosh(gamma * x) + Ir * Zc * np.sinh(gamma * x)
            V_mag.append(np.abs(Vx))

        fig3, ax3 = plt.subplots(figsize=(10, 4))
        style_dark_plot(ax3, fig3,
                        "Voltage Profile along the Transmission Line",
                        "Distance from Receiving End (km)",
                        "Voltage Magnitude |V| (kV)")

        ax3.plot(distance, V_mag, color='#60a5fa', linewidth=2.2)
        ax3.fill_between(distance, V_mag, alpha=0.15, color='#60a5fa')

        # Mark endpoints
        ax3.plot(0, V_mag[0], 'o', color='#1abc9c', markersize=9,
                 label=f'Receiving end: {V_mag[0]:.2f} kV')
        ax3.plot(distance[-1], V_mag[-1], 's', color='#f87171', markersize=9,
                 label=f'Sending end: {V_mag[-1]:.2f} kV')

        ax3.invert_xaxis()  # 0 = receiving, max = sending
        leg3 = ax3.legend(fontsize=9)
        fix_legend(leg3)
        st.pyplot(fig3)
        plt.close(fig3)
    else:
        st.info("Enter line length > 0 and valid parameters to display the voltage profile.")

# ============================================================
# TAB 4 — Detailed System Status
# ============================================================
with tab4:
    st.markdown("#### 📋 System Detailed Status")

    # Reactive condition
    if Qs > Qmax:
        q_status = "⚠️ Overexcited (Qs > Qmax)"
        q_color = "#f59e0b"
    elif Qs < Qmin:
        q_status = "⚠️ Underexcited (Qs < Qmin)"
        q_color = "#f59e0b"
    else:
        q_status = "✅ Within Reactive Limits"
        q_color = "#3fb950"

    c_s1, c_s2 = st.columns(2)
    with c_s1:
        st.markdown("**Phasors**")
        st.code(
            f"Vr = {polar_str(Vr, 'kV')}\n"
            f"Ir = {polar_str(Ir, 'kA')}\n"
            f"Vs = {polar_str(Vs, 'kV')}\n"
            f"Is = {polar_str(Is, 'kA')}"
        )
        st.markdown("**Power Summary**")
        st.code(
            f"Pr  = {Pr:.3f}  MW\n"
            f"Qr  = {Qr:.3f}  MVAR\n"
            f"Ps  = {Ps:.3f}  MW\n"
            f"Qs  = {Qs:.3f}  MVAR\n"
            f"Loss = {Loss_MW:.3f} MW\n"
            f"Eff  = {Efficiency:.3f} %\n"
            f"VReg = {Voltage_reg:.3f} %"
        )
        st.markdown(
            f"<span style='color:{q_color}; font-family:monospace;'>{q_status}</span>",
            unsafe_allow_html=True
        )

    with c_s2:
        st.markdown("**ABCD Parameters**")
        st.code(
            f"A = {polar_str(A)}\n"
            f"B = {polar_str(B, 'Ω')}\n"
            f"C = {polar_str(C, 'S')}\n"
            f"D = {polar_str(D)}"
        )
        st.markdown("**Power Circle Parameters**")
        st.code(
            f"β - α        = {np.degrees(beta_alpha):.4f}°\n"
            f"Circle Radius = {R_circle:.4f} MVA\n"
            f"Nrx (center) = {Nrx:.4f} MW\n"
            f"Nry (center) = {Nry:.4f} MVAR\n"
            f"Nsx (center) = {Nsx:.4f} MW\n"
            f"Nsy (center) = {Nsy:.4f} MVAR\n"
            f"Pmax         = {Pmax_calc:.4f} MW\n"
            f"Loading      = {Loading:.2f}%\n"
            f"Margin       = {Margin:.4f} MW"
        )

# ============================================================
# TAB 5 — Saved Cases
# ============================================================
with tab5:
    st.markdown("#### 💾 Saved Operating Conditions")
    if len(st.session_state.history) > 0:
        headers = st.columns([1.5, 1, 1, 1, 1, 1, 1, 1.5, 1])
        labels = ["Mode", "Pr (MW)", "PF", "Ps (MW)", "Qs (MVAR)", "Qr (MVAR)", "Load %", "Stability", "Del"]
        for i, h in enumerate(headers):
            h.markdown(f"**{labels[i]}**")
        st.markdown("<hr style='margin:6px 0;'>", unsafe_allow_html=True)

        for i, case in enumerate(st.session_state.history):
            cols = st.columns([1.5, 1, 1, 1, 1, 1, 1, 1.5, 1])
            cols[0].write(case["Mode"])
            cols[1].write(f"{case['Pr (MW)']}")
            cols[2].write(f"{case['PF']}")
            cols[3].write(f"{case['Ps (MW)']}")
            cols[4].write(f"{case['Qs (MVAR)']}")
            cols[5].write(f"{case['Qr (MVAR)']}")
            cols[6].write(f"{case['Loading %']}")
            s_color = "#3fb950" if case['Stability'] == "Stable" else "#f85149"
            cols[7].markdown(
                f"<span style='color:{s_color}; font-weight:bold; font-family:monospace;'>{case['Stability']}</span>",
                unsafe_allow_html=True
            )
            if cols[8].button("🗑️", key=f"del_{i}"):
                st.session_state.history.pop(i)
                st.rerun()
    else:
        st.info("No cases saved yet. Adjust parameters and click 'Save Condition'.")
