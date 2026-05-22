import streamlit as st
import numpy as np
import cmath
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# ======================================
# 1. Page Config & Styling
# ======================================
st.set_page_config(page_title="Power Transmission Analysis Pro", layout="wide", page_icon="⚡")

st.markdown("""
    <style>
    .stApp { background: linear-gradient(135deg, #151928 0%, #0b0f19 100%); color: #ffffff; }
    h1, h2, h3, h4 { color: #ffffff !important; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; text-align: center; }
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
    fig.patch.set_facecolor('#0e1117')
    ax.set_facecolor('#0e1117')
    ax.tick_params(colors='#ffffff')
    for spine in ax.spines.values():
        spine.set_color('#444a5e')
    ax.grid(True, linestyle=':', color='#444a5e', alpha=0.3)
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
st.title("⚡ Transmission Line Analysis & Performance Tool")
st.markdown("<hr>", unsafe_allow_html=True)

# ======================================
# 4. Inputs Panel
# ======================================
st.markdown("<span style='color:#00d2ff; font-weight:bold;'>1. Operating Conditions</span>", unsafe_allow_html=True)
c1, c2, c3, c4 = st.columns(4)
with c1:
    mode = st.selectbox("Operating Mode", ["Load", "No Load"])
with c2:
    Pr_mw = st.number_input("Receiving Power Pr (MW)", value=40.0) if mode == "Load" else 0.0
with c3:
    pfr_val = st.slider("Power Factor (PF r)", 0.0, 1.0, 0.8) if mode == "Load" else 1.0
with c4:
    pf_type = st.selectbox("PF Type", ["Lagging", "Leading"]) if mode == "Load" else "Lagging"

st.markdown("<br><span style='color:#00d2ff; font-weight:bold;'>2. Line Parameters</span>", unsafe_allow_html=True)
c5, c6, c7, c8 = st.columns(4)
with c5:
    method = st.selectbox("Parameter Input Method", ["Enter Z,Y,length (Per km)", "Enter A,B,C directly"])
    line_length = st.number_input("Line Length (km)", value=100.0, min_value=0.1)

with c6:
    if method == "Enter A,B,C directly":
        A_mag = st.number_input("|A|", value=0.986)
        A_ang = st.number_input("Angle A (°)", value=0.1)
    else:
        r_per_km = st.number_input("Resistance r (Ω/km)", value=0.1, format="%.4f")

with c7:
    if method == "Enter A,B,C directly":
        B_mag = st.number_input("|B| (Ω)", value=60.0)
        B_ang = st.number_input("Angle B (°)", value=83.0)
    else:
        xl_per_km = st.number_input("Reactance xl (Ω/km)", value=0.5, format="%.4f")

with c8:
    if method == "Enter A,B,C directly":
        C_mag = st.number_input("|C| (S)", value=0.0004, format="%.5f")
        C_ang = st.number_input("Angle C (°)", value=90.0)
    else:
        y_per_km = st.number_input("Admittance y (S/km)", value=0.000003, format="%.7f")

st.markdown("<br><span style='color:#00d2ff; font-weight:bold;'>3. Constraints & Actions</span>", unsafe_allow_html=True)
c9, c10, c11, c12 = st.columns(4)
with c9:
    vr_kv = st.number_input("Receiving Voltage Vr (kV L-L)", value=132.0)
with c10:
    Qmax = st.number_input("Qmax (MVAR)", value=300.0)
with c11:
    Qmin = st.number_input("Qmin (MVAR)", value=-300.0)
with c12:
    save_flag = st.button("💾 Save Condition")

# ======================================
# 5. Core Calculations (Rigorous Modeling)
# ======================================
if method == "Enter Z,Y,length (Per km)":
    r_total = r_per_km * line_length
    xl_total = xl_per_km * line_length
    y_total_val = y_per_km * line_length
    
    z_complex = complex(r_total, xl_total)
    y_complex = complex(0, y_total_val)

    z_mag, z_angle = cmath.polar(z_complex)
    y_mag, y_angle = cmath.polar(y_complex)
    theta_mag = np.sqrt(z_mag * y_mag)
    theta_angle_rad = (z_angle + y_angle) / 2
    theta = cmath.rect(theta_mag, theta_angle_rad)

    A = D = cmath.cosh(theta)
    B = z_complex * (cmath.sinh(theta) / theta) if theta_mag > 1e-10 else z_complex
    C = y_complex * (cmath.sinh(theta) / theta) if theta_mag > 1e-10 else y_complex
else:
    A = D = cmath.rect(A_mag, np.radians(A_ang))
    B = cmath.rect(B_mag, np.radians(B_ang))
    C = cmath.rect(C_mag, np.radians(C_ang))
    r_total, xl_total, y_total_val = B.real, B.imag, C.imag

vr_ph = (vr_kv * 1000) / np.sqrt(3) 

phi_r = np.arccos(pfr_val)
if pf_type == "Lagging":
    ir_vec = cmath.rect((Pr_mw * 1e6) / (3 * vr_ph * pfr_val), -phi_r) if pfr_val > 0 else 0j
    qr_val = -Pr_mw * np.tan(phi_r) if pfr_val > 0 else 0.0
else:
    ir_vec = cmath.rect((Pr_mw * 1e6) / (3 * vr_ph * pfr_val), phi_r) if pfr_val > 0 else 0j
    qr_val = Pr_mw * np.tan(phi_r) if pfr_val > 0 else 0.0

vs_vec = A * vr_ph + B * ir_vec
is_vec = C * vr_ph + D * ir_vec

vs_ll = (abs(vs_vec) * np.sqrt(3)) / 1000
vs_angle = np.degrees(cmath.phase(vs_vec))

s_send = 3 * vs_vec * is_vec.conjugate()
ps_mw = s_send.real / 1e6
qs_mvar_val = s_send.imag / 1e6

phi_s_rad = cmath.phase(vs_vec) - cmath.phase(is_vec)
pfs_val = np.cos(phi_s_rad)
pfs_type = "Lagging" if phi_s_rad > 0 else "Leading"

v_no_load = abs(vs_vec) / abs(A)
v_reg = (v_no_load - vr_ph) / vr_ph * 100
efficiency = (Pr_mw / ps_mw) * 100 if ps_mw > 1e-4 else 0.0
loss = ps_mw - Pr_mw

# ======================================
# Power Circle Diagram Geometry Calculations
# ======================================
Vr_complex = complex(vr_ph, 0)
const_r = 3 * ((A * (abs(Vr_complex)**2)) / B) / 1e6
nr_x, nr_y = -const_r.real, -const_r.imag
radius_r = np.sqrt((Pr_mw - nr_x)**2 + (qr_val - nr_y)**2)

Vs_mag_ph = abs(vs_vec)
const_s = 3 * ((D * (Vs_mag_ph**2)) / B) / 1e6
ns_x, ns_y = const_s.real, const_s.imag
qs_plot = -abs(qs_mvar_val) if pf_type == "Lagging" else abs(qs_mvar_val)
radius_s = np.sqrt((ps_mw - ns_x)**2 + (qs_plot - ns_y)**2)

Pmax_calc = nr_x + radius_r
stability = "Stable" if Pr_mw < Pmax_calc else "Unstable"
loading_pc = (Pr_mw / Pmax_calc * 100) if Pmax_calc > 0 else 0.0

# ======================================
# 6. Session History Save Execution
# ======================================
if save_flag:
    st.session_state.history.append({
        "Mode": mode, "Pr (MW)": round(Pr_mw, 2), "PF": round(pfr_val, 3),
        "Ps (MW)": round(ps_mw, 2), "Qs (MVAR)": round(qs_mvar_val, 2),
        "Loading %": round(loading_pc, 1), "Stability": stability
    })

# ======================================
# 7. Main KPI Dashboards
# ======================================
st.markdown("<hr>", unsafe_allow_html=True)
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Sending Power (Ps)", f"{ps_mw:.2f} MW")
m2.metric("Power Losses", f"{loss:.2f} MW")
m3.metric("Efficiency (η)", f"{efficiency:.2f} %")
m4.metric("Voltage Regulation", f"{v_reg:.2f} %")

stab_color = "#5cb85c" if stability == "Stable" else "#d9534f"
m5.markdown(f"""
    <div data-testid="metric-container" style="border-left:5px solid {stab_color};">
        <label style="color:#a0a5b5;font-size:15px;font-weight:600;">System Stability</label>
        <div style="color:{stab_color};font-weight:800;font-size:1.8rem;margin-top:-5px;">{stability}</div>
    </div>
""", unsafe_allow_html=True)

# ======================================
# 8. Analytical Workspace Tabs
# ======================================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "📊 Power Circle Diagrams", "📈 Voltage Stability (PV Curve)",
    "⚡ Distributed Voltage Profile", "📋 Complete Mathematical Status", "💾 Saved Simulation Cases"
])

# ─────────────────────────────────────────────────────────
# TAB 1 — Power Circle Diagrams Workspace
# ─────────────────────────────────────────────────────────
with tab1:
    # 1. Interactive Overlaid Plotly Control Center
    st.subheader("Interactive Power Circle Diagram (Plotly Engine)")
    fig_plotly = go.Figure()

    def add_plotly_circle(fig, cx, cy, r, name, color):
        t_arr = np.linspace(0, 2*np.pi, 200)
        fig.add_trace(go.Scatter(
            x=cx + r*np.cos(t_arr), y=cy + r*np.sin(t_arr),
            mode='lines', name=name, line=dict(color=color, dash='dash', width=1.5),
            hoverinfo='skip'
        ))

    add_plotly_circle(fig_plotly, nr_x, nr_y, radius_r, "Receiving Circle", "purple")
    add_plotly_circle(fig_plotly, ns_x, ns_y, radius_s, "Sending Circle", "orange")

    # Receiving End Geometry
    fig_plotly.add_trace(go.Scatter(x=[nr_x, 0], y=[nr_y, 0], mode='lines+markers', name="R-Base Line ((A/B)Vr²)", line=dict(color='green', width=2.5)))
    fig_plotly.add_trace(go.Scatter(x=[nr_x, Pr_mw], y=[nr_y, qr_val], mode='lines', name="R-Radius Line (VrVs/B)", line=dict(color='purple', width=2)))
    fig_plotly.add_trace(go.Scatter(x=[0, Pr_mw], y=[0, qr_val], mode='lines+markers', name="Sr Vector (Load)", line=dict(color='red', width=3.5)))

    # Sending End Geometry
    fig_plotly.add_trace(go.Scatter(x=[ns_x, 0], y=[ns_y, 0], mode='lines+markers', name="S-Base Line ((D/B)Vs²)", line=dict(color='blue', width=2.5)))
    fig_plotly.add_trace(go.Scatter(x=[ns_x, ps_mw], y=[ns_y, qs_plot], mode='lines', name="S-Radius Line (VsVr/B)", line=dict(color='orange', width=2)))
    fig_plotly.add_trace(go.Scatter(x=[0, ps_mw], y=[0, qs_plot], mode='lines+markers', name="Ss Vector (Source)", line=dict(color='darkgreen', width=3.5)))

    # Intersections & Key Reference Node Marks
    fig_plotly.add_trace(go.Scatter(x=[nr_x, ns_x], y=[nr_y, ns_y], mode='markers+text', 
                                    text=["n (R-Center)", "n' (S-Center)"], textposition="top center",
                                    marker=dict(color='yellow', size=10), name="Coordinates Centers"))

    fig_plotly.update_layout(
        template="plotly_dark", width=1000, height=750,
        xaxis=dict(title="Active Power P (MW)", gridcolor='gray', zerolinecolor='white'),
        yaxis=dict(title="Reactive Power Q (MVAR)", gridcolor='gray', zerolinecolor='white', scaleanchor="x", scaleratio=1),
        hovermode="closest"
    )
    st.plotly_chart(fig_plotly, use_container_width=True)

    st.markdown("---")
    
    # 2. Complete High-Resolution Matplotlib Visualizations
    st.subheader("Analytical Power Circle Plots (Matplotlib Engine)")
    p_choice = st.selectbox("Select Display Topology view:", ["Overlaid Combined Layout", "Isolated Receiving-End Chart", "Isolated Sending-End Chart"])
    
    circle_theta = np.linspace(0, 2 * np.pi, 500)
    
    if p_choice == "Isolated Receiving-End Chart":
        fig_m, ax_m = plt.subplots(figsize=(8, 8))
        style_dark_plot(ax_m, fig_m, "Receiving-End Power Circle Diagram", "Active Power P (MW)", "Reactive Power Q (MVAR)")
        ax_m.plot(nr_x + radius_r * np.cos(circle_theta), nr_y + radius_r * np.sin(circle_theta), color='purple', linestyle='--', linewidth=1.5, label='Power Circle')
        ax_m.plot([nr_x, 0], [nr_y, 0], color='green', linewidth=2.5, label=r'Base Line $\frac{A}{B}V_r^2$')
        ax_m.plot([nr_x, Pr_mw], [nr_y, qr_val], color='purple', linewidth=2.5, label=r'Radius Line $\frac{V_r V_s}{B}$')
        ax_m.annotate('', xy=(Pr_mw, qr_val), xytext=(0, 0), arrowprops=dict(arrowstyle='->', color='darkred', lw=3, label='$S_r$ vector'))
        ax_m.scatter([nr_x, 0, Pr_mw], [nr_y, 0, qr_val], color='yellow', zorder=5)
        ax_m.text(nr_x, nr_y, ' n (Center)', color='yellow', fontweight='bold', ha='right')
        ax_m.text(Pr_mw, qr_val, f' k ({Pr_mw}MW, {qr_val:.1f}MVAR)', color='cyan', fontweight='bold')
        ax_m.axhline(0, color='red', alpha=0.5); ax_m.axvline(0, color='red', alpha=0.5)
        ax_m.set_aspect('equal', adjustable='box')
        leg = ax_m.legend(loc='upper left'); fix_legend_color(leg)
        st.pyplot(fig_m)

    elif p_choice == "Isolated Sending-End Chart":
        fig_m, ax_m = plt.subplots(figsize=(8, 8))
        style_dark_plot(ax_m, fig_m, "Sending-End Power Circle Diagram", "Sending Power Ps (MW)", "Sending Reactive Power Qs (MVAR)")
        ax_m.plot(ns_x + radius_s * np.cos(circle_theta), ns_y + radius_s * np.sin(circle_theta), color='orange', linestyle='--', linewidth=1.5, label='Sending Circle')
        ax_m.plot([ns_x, 0], [ns_y, 0], color='blue', linewidth=2.5, label=r'Base Line $\frac{D}{B}V_s^2$')
        ax_m.plot([ns_x, ps_mw], [ns_y, qs_plot], color='orange', linewidth=2.5, label=r'Radius Line $\frac{V_s V_r}{B}$')
        ax_m.annotate('', xy=(ps_mw, qs_plot), xytext=(0, 0), arrowprops=dict(arrowstyle='->', color='darkgreen', lw=3))
        ax_m.scatter([ns_x, 0, ps_mw], [ns_y, 0, qs_plot], color='yellow', zorder=5)
        ax_m.text(ns_x, ns_y, " n' (Center)", color='yellow', fontweight='bold')
        ax_m.text(ps_mw, qs_plot, f" k' ({ps_mw:.1f}MW, {qs_plot:.1f}MVAR)", color='lime', fontweight='bold')
        ax_m.axhline(0, color='red', alpha=0.5); ax_m.axvline(0, color='red', alpha=0.5)
        ax_m.set_aspect('equal', adjustable='box')
        leg = ax_m.legend(loc='upper right'); fix_legend_color(leg)
        st.pyplot(fig_m)

    else:
        fig_m, ax_m = plt.subplots(figsize=(10, 10))
        style_dark_plot(ax_m, fig_m, "Combined Power Circle Diagram", "Active Power P (MW)", "Reactive Power Q (MVAR)")
        ax_m.axhline(0, color='red', linewidth=1.2); ax_m.axvline(0, color='red', linewidth=1.2)
        
        # Overlay Trajectories
        ax_m.plot([nr_x, 0], [nr_y, 0], color='green', linewidth=2, label=r'Receiving Base: $\frac{A}{B}V_r^2$')
        ax_m.plot([nr_x, Pr_mw], [nr_y, qr_val], color='purple', linewidth=2)
        ax_m.plot([0, Pr_mw], [0, qr_val], color='darkred', linewidth=3, label='$S_r$ (Load Vector)')
        
        ax_m.plot([ns_x, 0], [ns_y, 0], color='blue', linewidth=2, label=r'Sending Base: $\frac{D}{B}V_s^2$')
        ax_m.plot([ns_x, ps_mw], [ns_y, qs_plot], color='orange', linewidth=2)
        ax_m.plot([0, ps_mw], [0, qs_plot], color='darkgreen', linewidth=3, label='$S_s$ (Sending Vector)')
        
        ax_m.scatter([nr_x, ns_x, Pr_mw, ps_mw], [nr_y, ns_y, qr_val, qs_plot], color='yellow', zorder=5)
        ax_m.text(nr_x, nr_y, ' n (R-Center)', color='white', fontweight='bold', ha='right')
        ax_m.text(ns_x, ns_y, " n' (S-Center)", color='white', fontweight='bold', ha='left')
        
        # Highlight Boundaries Constraint Lines
        ax_m.axhline(Qmax, ls=':', color='#e14eca', alpha=0.7, label=f'Qmax Limit ({Qmax} MVAR)')
        ax_m.axhline(Qmin, ls=':', color='#e14eca', alpha=0.7, label=f'Qmin Limit ({Qmin} MVAR)')
        
        # Max Load Point Reference
        ax_m.scatter(Pmax_calc, 0, color='#e14eca', s=120, marker='^', zorder=6)
        ax_m.text(Pmax_calc, 5, f'Pmax = {Pmax_calc:.1f} MW', color='#e14eca', fontweight='bold')
        
        ax_m.set_aspect('equal', adjustable='box')
        max_dim = max(abs(nr_x), abs(ns_x), ps_mw, abs(nr_y), abs(ns_y), abs(qs_plot)) * 1.4
        ax_m.set_xlim(-max_dim, max_dim)
        ax_m.set_ylim(-max_dim, max_dim)
        leg = ax_m.legend(loc='best', prop={'size': 9}); fix_legend_color(leg)
        st.pyplot(fig_m)

# ─────────────────────────────────────────────────────────
# TAB 2 — Power Voltage (PV Curve) Stability Workspace
# ─────────────────────────────────────────────────────────
with tab2:
    st.markdown("#### 📈 Voltage Stability Analysis (PV Curve)")
    if Pmax_calc > 0:
        P_range = np.linspace(0, Pmax_calc * 1.15, 80)
        V_list = []
        for P_t in P_range:
            Q_t = P_t * np.tan(phi_r) if pf_type == "Lagging" else -P_t * np.tan(phi_r)
            Sr_t_ph = ((P_t / 3) + 1j * (Q_t / 3)) * 1e6
            Ir_t_ph = np.conj(Sr_t_ph) / np.conj(vr_ph) if vr_ph > 0 else 0j
            Vs_t_ph = A * vr_ph + B * Ir_t_ph
            V_list.append(abs(Vs_t_ph) * np.sqrt(3) / 1000)

        fig_pv, ax_pv = plt.subplots(figsize=(10, 4))
        style_dark_plot(ax_pv, fig_pv, "Voltage Profile Tracking Curve", "Load Power P (MW)", "Sending Voltage |Vs| L-L (kV)")
        ax_pv.plot(P_range, V_list, color='#e14eca', lw=2.5)
        ax_pv.axvline(Pmax_calc, ls='--', color='#ff8d72', label=f"Critical Power Point (Pmax = {Pmax_calc:.1f} MW)")
        ax_pv.scatter(Pr_mw, vs_ll, color='#00d2ff', s=100, zorder=5, label="Operating Point")
        leg_pv = ax_pv.legend(); fix_legend_color(leg_pv)
        st.pyplot(fig_pv)
    else:
        st.error("Error computing stability path trace.")

# ─────────────────────────────────────────────────────────
# TAB 3 — Hyperbolic Line Profile Analysis Workspace
# ─────────────────────────────────────────────────────────
with tab3:
    st.markdown("#### ⚡ Distributed Voltage Profile Along the Line")
    if method == "Enter Z,Y,length (Per km)" and line_length > 0:
        gamma = cmath.sqrt((r_per_km + 1j * xl_per_km) * (1j * y_per_km))
        Zc = cmath.sqrt((r_per_km + 1j * xl_per_km) / (1j * y_per_km)) if y_per_km > 0 else 1.0
        
        distance_line = np.linspace(0, line_length, 150)
        V_profile_array = []
        for space in distance_line:
            V_local = vr_ph * cmath.cosh(gamma * space) + ir_vec * Zc * cmath.sinh(gamma * space)
            V_profile_array.append(abs(V_local) * np.sqrt(3) / 1000)
            
        fig_prof, ax_prof = plt.subplots(figsize=(10, 4))
        style_dark_plot(ax_prof, fig_prof, "Voltage Magnitude Profile", "Distance from Receiving End (km)", "Voltage Magnitude L-L (kV)")
        ax_prof.plot(distance_line, V_profile_array, color='#00d2ff', lw=2.5)
        ax_prof.invert_xaxis() 
        st.pyplot(fig_prof)
    else:
        st.info("Distributed profiles tracking requires 'Per-km' parameters modeling method.")

# ─────────────────────────────────────────────────────────
# TAB 4 — Mathematical Verification Summary Reports 
# ─────────────────────────────────────────────────────────
with tab4:
    st.markdown("#### 📋 Internal Mathematical Status & Phasor Matrices")
    cs1, cs2 = st.columns(2)
    with cs1:
        st.markdown("**Electrical Phasor Quantities:**")
        st.code(
            f"Vr (Phase Base)  = {polar_str(vr_ph, 'V')}\n"
            f"Ir (Line Vector) = {polar_str(ir_vec, 'A')}\n"
            f"Vs (Phase Vector)= {polar_str(vs_vec, 'V')}\n"
            f"Is (Line Current)= {polar_str(is_vec, 'A')}\n"
            f"Vs (Line-To-Line)= {vs_ll:.4f} kV\n"
            f"Vr (Line-To-Line)= {vr_kv:.4f} kV"
        )
        st.markdown("**Transmission Performance Metrics:**")
        st.code(
            f"Transmission Angle (δ)   = {vs_angle:.4f}°\n"
            f"Voltage Regulation (VR%) = {v_reg:.4f} %\n"
            f"Efficiency (η%)          = {efficiency:.4f} %\n"
            f"Active Power Losses (ΔP) = {loss:.4f} MW"
        )
    with cs2:
        st.markdown("**ABCD Parameters Network Verification:**")
        st.code(
            f"A = D = {polar_str(A)}\n"
            f"B     = {polar_str(B, 'Ω')}\n"
            f"C     = {polar_str(C, 'S')}"
        )
        st.markdown("**Power Circle Geometry Center Matrix:**")
        st.code(
            f"Center Coordinates n_r  = ({nr_x:.3f} MW, {nr_y:.3f} MVAR)\n"
            f"Center Coordinates n_s  = ({ns_x:.3f} MW, {ns_y:.3f} MVAR)\n"
            f"Calculated Radius       = {radius_r:.3f} MVA\n"
            f"Maximum Capability (Pmax)= {Pmax_calc:.4f} MW\n"
            f"Loading Level           = {loading_pc:.1f} %"
        )

# ─────────────────────────────────────────────────────────
# TAB 5 — Saved Simulation Engineering Logs 
# ─────────────────────────────────────────────────────────
with tab5:
    st.markdown("#### 💾 Historic Parameters Log")
    if st.session_state.history:
        headers_col = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
        for slot, label in zip(headers_col, ["Mode","Pr (MW)","PF r","Ps (MW)","Qs (MVAR)","Loading %","Stability Status","Action"]):
            slot.markdown(f"**{label}**")
        st.markdown("<hr style='margin:5px 0'>", unsafe_allow_html=True)
        for idx, case in enumerate(st.session_state.history):
            slots = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
            slots[0].write(case["Mode"])
            slots[1].write(case["Pr (MW)"])
            slots[2].write(case["PF"])
            slots[3].write(case["Ps (MW)"])
            slots[4].write(case["Qs (MVAR)"])
            slots[5].write(case["Loading %"])
            status_color = "#5cb85c" if case["Stability"] == "Stable" else "#d9534f"
            slots[6].markdown(f"<span style='color:{status_color}; font-weight:bold;'>{case['Stability']}</span>", unsafe_allow_html=True)
            if slots[7].button("🗑️", key=f"delete_node_{idx}"):
                st.session_state.history.pop(idx)
                st.rerun()
    else:
        st.info("No active configurations captured in logs database records yet.")
