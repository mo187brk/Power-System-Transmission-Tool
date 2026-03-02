import streamlit as st
import numpy as np
import matplotlib.pyplot as plt


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
        color: #ffffff !important;
        border-radius: 8px !important;
        border: none !important;
        font-weight: 700 !important;
        font-size: 15px !important;
        height: 44px !important;
        width: 100% !important;
        margin-top: 27px !important;
        box-shadow: 0 4px 15px rgba(216, 60, 184, 0.4) !important;
        transition: all 0.3s ease !important;
    }
    div[data-testid="stButton"] button:hover {
        box-shadow: 0 6px 20px rgba(216, 60, 184, 0.6) !important;
        transform: translateY(-2px) !important;
        background: linear-gradient(90deg, #f87462 0%, #d83cb8 100%) !important;
    }
    /* ========================================= */
    
    hr { border-color: #2d3246; margin-top: 15px; margin-bottom: 15px; }
    </style>
""", unsafe_allow_html=True)

# Helper Function
def polar_str(value, unit=""):
    mag = np.abs(value)
    ang = np.degrees(np.angle(value))
    return f"{mag:.4f} ∠ {ang:.2f}° {unit}"

# ======================================
# 2. Session State Management
# ======================================
if 'history' not in st.session_state:
    st.session_state.history = []

# ======================================
# 3. Header & Control Panel
# ======================================
col_logo, col_title = st.columns([1, 11])

with col_title:
    st.title("⚡ Power Transmission Analysis")


st.markdown("<hr>", unsafe_allow_html=True)

st.markdown("<span style='color:#00d2ff; font-weight:bold;'>1. Operating Conditions</span>", unsafe_allow_html=True)
c1, c2, c3, c4 = st.columns(4)
with c1:
    mode = st.selectbox("Operating Mode", ["Load", "No Load"])
with c2:
    if mode == "Load":
        Pr = st.number_input("Receiving Power Pr (MW)", value=100.0)
    else:
        Pr = 0.0
with c3:
    if mode == "Load":
        PF = st.slider("Power Factor", 0.0, 1.0, 0.9)
    else:
        PF = 0.0
with c4:
    if mode == "Load":
        pf_type = st.selectbox("PF Type", ["Lagging", "Leading"])
    else:
        pf_type = "Lagging"

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
    Vr_mag = st.number_input("Receiving Voltage Vr (kV)", value=220.0)
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
# 4. Core Mathematical Calculations
# ======================================
phi = np.arccos(PF) if PF != 0 else 0
if pf_type == "Leading": phi = -phi

Pr_phase = Pr / 3
Qr = Pr * np.tan(phi) if PF != 0 else 0
Qr_phase = Qr / 3
Sr = Pr_phase + 1j * Qr_phase

Z = None; Y = None
if method == "Enter A,B,C directly":
    A = A_mag*np.exp(1j*np.radians(A_ang)); B = B_mag*np.exp(1j*np.radians(B_ang))
    C = C_mag*np.exp(1j*np.radians(C_ang)); D = A
elif method == "Enter Z,Y,length":
    Z = (R + 1j*X)*line_length; Y = (G + 1j*B_sh)*line_length
    theta = np.sqrt(Z*Y)
    A = np.cosh(theta); B = Z*np.sinh(theta)/theta; C = Y*np.sinh(theta)/theta; D = A
else:
    w = 2*np.pi*50
    Z = (r + 1j*w*L)*line_length; Y = (1j*w*C_val)*line_length
    theta = np.sqrt(Z*Y)
    A = np.cosh(theta); B = Z*np.sinh(theta)/theta; C = Y*np.sinh(theta)/theta; D = A

Vr = Vr_mag + 0j
Ir = np.conj(Sr) / np.conj(Vr) if Vr != 0 else 0

if fault == "3-Phase Fault at Receiving End":
    Vr = fault_impedance * Ir

Vs = A*Vr + B*Ir
Is = C*Vr + D*Ir

Ss = 3*Vs*np.conj(Is)
Ps = np.real(Ss); Qs = np.imag(Ss)
Loss = Ps - Pr
Efficiency = (Pr/Ps)*100 if Ps != 0 else 0
Voltage_reg = ((np.abs(Vs)-Vr_mag)/Vr_mag)*100

Pmax_calc = (np.abs(Vs)*Vr_mag)/np.abs(B) if B != 0 else 0
Loading = (Pr/Pmax_calc)*100 if Pmax_calc !=0 else 0
Margin = Pmax_calc - Pr
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
# 7. Visualizations & Tabs
# ======================================
st.markdown("<br>", unsafe_allow_html=True)
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📊 Power Circle Diagram", "📈 PV Curve", "⚡ Voltage Profile", "📋 Detailed Status", "💾 Saved Cases"])

def style_dark_plot(ax, fig, title, xlabel, ylabel):
    fig.patch.set_alpha(0.0)
    ax.patch.set_alpha(0.0)
    ax.tick_params(colors='#ffffff') 
    for spine in ax.spines.values(): spine.set_color('#444a5e')
    ax.grid(True, linestyle=':', color='#444a5e')
    ax.set_title(title, fontsize=14, fontweight='bold', color='#ffffff', pad=15)
    ax.set_xlabel(xlabel, color='#ffffff', fontweight='bold')
    ax.set_ylabel(ylabel, color='#ffffff', fontweight='bold')
    
def fix_legend_color(legend):
    if legend:
        for text in legend.get_texts():
            text.set_color('white')

with tab1:
    col_plot, col_info = st.columns([3, 1])
    with col_plot:
        alpha_ang = np.angle(A); beta_ang = np.angle(B); ratio = np.abs(A/B) if B != 0 else 0
        Rr = ratio*(Vr_mag**2); theta_r = beta_ang - alpha_ang
        Crx = Rr*np.cos(theta_r); Cry = Rr*np.sin(theta_r)
        Vs_mag = np.abs(Vs); Rs = ratio*(Vs_mag**2)
        Csx = Rs*np.cos(beta_ang); Csy = Rs*np.sin(beta_ang)
        
        theta_plot = np.linspace(0, 2*np.pi, 400)
        Pr_circle = Crx + Rr*np.cos(theta_plot); Qr_circle = Cry + Rr*np.sin(theta_plot)
        Ps_circle = Csx + Rs*np.cos(theta_plot); Qs_circle = Csy + Rs*np.sin(theta_plot)

        if run_anim:
            st.markdown("<h5 style='color:white; text-align:center;'>Operating Point Movement (Stability Animation)</h5>", unsafe_allow_html=True)
            load_values = np.linspace(0, Pmax_calc*1.2, 40)
            fig_anim, ax_anim = plt.subplots(figsize=(8,8))
            style_dark_plot(ax_anim, fig_anim, "", "Active Power P (MW)", "Reactive Power Q (MVAR)")
            
            ax_anim.plot(Pr_circle, Qr_circle, '#00d2ff', linewidth=2, label='Receiving', alpha=0.6)
            ax_anim.plot(Ps_circle, Qs_circle, '#ff8d72', linewidth=2, label='Sending', alpha=0.6)
            
            for P_test in load_values:
                Q_test = P_test * np.tan(phi) if PF != 0 else 0
                color = '#5cb85c' if P_test <= Pmax_calc else '#d9534f'
                ax_anim.scatter(P_test, Q_test, color=color, s=25, zorder=5)
                
            ax_anim.scatter(Pmax_calc, 0, color='#e14eca', s=100, zorder=6)
            ax_anim.text(Pmax_calc, 0, f' Pmax\n{Pmax_calc:.0f}', color='#e14eca', fontsize=12)
            
            focus_P = max(abs(Pr), abs(Ps), Pmax_calc)
            focus_Q = max(abs(Qr), abs(Qs))
            window_anim = max(focus_P, focus_Q) * 1.5
            ax_anim.set_xlim(-window_anim, window_anim)
            ax_anim.set_ylim(-window_anim, window_anim)
            ax_anim.set_aspect('equal')
            ax_anim.axhline(0, color='#626b82', linewidth=1.2)
            ax_anim.axvline(0, color='#626b82', linewidth=1.2)
            
            leg = ax_anim.legend(facecolor='#1e2233', edgecolor='#444a5e')
            fix_legend_color(leg)
            st.pyplot(fig_anim)
            
        else:
            fig, ax = plt.subplots(figsize=(8,8))
            style_dark_plot(ax, fig, "Combined Power Circle with Vectors", "Active Power P (MW)", "Reactive Power Q (MVAR)")
            
            ax.plot(Pr_circle, Qr_circle, '#00d2ff', linewidth=2, label="Receiving")
            ax.plot(Ps_circle, Qs_circle, '#ff8d72', linewidth=2, label="Sending")
            ax.axhline(Qmax, linestyle='--', color='#e14eca', label="Qmax", alpha=0.7)
            ax.axhline(Qmin, linestyle='--', color='#e14eca', label="Qmin", alpha=0.7)

            ax.scatter(Pr, Qr, color='#00d2ff', s=120, edgecolors='white', zorder=5)
            ax.scatter(Ps, Qs, color='#ff8d72', s=120, edgecolors='white', zorder=5)

            scale = 0.06 * max(Rr, Rs) if max(Rr, Rs) > 0 else 1
            ax.arrow(0, 0, Pr, Qr, head_width=scale, color='#00d2ff', length_includes_head=True)
            ax.arrow(0, 0, Ps, Qs, head_width=scale, color='#ff8d72', length_includes_head=True)

            ax.text(Pr, Qr, "  Sr", fontsize=12, color='#00d2ff')
            ax.text(Ps, Qs, "  Ss", fontsize=12, color='#ff8d72')
            
            ax.plot([Pr, Ps], [Qr, Qs], '#5cb85c', linewidth=2, label="Power Transfer")
            ax.scatter(Pmax_calc, 0, color='#e14eca', s=120)
            ax.text(Pmax_calc, 0, f" Pmax\n{Pmax_calc:.0f}", fontsize=11, color='#e14eca')

            delta_ang = np.degrees(np.angle(Vs) - np.angle(Vr))
            ax.text(Pr*0.5, Qr*0.5, f"δ = {delta_ang:.1f}°", color='white', fontsize=11, bbox=dict(facecolor='#1e2233', alpha=0.8, edgecolor='#444a5e'))

            focus = max(abs(Pr), abs(Ps), abs(Qr), abs(Qs))
            window = max(focus*2.2, 0.25*max(Rr, Rs)) if focus > 0 else 100
            ax.set_xlim(-window, window)
            ax.set_ylim(-window, window)
            ax.set_aspect('equal')
            ax.axhline(0, color='#626b82', linewidth=1.2)
            ax.axvline(0, color='#626b82', linewidth=1.2)
            
            leg = ax.legend(facecolor='#1e2233', edgecolor='#444a5e')
            fix_legend_color(leg)
            st.pyplot(fig)

    with col_info:
        st.markdown("#### Point Details")
        st.info(f"**Pr:** {Pr:.2f} MW")
        st.info(f"**Qr:** {Qr:.2f} MVAR")
        st.info(f"**Ps:** {Ps:.2f} MW")
        st.info(f"**Qs:** {Qs:.2f} MVAR")
        st.info(f"**Pmax:** {Pmax_calc:.2f} MW")
        st.info(f"**Margin:** {Margin:.2f} MW")

with tab2:
    st.markdown("#### PV Curve (Voltage Stability)")
    P_range = np.linspace(0, Pmax_calc*1.2, 50)
    V_list = []
    for P_test in P_range:
        Q_test = P_test * np.tan(phi) if PF != 0 else 0
        Sr_test = (P_test + 1j*Q_test)/3
        Ir_test = np.conj(Sr_test)/np.conj(Vr) if Vr !=0 else 0
        Vs_test = A*Vr + B*Ir_test
        V_list.append(np.abs(Vs_test))

    fig_pv, ax_pv = plt.subplots(figsize=(10,4))
    style_dark_plot(ax_pv, fig_pv, "PV Curve", "Load P (MW)", "Sending Voltage |Vs| (kV)")
    ax_pv.plot(P_range, V_list, color='#e14eca', linewidth=2)
    ax_pv.axvline(Pmax_calc, linestyle='--', color='#ff8d72', label=f"Pmax ({Pmax_calc:.1f} MW)")
    leg_pv = ax_pv.legend(facecolor='#1e2233', edgecolor='#444a5e')
    fix_legend_color(leg_pv)
    st.pyplot(fig_pv)

with tab3:
    st.markdown("#### Voltage Profile (Long Line Model)")
    Z_total = Z if Z is not None else B
    Y_total = Y if Y is not None else C

    Z_per_km = Z_total / line_length if line_length != 0 else 0
    Y_per_km = Y_total / line_length if line_length != 0 else 0

    gamma = np.sqrt(Z_per_km * Y_per_km)
    Zc = np.sqrt(Z_per_km / Y_per_km) if Y_per_km != 0 else 1

    distance = np.linspace(0, line_length, 200)
    V_profile = []
    for x in distance:
        Vx = Vr*np.cosh(gamma*x) + Ir*Zc*np.sinh(gamma*x)
        V_profile.append(np.abs(Vx))

    fig2, ax2 = plt.subplots(figsize=(10,4))
    style_dark_plot(ax2, fig2, "Voltage Profile along the Transmission Line", "Distance from Receiving End (km)", "Voltage (kV)")
    ax2.plot(distance, V_profile, color='#00d2ff', linewidth=2)
    ax2.invert_xaxis()
    st.pyplot(fig2)

with tab4:
    st.markdown("#### 📋 System Detailed Status")
    c_status1, c_status2 = st.columns(2)
    with c_status1:
        st.markdown("**Voltages & Currents**")
        st.code(f"Vr = {polar_str(Vr, 'kV')}\nIr = {polar_str(Ir, 'kA')}\nVs = {polar_str(Vs, 'kV')}\nIs = {polar_str(Is, 'kA')}")
        st.markdown("**Performance Metrics**")
        st.code(f"Voltage Reg  = {Voltage_reg:.2f} %\nEfficiency   = {Efficiency:.2f} %\nReactive Stat= {'Overexcited' if Qs > Qmax else 'Underexcited' if Qs < Qmin else 'Within Limits'}")
    with c_status2:
        st.markdown("**ABCD Parameters**")
        st.code(f"A = {polar_str(A)}\nB = {polar_str(B, 'Ω')}\nC = {polar_str(C, 'S')}\nD = {polar_str(D)}")

with tab5:
    st.markdown("#### 💾 Saved Cases Management")
    if len(st.session_state.history) > 0:
        headers = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
        headers[0].markdown("**Mode**"); headers[1].markdown("**Pr (MW)**")
        headers[2].markdown("**PF**"); headers[3].markdown("**Ps (MW)**")
        headers[4].markdown("**Qs (MVAR)**"); headers[5].markdown("**Load %**")
        headers[6].markdown("**Stability**"); headers[7].markdown("**Action**")
        st.markdown("<hr style='margin: 5px 0;'>", unsafe_allow_html=True)
        
        for i, case in enumerate(st.session_state.history):
            cols = st.columns([1.5, 1, 1, 1, 1, 1, 1.5, 1])
            cols[0].write(case["Mode"]); cols[1].write(f"{case['Pr (MW)']}")
            cols[2].write(f"{case['PF']}"); cols[3].write(f"{case['Ps (MW)']}")
            cols[4].write(f"{case['Qs (MVAR)']}"); cols[5].write(f"{case['Loading %']}")
            
            s_color = "#5cb85c" if case['Stability'] == "Stable" else "#d9534f"
            cols[6].markdown(f"<span style='color:{s_color}; font-weight:bold;'>{case['Stability']}</span>", unsafe_allow_html=True)
            
            if cols[7].button("🗑️", key=f"del_{i}"):
                st.session_state.history.pop(i)
                st.rerun()
    else:
        st.info("No cases saved yet. Adjust parameters and click 'Save Condition'.")
