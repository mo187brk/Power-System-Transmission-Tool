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

if run_anim:
    st.subheader("Operating Point Movement (Stability Animation)")

    # قيم الحمل للتشغيل التدريجي
    load_values = np.linspace(0, Pmax_calc*1.2, 40)

    fig_anim, ax_anim = plt.subplots(figsize=(8,8))

    # رسم الدوائر
    ax_anim.plot(Pr_circle, Qr_circle, 'b', linewidth=2, label='Receiving')
    ax_anim.plot(Ps_circle, Qs_circle, 'r', linewidth=2, label='Sending')

    # النقاط على المسار + تلوين حسب الاستقرار
    for P_test in load_values:
        Q_test = P_test * np.tan(phi) if PF != 0 else 0
        color = 'green' if P_test <= Pmax_calc else 'red'
        ax_anim.scatter(P_test, Q_test, color=color, s=25, zorder=5)

    # Pmax كنقطة مميزة
    ax_anim.scatter(Pmax_calc, 0, color='purple', s=80, zorder=6)
    ax_anim.text(Pmax_calc, 0, f'  Pmax\n{Pmax_calc:.0f} MW', color='purple', fontsize=12)

    # Power vectors (نقطة التشغيل الحالية)
    scale = 0.04 * max(Rr, Rs)
    ax_anim.arrow(0, 0, Pr, Qr, head_width=scale, length_includes_head=True, color='blue')
    ax_anim.arrow(0, 0, Ps, Qs, head_width=scale, length_includes_head=True, color='red')

    # خطوط الإسقاط projection lines
    ax_anim.plot([Pr, Pr], [0, Qr], 'b--')
    ax_anim.plot([0, Pr], [Qr, Qr], 'b--')
    ax_anim.plot([Ps, Ps], [0, Qs], 'r--')
    ax_anim.plot([0, Ps], [Qs, Qs], 'r--')

    # خط الربط Power Transfer Line
    ax_anim.plot([Pr, Ps], [Qr, Qs], 'g', linewidth=2, label="Power Transfer Line")

    # Smart zoom
    focus_P = max(abs(Pr), abs(Ps), Pmax_calc)
    focus_Q = max(abs(Qr), abs(Qs))
    window = max(focus_P, focus_Q) * 1.5
    ax_anim.set_xlim(-window, window)
    ax_anim.set_ylim(-window, window)

    # Formatting
    ax_anim.set_aspect('equal')
    ax_anim.grid(True, linestyle='--', alpha=0.3)
    ax_anim.set_title("Stability Path Animation", fontsize=14)
    ax_anim.set_xlabel("Active Power P (MW)")
    ax_anim.set_ylabel("Reactive Power Q (MVAR)")
    ax_anim.legend()

    st.pyplot(fig_anim)
# ======================================
# Professional Plot with arrows
# ======================================
fig, ax = plt.subplots(figsize=(8,8))
ax.plot(Pr_circle, Qr_circle, 'b', label="Receiving")
ax.plot(Ps_circle, Qs_circle, 'r', label="Sending")
ax.axhline(Qmax, linestyle='--', color='orange', label="Qmax")
ax.axhline(Qmin, linestyle='--', color='orange', label="Qmin")

# Operating points
ax.scatter(Pr, Qr, color='blue', s=80, label="Pr,Qr")
ax.scatter(Ps, Qs, color='red', s=80, label="Ps,Qs")

# Vectors (arrows)
scale = 0.04 * max(Rr, Rs)
ax.arrow(0, 0, Pr, Qr, head_width=scale, color='blue', length_includes_head=True)
ax.arrow(0, 0, Ps, Qs, head_width=scale, color='red', length_includes_head=True)

# Projection lines
ax.plot([Pr, Pr], [0, Qr], 'b--')
ax.plot([0, Pr], [Qr, Qr], 'b--')
ax.plot([Ps, Ps], [0, Qs], 'r--')
ax.plot([0, Ps], [Qs, Qs], 'r--')

# Transfer line
ax.plot([Pr, Ps], [Qr, Qs], 'g', label="Power Transfer")

# Stability limit
ax.scatter(Pmax_calc, 0, color='purple', s=80)
ax.text(Pmax_calc, 0, "Pmax")

# Stability angle
delta = np.degrees(np.angle(Vs) - np.angle(Vr))
ax.text(Pr*0.6, Qr*0.6, f"δ={delta:.1f}°")

ax.set_xlim(-window, window)
ax.set_ylim(-window, window)
ax.grid()
ax.legend()
ax.set_aspect('equal')
ax.set_xlabel("P (MW)")
ax.set_ylabel("Q (MVAR)")
ax.set_title("Combined Power Circle with Vectors")
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





