# ─────────────────────────────────────────────────────────
# TAB 1 — Power Circle Diagram  (ACADEMIC STYLE FIXED)
# ─────────────────────────────────────────────────────────
with tab1:

    col_plot, col_info = st.columns([3, 1])

    with col_plot:

        # ==========================================
        # Circle Coordinates
        # ==========================================
        theta_arr = np.linspace(0, 2*np.pi, 800)

        Pr_circle = n_rx + R_circle * np.cos(theta_arr)
        Qr_circle = n_ry + R_circle * np.sin(theta_arr)

        Ps_circle = n_sx + R_circle * np.cos(theta_arr)
        Qs_circle = n_sy + R_circle * np.sin(theta_arr)

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
        # Operating Labels
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
        # Radius Vectors
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
        # Origin → Operating Point Vectors
        # ==========================================
        ax.annotate(
            "",
            xy=(Pr, Qr),
            xytext=(0,0),
            arrowprops=dict(
                arrowstyle='->',
                linewidth=2.8,
                color='#00d2ff'
            )
        )

        ax.annotate(
            "",
            xy=(Ps, Qs),
            xytext=(0,0),
            arrowprops=dict(
                arrowstyle='->',
                linewidth=2.8,
                color='#ff8d72'
            )
        )

        # ==========================================
        # Vertical Projection
        # ==========================================
        ax.plot(
            [Pr, Pr],
            [0, Qr],
            linestyle='--',
            linewidth=1,
            color='white',
            alpha=0.25
        )

        # ==========================================
        # Pmax Construction
        # ==========================================
        ax.plot(
            [n_rx, Pmax_calc],
            [n_ry, 0],
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
        # PF Angle Arc
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
        # Delta Label
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
        # β−α Label
        # ==========================================
        ax.text(
            n_rx + R_circle*0.04,
            n_ry - R_circle*0.10,
            f'β−α = {beta-alpha:.1f}°',
            color='#bbbbbb',
            fontsize=11
        )

        # ==========================================
        # Window Limits
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

        # ==========================================
        # IMPORTANT FIX
        # ==========================================
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
