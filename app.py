import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# إعدادات الصفحة
st.set_page_config(page_title="Power Circle Diagram - Pro Version", layout="wide")

st.markdown("""
    <style>
    .stApp { background-color: #0e1117; color: white; }
    h1 { color: #00d2ff; text-align: center; }
    </style>
    """, unsafe_allow_html=True)

st.title("⚡ Advanced Power Circle Analysis")

# --- المدخلات (Sidebar) ---
st.sidebar.header("Line Parameters (ABCD)")
Vr_kv = st.sidebar.number_input("Receiving Voltage Vr (Line-Line kV)", value=220.0)
Vs_kv = st.sidebar.number_input("Sending Voltage Vs (Line-Line kV)", value=240.0)

# ABCD Constants
A_mag = st.sidebar.number_input("|A|", value=0.98, format="%.4f")
A_ang = st.sidebar.number_input("Angle alpha (deg)", value=0.2)
B_mag = st.sidebar.number_input("|B| (Ohm)", value=110.0)
B_ang = st.sidebar.number_input("Angle beta (deg)", value=75.0)

# Operating Point
Pr_mw = st.sidebar.number_input("Load Pr (MW)", value=150.0)
pf = st.sidebar.slider("Power Factor", 0.7, 1.0, 0.9)
pf_type = st.sidebar.selectbox("Type", ["Lagging", "Leading"])

# --- الحسابات الرياضية (مطابقة للمحاضرة) ---
alpha = np.radians(A_ang)
beta = np.radians(B_ang)
phi = np.arccos(pf)
if pf_type == "Lagging": phi = -phi # الزاوية سالبة للـ Lagging في نظام الإحداثيات

# حساب نصف القطر والمراكز (Receiving End)
# Radius = |Vs|*|Vr| / |B|
radius_r = (Vs_kv * Vr_kv) / B_mag
# Center = (|A|*|Vr|^2 / |B|) at angle (beta - alpha)
center_mag_r = (A_mag * (Vr_kv**2)) / B_mag
cx_r = -center_mag_r * np.cos(beta - alpha)
cy_r = -center_mag_r * np.sin(beta - alpha)

# نقطة التشغيل الحالية
Qr_mvar = Pr_mw * np.tan(np.arccos(pf))
if pf_type == "Lagging": Qr_mvar = -Qr_mvar # تمثيل الأحمال الحثية بالأسفل

# --- الرسم البياني المحسن ---
fig, ax = plt.subplots(figsize=(10, 10))
fig.patch.set_facecolor('#0e1117')
ax.set_facecolor('#0e1117')

# رسم الدائرة
theta = np.linspace(0, 2*np.pi, 500)
x_circle = cx_r + radius_r * np.cos(theta)
y_circle = cy_r + radius_r * np.sin(theta)
ax.plot(x_circle, y_circle, color='#00d2ff', linewidth=2, label="Receiving-end Circle")

# رسم المحاور
ax.axhline(0, color='white', linewidth=1)
ax.axvline(0, color='white', linewidth=1)

# رسم متجه مركز الدائرة (من الأصل للمركز)
ax.annotate('', xy=(cx_r, cy_r), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color='yellow', linestyle='--'))
ax.text(cx_r/2, cy_r/2, '  Center Vector', color='yellow')

# رسم نقطة التشغيل (Load Point)
ax.scatter(Pr_mw, Qr_mvar, color='red', s=100, zorder=5, label="Operating Point")
ax.plot([0, Pr_mw], [0, Qr_mvar], color='red', linestyle=':')

# رسم حدود الاستقرار (الخط الأفقي المار بالمركز)
ax.axhline(cy_r, color='gray', linestyle='--', alpha=0.5, label="Theoretical Limit Line")

# تنسيق المحاور كما طلبت (P على X و Q على Y)
ax.set_xlabel("Active Power P (MW)", color='white', fontsize=12)
ax.set_ylabel("Reactive Power Q (MVAR)", color='white', fontsize=12)
ax.tick_params(colors='white')
ax.grid(True, linestyle=':', alpha=0.3)
ax.set_aspect('equal')
ax.legend()

# إضافة بيانات توضيحية داخل الرسم
ax.text(cx_r, cy_r, ' Center O\'', color='white', verticalalignment='bottom')

st.pyplot(fig)

# --- عرض النتائج الرقمية ---
st.markdown("---")
c1, c2, c3 = st.columns(3)
with c1:
    st.metric("Radius (MW)", f"{radius_r:.2f}")
with c2:
    st.metric("Center (X, Y)", f"({cx_r:.1f}, {cy_r:.1f})")
with c3:
    P_max = radius_r + cx_r
    st.metric("Static Stability Limit", f"{P_max:.2f} MW")

st.info("💡 ملاحظة للدكتور: الرسم يوضح دائرة الاستقبال حيث تقع الإحداثيات في الربع الثالث للمركز، والقدرة الفعالة موجبة (استقبال).")
