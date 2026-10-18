import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math

# Intentar importar el modelo UTCI oficial
try:
    from pythermalcomfort.models import utci as calc_utci
    PYTHERMALCOMFORT_OK = True
except ImportError:
    PYTHERMALCOMFORT_OK = False

# Configuración técnica de la interfaz
st.set_page_config(
    page_title="Analizador Bioclimático: WBGT & UTCI",
    page_icon="🌡️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- FUNCIONES CIENTÍFICAS AUXILIARES ---

def calcular_tmrt_iso7726(tg: float, ta: float, va: float, d: float = 0.15, emissivity: float = 0.95) -> float:
    """
    Calcula rigurosamente la Temperatura Media Radiante (Tmrt) según la norma ISO 7726
    a partir de la temperatura de globo (Tg), temperatura de aire (Ta) y velocidad del aire (va).
    d = diámetro del globo estándar (0.15 m = 150 mm)
    emissivity = emisividad del negro mate (0.95)
    """
    tg_k = tg + 273.15
    ta_k = ta + 273.15

    # Convección forzada vs natural según ISO 7726
    h_forced = 1.1e8 * (va ** 0.6) / (emissivity * (d ** 0.4))
    h_natural = (0.25e8 / emissivity) * (abs(tg - ta) / d) ** 0.25
    h_conv = max(h_forced, h_natural)

    val = (tg_k ** 4) + h_conv * (tg - ta)
    if val > 0:
        return (val ** 0.25) - 273.15
    return tg

def clasificar_utci(valor_utci: float):
    """
    Clasificación oficial de estrés térmico de la escala UTCI (WMO / Cost Action 730).
    """
    if valor_utci > 46.0:
        return "Estrés térmico extremo por calor", "#7B1FA2", "🔴 Peligro Vital Inminente"
    elif 38.0 < valor_utci <= 46.0:
        return "Estrés térmico muy fuerte por calor", "#D32F2F", "🔴 Muy Fuerte"
    elif 32.0 < valor_utci <= 38.0:
        return "Estrés térmico fuerte por calor", "#F57C00", "🟠 Fuerte"
    elif 26.0 < valor_utci <= 32.0:
        return "Estrés térmico moderado por calor", "#FBC02D", "🟡 Moderado"
    elif 9.0 <= valor_utci <= 26.0:
        return "Sin estrés térmico (Zona de Confort)", "#388E3C", "🟢 Confort Térmico"
    elif 0.0 <= valor_utci < 9.0:
        return "Estrés térmico ligero por frío", "#0288D1", "🔵 Ligero Frío"
    elif -13.0 <= valor_utci < 0.0:
        return "Estrés térmico moderado por frío", "#1976D2", "🔵 Moderado Frío"
    elif -27.0 <= valor_utci < -13.0:
        return "Estrés térmico fuerte por frío", "#303F9F", "🔷 Fuerte Frío"
    elif -40.0 <= valor_utci < -27.0:
        return "Estrés térmico muy fuerte por frío", "#512DA8", "🟣 Muy Fuerte Frío"
    else:
        return "Estrés térmico extremo por frío", "#311B92", "⬛ Peligro Extremo Frío"

# Base normativa ISO 7243 para WBGT
ISO_LIMITS_WBGT = {
    "Clase 0: Reposo (< 115 W)": {"aclimatado": 33.0, "no_aclimatado": 32.0},
    "Clase 1: Ligero (115 - 180 W)": {"aclimatado": 30.0, "no_aclimatado": 29.0},
    "Clase 2: Moderado (180 - 300 W)": {"aclimatado": 28.0, "no_aclimatado": 26.0},
    "Clase 3: Pesado (300 - 390 W)": {"aclimatado": 26.0, "no_aclimatado": 23.0},
    "Clase 4: Muy pesado (> 390 W)": {"aclimatado": 25.0, "no_aclimatado": 20.0}
}

# --- BARRA LATERAL: ENTRADA DE DATOS EXPERIMENTALES ---
with st.sidebar:
    st.header("🎛️ Datos Experimentales")

    st.markdown("#### 1. Temperaturas Básicas (°C)")
    t_a = st.number_input("Temp. Aire / Bulbo Seco ($T_a$):", value=33.5, step=0.1, format="%.2f")
    t_g = st.number_input("Temp. de Globo ($T_g$):", value=42.0, step=0.1, format="%.2f")
    t_nw = st.number_input("Temp. Bulbo Húmedo Natural ($T_{nw}$):", value=26.0, step=0.1, format="%.2f")

    st.markdown("#### 2. Variables Adicionales (UTCI / ISO 7726)")
    rh = st.slider("Humedad Relativa ($RH$, %):", min_value=5.0, max_value=100.0, value=50.0, step=1.0)
    v_a = st.number_input("Velocidad del aire a 1.2 m ($v_a$, m/s):", min_value=0.1, max_value=25.0, value=1.2, step=0.1)

    st.markdown("#### 3. Parámetros Operacionales WBGT")
    con_sol = st.checkbox("Exposición solar directa (Carga solar)", value=True)
    clase_metabolica = st.selectbox("Actividad laboral (ISO 7243):", list(ISO_LIMITS_WBGT.keys()), index=2)
    aclimatado = st.toggle("Personal Aclimatado", value=True)

# --- CÁLCULOS RIGUROSOS ---

# 1. WBGT (ISO 7243)
if con_sol:
    wbgt = 0.7 * t_nw + 0.2 * t_g + 0.1 * t_a
else:
    wbgt = 0.7 * t_nw + 0.3 * t_g
wbgt = round(wbgt, 2)

estado_str = "aclimatado" if aclimatado else "no_aclimatado"
wbgt_limite = ISO_LIMITS_WBGT[clase_metabolica][estado_str]
delta_wbgt = round(wbgt - wbgt_limite, 2)

# 2. Temperatura Media Radiante (Tmrt) según ISO 7726
t_mrt = calcular_tmrt_iso7726(tg=t_g, ta=t_a, va=v_a)

# 3. Viento a 10 m (v_10) según ley logarítmica normalizada UTCI
v_10 = max(0.5, v_a * (math.log(10.0 / 0.01) / math.log(1.2 / 0.01)))

# 4. UTCI
utci_val = None
utci_cat = ""
utci_color = "#333333"
utci_badge = ""

if PYTHERMALCOMFORT_OK:
    try:
        res = calc_utci(tdb=t_a, tr=t_mrt, v=v_10, rh=rh, limit_inputs=False)
        utci_val = round(float(res.utci if hasattr(res, 'utci') else res), 1)
        utci_cat, utci_color, utci_badge = clasificar_utci(utci_val)
    except Exception as e:
        st.error(f"Error al evaluar UTCI: {e}")
else:
    st.sidebar.warning("⚠️ Falta `pythermalcomfort`. Ejecuta en consola: `python -m pip install pythermalcomfort`")

# --- CABECERA PRINCIPAL ---
st.title("🔬 Monitor Dual de Estrés Térmico: WBGT & UTCI")
st.caption("Comparación de índices normativos de estrés térmico laboral (ISO 7243) y biometeorológico universal (COST 730 / WMO)")

# Panel de variables físicas calculadas
col_m1, col_m2, col_m3, col_m4 = st.columns(4)
col_m1.metric("Temp. Aire ($T_a$)", f"{t_a:.1f} °C")
col_m2.metric("Temp. Globo ($T_g$)", f"{t_g:.1f} °C")
col_m3.metric("Temp. Media Radiante ($T_{mrt}$)", f"{t_mrt:.1f} °C", help="Deducida por balance ISO 7726 a partir de Tg, Ta y va")
col_m4.metric("Viento a 10m ($v_{10}$)", f"{v_10:.1f} m/s", help="Proyectado desde va a 1.2m según perfil de rugosidad UTCI")

st.markdown("---")

# --- COMPARACIÓN GRÁFICA PARALELA ---
col_gauge_wbgt, col_gauge_utci = st.columns(2)

# --- COLUMNA 1: TACÓMETRO WBGT ---
with col_gauge_wbgt:
    st.subheader("1. Índice WBGT (ISO 7243)")
    fig_wbgt = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=wbgt,
        delta={'reference': wbgt_limite, 'suffix': " °C (Límite)", 'increasing': {'color': "#FF4B4B"}, 'decreasing': {'color': "#00CC96"}},
        gauge={
            'axis': {'range': [10, 42], 'tickwidth': 1.5},
            'bar': {'color': "#1565C0", 'thickness': 0.3},
            'steps': [
                {'range': [10, wbgt_limite - 2], 'color': "rgba(46, 204, 113, 0.35)"},
                {'range': [wbgt_limite - 2, wbgt_limite], 'color': "rgba(241, 196, 15, 0.45)"},
                {'range': [wbgt_limite, 42], 'color': "rgba(231, 76, 60, 0.45)"}
            ],
            'threshold': {'line': {'color': "#e74c3c", 'width': 4}, 'thickness': 0.8, 'value': wbgt_limite}
        },
        number={'suffix': " °C", 'font': {'size': 32}}
    ))
    fig_wbgt.update_layout(height=280, margin=dict(l=20, r=20, t=30, b=10))
    st.plotly_chart(fig_wbgt, use_container_width=True)

    if wbgt <= (wbgt_limite - 2.0):
        st.success(f"**CONFORME (Bajo Riesgo)**: Por debajo del límite normativo ({wbgt_limite} °C).")
    elif (wbgt_limite - 2.0) < wbgt <= wbgt_limite:
        st.warning(f"**ALERTA PREVENTIVA**: Cerca del límite ({wbgt_limite} °C). Hidratación activa.")
    else:
        st.error(f"**SOBRECARGA TÉRMICA**: Supera el umbral laboral por {abs(delta_wbgt)} °C.")

# --- COLUMNA 2: TACÓMETRO ESCALA UTCI ---
with col_gauge_utci:
    st.subheader("2. Índice UTCI (Biometeorología)")
    if utci_val is not None:
        fig_utci = go.Figure(go.Indicator(
            mode="gauge+number",
            value=utci_val,
            gauge={
                'axis': {'range': [-10, 52], 'tickwidth': 1.5, 'tickmode': 'array', 'tickvals': [-10, 0, 9, 26, 32, 38, 46, 52]},
                'bar': {'color': "#212121", 'thickness': 0.28},
                'steps': [
                    {'range': [-10, 0], 'color': "#1976D2"},     # Moderado frío
                    {'range': [0, 9], 'color': "#4FC3F7"},       # Ligero frío
                    {'range': [9, 26], 'color': "#4CAF50"},      # Confort
                    {'range': [26, 32], 'color': "#FDD835"},     # Moderado calor
                    {'range': [32, 38], 'color': "#FB8C00"},     # Fuerte calor
                    {'range': [38, 46], 'color': "#E53935"},     # Muy fuerte calor
                    {'range': [46, 52], 'color': "#8E24AA"}      # Extremo calor
                ],
                'threshold': {'line': {'color': "#212121", 'width': 4}, 'thickness': 0.8, 'value': utci_val}
            },
            number={'suffix': " °C", 'font': {'size': 32}}
        ))
        fig_utci.update_layout(height=280, margin=dict(l=20, r=20, t=30, b=10))
        st.plotly_chart(fig_utci, use_container_width=True)

        st.info(f"**Categoría:** {utci_cat} ({utci_badge})")
    else:
        st.error("Instala `pythermalcomfort` para ver el tacómetro UTCI.")

st.markdown("---")

# --- GUÍA VISUAL OFICIAL DE LA ESCALA UTCI ---
st.subheader("📊 Escala Normativa Oficial de Estrés Térmico UTCI (WMO / COST 730)")

escala_data = [
    {"Rango UTCI (°C)": "> +46 °C", "Categoría de Estrés": "Estrés térmico extremo por calor", "Efecto Fisiológico": "Fallo inminente de termorregulación, golpe de calor letal", "Color": "🟣 Púrpura"},
    {"Rango UTCI (°C)": "+38 a +46 °C", "Categoría de Estrés": "Estrés térmico muy fuerte por calor", "Efecto Fisiológico": "Sudoración masiva, aumento pronunciado de la temperatura interna", "Color": "🔴 Rojo"},
    {"Rango UTCI (°C)": "+32 a +38 °C", "Categoría de Estrés": "Estrés térmico fuerte por calor", "Efecto Fisiológico": "Sobrecarga cardiovascular moderada/fuerte", "Color": "🟠 Naranja"},
    {"Rango UTCI (°C)": "+26 a +32 °C", "Categoría de Estrés": "Estrés térmico moderado por calor", "Efecto Fisiológico": "Sensación de bochorno, incremento de sudoración", "Color": "🟡 Amarillo"},
    {"Rango UTCI (°C)": "+9 a +26 °C", "Categoría de Estrés": "Sin estrés térmico (Confort)", "Efecto Fisiológico": "Equilibrio térmico sin esfuerzo fisiológico", "Color": "🟢 Verde"},
    {"Rango UTCI (°C)": "0 a +9 °C", "Categoría de Estrés": "Estrés térmico ligero por frío", "Efecto Fisiológico": "Vasoconstricción periférica leve", "Color": "🔵 Celeste"},
    {"Rango UTCI (°C)": "-13 a 0 °C", "Categoría de Estrés": "Estrés térmico moderado por frío", "Efecto Fisiológico": "Tiritona involuntaria, extremidades frías", "Color": "🔵 Azul"},
    {"Rango UTCI (°C)": "< -13 °C", "Categoría de Estrés": "Estrés fuerte / extremo por frío", "Efecto Fisiológico": "Hipotermia y riesgo severo de congelación", "Color": "🔷 Azul Oscuro"}
]

df_escala = pd.DataFrame(escala_data)
st.dataframe(df_escala, use_container_width=True, hide_index=True)
