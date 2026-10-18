import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math

# --- CONFIGURACIÓN GENERAL ---
st.set_page_config(
    page_title="Monitor Bioclimático Integral | PMV · WBGT · UTCI",
    page_icon="🌡️",
    layout="wide"
)

# Inyección de tipografía técnica Arial Narrow y tarjetas de instrumentación
st.markdown("""
<style>
    html, body, [class*="css"], .stMarkdown, .stText, div, p, span {
        font-family: 'Arial Narrow', 'Segoe UI', -apple-system, sans-serif !important;
    }
    .block-header {
        background-color: #1a222d;
        border-left: 4px solid #00a8ff;
        padding: 6px 12px;
        margin-bottom: 8px;
        color: #f1f2f6;
        font-weight: 700;
        letter-spacing: 0.5px;
        text-transform: uppercase;
        font-size: 0.95rem;
    }
    .metric-box {
        background: rgba(255, 255, 255, 0.04);
        border: 1px solid rgba(255, 255, 255, 0.12);
        border-radius: 4px;
        padding: 8px 12px;
        text-align: center;
        margin-top: -15px;
        margin-bottom: 10px;
    }
    .metric-title {
        font-size: 0.8rem;
        color: #a4b0be;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    .metric-val {
        font-size: 1.85rem;
        font-weight: bold;
        color: #ffffff;
        line-height: 1.2;
    }
    .metric-sub {
        font-size: 0.85rem;
        color: #00d2d3;
        font-weight: 600;
    }
    .scale-table {
        font-size: 0.78rem;
        width: 100%;
        border-collapse: collapse;
    }
    .scale-table td, .scale-table th {
        padding: 4px 6px;
        border-bottom: 1px solid rgba(255,255,255,0.06);
    }
</style>
""", unsafe_allow_html=True)

# --- MODELOS CIENTÍFICOS RIGUROSOS ---

def calcular_tmrt(tg: float, ta: float, va: float, d: float = 0.15, emissivity: float = 0.95) -> float:
    """Calcula Tmrt según balance térmico de globo (ISO 7726)"""
    tg_k, ta_k = tg + 273.15, ta + 273.15
    h_forced = 1.1e8 * (va ** 0.6) / (emissivity * (d ** 0.4))
    h_natural = (0.25e8 / emissivity) * (abs(tg - ta) / d) ** 0.25
    h_conv = max(h_forced, h_natural)
    val = (tg_k ** 4) + h_conv * (tg - ta)
    return round((val ** 0.25) - 273.15, 2) if val > 0 else tg

def calcular_pmv_ppd_fanger(ta: float, tr: float, va: float, rh: float, met: float, clo: float):
    """Implementación rigurosa del modelo PMV/PPD de Fanger (ISO 7730)"""
    pa = (rh / 100.0) * 10.0 * 0.1 * math.exp(16.6536 - 4030.183 / (ta + 235.0))  # Presión vapor (Pa)
    m = met * 58.15  # W/m2
    w = 0.0          # Trabajo mecánico
    icl = 0.155 * clo
    fcl = 1.05 + 0.645 * icl if icl > 0.078 else 1.0 + 1.29 * icl
    
    # Iteración temperatura exterior de la ropa (tcl)
    tcla = ta + (35.5 - ta) / (3.5 * (6.45 * icl + 0.1))
    tcl = tcla
    for _ in range(150):
        hc = max(2.38 * abs(tcl - ta) ** 0.25, 12.1 * math.sqrt(max(va, 0.05)))
        tcl_new = 35.7 - 0.0275 * (m - w) - icl * (
            3.96e-8 * fcl * ((tcl + 273.15) ** 4 - (tr + 273.15) ** 4) + fcl * hc * (tcl - ta)
        )
        if abs(tcl_new - tcl) < 0.001:
            break
        tcl = (tcl + tcl_new) / 2.0

    hc = max(2.38 * abs(tcl - ta) ** 0.25, 12.1 * math.sqrt(max(va, 0.05)))
    hl1 = 3.05 * 0.001 * (5733.0 - 6.99 * (m - w) - pa)
    hl2 = 0.42 * ((m - w) - 58.15) if (m - w) > 58.15 else 0.0
    hl3 = 1.7e-5 * m * (5867.0 - pa)
    hl4 = 0.0014 * m * (34.0 - ta)
    hl5 = 3.96e-8 * fcl * ((tcl + 273.15) ** 4 - (tr + 273.15) ** 4)
    hl6 = fcl * hc * (tcl - ta)
    
    l = (m - w) - (hl1 + hl2 + hl3 + hl4 + hl5 + hl6)
    ts = 0.303 * math.exp(-0.036 * m) + 0.028
    pmv = ts * l
    ppd = 100.0 - 95.0 * math.exp(-0.03353 * (pmv ** 4) - 0.2179 * (pmv ** 2))
    return round(pmv, 2), round(min(ppd, 100.0), 1)

def calcular_utci_aprox(ta: float, tr: float, va: float, rh: float) -> float:
    """Calcula UTCI utilizando pythermalcomfort si está instalado o aproximación polinomial"""
    try:
        from pythermalcomfort.models import utci
        v10 = max(0.5, va * (math.log(10.0 / 0.01) / math.log(1.2 / 0.01)))
        res = utci(tdb=ta, tr=tr, v=v10, rh=rh, limit_inputs=False)
        return round(float(res.utci if hasattr(res, 'utci') else res), 1)
    except:
        # Aproximación polinomial estándar de apoyo
        d_tr = tr - ta
        pa_kpa = (rh / 100.0) * 0.61078 * math.exp((17.27 * ta) / (ta + 237.3))
        v10 = max(0.5, va * 1.44)
        utci_aprox = ta + (0.6075 * d_tr) - (0.0288 * d_tr * v10) - (0.0036 * ta * v10) + (0.09 * (pa_kpa - 1.0))
        return round(utci_aprox, 1)

# ==============================================================================
# BLOQUE 1: ENTRADA DE DATOS (DISEÑO TÉCNICO E INTUITIVO)
# ==============================================================================
st.markdown('<div class="block-header">1. Variables de Entrada | Monitorización Experimental y Humana</div>', unsafe_allow_html=True)

with st.container():
    c_env1, c_env2, c_env3, c_hum1, c_hum2 = st.columns([1.1, 1.1, 1.1, 1.2, 1.3])
    
    with c_env1:
        st.markdown("**Temperaturas Ambientales (°C)**")
        t_a = st.number_input("Aire / Bulbo Seco (Ta):", value=28.0, step=0.1, format="%.2f")
        t_g = st.number_input("Termómetro Globo (Tg):", value=34.5, step=0.1, format="%.2f")

    with c_env2:
        st.markdown("**Humedad y Bulbo Húmedo**")
        t_nw = st.number_input("Bulbo Húmedo Nat. (Tnw) [°C]:", value=23.0, step=0.1, format="%.2f")
        rh = st.slider("Humedad Relativa (RH) [%]:", min_value=5.0, max_value=100.0, value=55.0, step=1.0)

    with c_env3:
        st.markdown("**Dinámica de Aire y Sol**")
        v_a = st.number_input("Velocidad Aire (va) [m/s]:", value=0.4, min_value=0.05, max_value=20.0, step=0.05)
        con_sol = st.checkbox("Carga solar directa (WBGT exterior)", value=False)

    with c_hum1:
        st.markdown("**Factores Personales (PMV)**")
        met = st.number_input("Metabolismo (met):", value=1.2, min_value=0.7, max_value=4.0, step=0.1, help="Oficina=1.2, Ligero=1.6, Moderado=2.0, Pesado=3.0")
        clo = st.number_input("Aislamiento Ropa (clo):", value=0.6, min_value=0.0, max_value=2.5, step=0.05, help="Verano=0.5, Oficina=0.7, Invierno=1.0")

    with c_hum2:
        st.markdown("**Límites Normativos (WBGT)**")
        clase_iso = st.selectbox("Actividad laboral (ISO 7243):", [
            "Clase 1: Ligera (115-180 W)",
            "Clase 2: Moderada (180-300 W)",
            "Clase 3: Pesada (300-390 W)",
            "Clase 4: Muy pesada (>390 W)"
        ], index=1)
        aclimatado = st.toggle("Operario Aclimatado", value=True)

# --- PROCESAMIENTO CIENTÍFICO DE ÍNDICES ---
t_mrt = calcular_tmrt(t_g, t_a, v_a)
pmv, ppd = calcular_pmv_ppd_fanger(t_a, t_mrt, v_a, rh, met, clo)

# WBGT
if con_sol:
    wbgt = round(0.7 * t_nw + 0.2 * t_g + 0.1 * t_a, 2)
else:
    wbgt = round(0.7 * t_nw + 0.3 * t_g, 2)

limites_wbgt = {
    "Clase 1: Ligera (115-180 W)": 30.0 if aclimatado else 29.0,
    "Clase 2: Moderada (180-300 W)": 28.0 if aclimatado else 26.0,
    "Clase 3: Pesada (300-390 W)": 26.0 if aclimatado else 23.0,
    "Clase 4: Muy pesada (>390 W)": 25.0 if aclimatado else 20.0
}
limite_actual = limites_wbgt[clase_iso]

# UTCI
utci = calcular_utci_aprox(t_a, t_mrt, v_a, rh)

st.write("") # Espaciador

# ==============================================================================
# FILA DE SALIDAS: BLOQUE 2, BLOQUE 3 Y BLOQUE 4
# ==============================================================================
col_b2, col_b3, col_b4 = st.columns(3)

# ------------------------------------------------------------------------------
# BLOQUE 2: PMV / PPD (ISO 7730)
# ------------------------------------------------------------------------------
with col_b2:
    st.markdown('<div class="block-header">2. Confort Térmico: PMV / PPD</div>', unsafe_allow_html=True)
    
    # Gráfico Semicírculo
    fig_pmv = go.Figure(go.Indicator(
        mode="gauge",
        value=pmv,
        gauge={
            'axis': {'range': [-3, 3], 'tickwidth': 1, 'tickvals': [-3, -2, -1, 0, 1, 2, 3]},
            'bar': {'color': "#ffffff", 'thickness': 0.28},
            'steps': [
                {'range': [-3, -2], 'color': '#1e3799'},
                {'range': [-2, -0.7], 'color': '#4a69bd'},
                {'range': [-0.7, 0.7], 'color': '#2ed573'}, # Zona confort ISO
                {'range': [0.7, 2], 'color': '#eccc68'},
                {'range': [2, 3], 'color': '#eb4d4b'}
            ],
            'threshold': {'line': {'color': '#2ed573', 'width': 3}, 'thickness': 0.75, 'value': 0}
        }
    ))
    fig_pmv.update_layout(height=180, margin=dict(l=15, r=15, t=15, b=0), paper_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_pmv, use_container_width=True)

    # Valor numérico justo debajo
    st.markdown(f"""
    <div class="metric-box">
        <div class="metric-title">Índices de Confort ISO 7730</div>
        <div class="metric-val">PMV: {pmv:+.2f}</div>
        <div class="metric-sub">PPD: {ppd:.1f}% insatisfechos (Tmrt: {t_mrt:.1f} °C)</div>
    </div>
    """, unsafe_allow_html=True)

    # Escala normativa inferior
    st.markdown("""
    <table class="scale-table">
        <tr style="color:#a4b0be;"><th>PMV</th><th>Sensación</th><th>PPD [%]</th><th>Categoría ISO</th></tr>
        <tr><td>+2 a +3</td><td>Calor severo</td><td>> 75%</td><td>Fuera norma</td></tr>
        <tr><td>+0.7 a +2</td><td>Lig. Cálido</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr style="color:#2ed573; font-weight:bold;"><td>-0.5 a +0.5</td><td>Neutro / Confort</td><td>< 10%</td><td>Cat. B (Estándar)</td></tr>
        <tr><td>-2 a -0.7</td><td>Lig. Frío</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr><td>-3 a -2</td><td>Frío severo</td><td>> 75%</td><td>Fuera norma</td></tr>
    </table>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 3: WBGT (ISO 7243)
# ------------------------------------------------------------------------------
with col_b3:
    st.markdown('<div class="block-header">3. Estrés Térmico: WBGT</div>', unsafe_allow_html=True)

    # Gráfico Semicírculo
    fig_wbgt = go.Figure(go.Indicator(
        mode="gauge",
        value=wbgt,
        gauge={
            'axis': {'range': [10, 40], 'tickwidth': 1, 'tickvals': [10, 20, limite_actual, 40]},
            'bar': {'color': "#ffffff", 'thickness': 0.28},
            'steps': [
                {'range': [10, limite_actual - 2], 'color': '#2ed573'},
                {'range': [limite_actual - 2, limite_actual], 'color': '#ffa502'},
                {'range': [limite_actual, 40], 'color': '#ff4757'}
            ],
            'threshold': {'line': {'color': '#ffffff', 'width': 3}, 'thickness': 0.8, 'value': limite_actual}
        }
    ))
    fig_wbgt.update_layout(height=180, margin=dict(l=15, r=15, t=15, b=0), paper_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_wbgt, use_container_width=True)

    # Valor numérico justo debajo
    delta_w = round(wbgt - limite_actual, 2)
    estado_w = f"+{delta_w} °C sobre límite" if delta_w > 0 else f"{abs(delta_w)} °C bajo límite seguro"
    st.markdown(f"""
    <div class="metric-box">
        <div class="metric-title">Índice WBGT Ocupacional (ISO 7243)</div>
        <div class="metric-val">{wbgt:.2f} °C</div>
        <div class="metric-sub">Límite Normativo: {limite_actual:.1f} °C ({estado_w})</div>
    </div>
    """, unsafe_allow_html=True)

    # Escala normativa inferior
    st.markdown(f"""
    <table class="scale-table">
        <tr style="color:#a4b0be;"><th>WBGT (°C)</th><th>Criterio Preventivo (Límite: {limite_actual:.1f} °C)</th><th>Acción</th></tr>
        <tr style="color:#2ed573;"><td>&lt; {limite_actual - 2:.1f}</td><td>Zona Segura / Conforme</td><td>Jornada continua normal</td></tr>
        <tr style="color:#ffa502;"><td>{limite_actual - 2:.1f} - {limite_actual:.1f}</td><td>Alerta de Sobrecarga</td><td>Pausas / Hidratación</td></tr>
        <tr style="color:#ff4757; font-weight:bold;"><td>&gt; {limite_actual:.1f}</td><td>Excede Límite Permisible</td><td>Descansos 50% / Suspender</td></tr>
    </table>
    <div style="font-size:0.73rem; color:#888; margin-top:8px;">*Límites calculados en base a tasa metabólica y aclimatación.</div>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 4: UTCI (COST ACTION 730 / WMO)
# ------------------------------------------------------------------------------
with col_b4:
    st.markdown('<div class="block-header">4. Bioclima Universal: UTCI</div>', unsafe_allow_html=True)

    # Gráfico Semicírculo
    fig_utci = go.Figure(go.Indicator(
        mode="gauge",
        value=utci,
        gauge={
            'axis': {'range': [-15, 50], 'tickwidth': 1, 'tickvals': [-13, 0, 9, 26, 32, 38, 46]},
            'bar': {'color': "#ffffff", 'thickness': 0.28},
            'steps': [
                {'range': [-15, 0], 'color': '#1e3799'},
                {'range': [0, 9], 'color': '#00a8ff'},
                {'range': [9, 26], 'color': '#2ed573'},   # Confort
                {'range': [26, 32], 'color': '#fbc531'},  # Moderado
                {'range': [32, 38], 'color': '#e67e22'},  # Fuerte
                {'range': [38, 46], 'color': '#e84118'},  # Muy fuerte
                {'range': [46, 50], 'color': '#8c7ae6'}   # Extremo
            ],
            'threshold': {'line': {'color': '#ffffff', 'width': 3}, 'thickness': 0.8, 'value': utci}
        }
    ))
    fig_utci.update_layout(height=180, margin=dict(l=15, r=15, t=15, b=0), paper_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_utci, use_container_width=True)

    # Valor numérico justo debajo
    if utci > 46: cat_u = "Estrés Extremo Calor"
    elif utci > 38: cat_u = "Estrés Muy Fuerte Calor"
    elif utci > 32: cat_u = "Estrés Fuerte Calor"
    elif utci > 26: cat_u = "Estrés Moderado Calor"
    elif utci >= 9: cat_u = "Sin Estrés (Confort Térmico)"
    elif utci >= 0: cat_u = "Estrés Ligero Frío"
    else: cat_u = "Estrés Fuerte Frío"

    st.markdown(f"""
    <div class="metric-box">
        <div class="metric-title">Índice Térmico Universal (UTCI)</div>
        <div class="metric-val">{utci:.1f} °C</div>
        <div class="metric-sub">{cat_u}</div>
    </div>
    """, unsafe_allow_html=True)

    # Escala normativa inferior
    st.markdown("""
    <table class="scale-table">
        <tr style="color:#a4b0be;"><th>Rango UTCI</th><th>Categoría de Estrés</th><th>Criterio Fisiológico</th></tr>
        <tr style="color:#8c7ae6;"><td>&gt; +46 °C</td><td>Extremo por calor</td><td>Peligro térmico inminente</td></tr>
        <tr style="color:#e84118;"><td>+38 a +46 °C</td><td>Muy fuerte por calor</td><td>Sudoración y taquicardia severa</td></tr>
        <tr style="color:#e67e22;"><td>+32 a +38 °C</td><td>Fuerte por calor</td><td>Sobrecarga cardiovascular</td></tr>
        <tr style="color:#fbc531;"><td>+26 a +32 °C</td><td>Moderado por calor</td><td>Bochorno / Sudoración leve</td></tr>
        <tr style="color:#2ed573; font-weight:bold;"><td>+9 a +26 °C</td><td>Confort Térmico</td><td>Equilibrio fisiológico óptimo</td></tr>
        <tr style="color:#00a8ff;"><td>&lt; +9 °C</td><td>Estrés por frío</td><td>Vasoconstricción / Tiritona</td></tr>
    </table>
    """, unsafe_allow_html=True)
