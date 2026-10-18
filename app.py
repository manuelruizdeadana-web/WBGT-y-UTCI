import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math

# Configuración de página
st.set_page_config(
    page_title="Monitor Bioclimático de Alta Precisión | ISO & WMO",
    page_icon="⚡",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- CSS TÉCNICO PROFESIONAL: TEMA SALA DE CONTROL / INSTRUMENTACIÓN ---
st.markdown("""
<style>
    /* Forzar tipografía técnica Arial Narrow / Sans condensada */
    * {
        font-family: 'Arial Narrow', 'Franklin Gothic Medium', 'Helvetica Neue', sans-serif !important;
    }
    
    /* Fondo global técnico de alto contraste */
    .stApp {
        background-color: #0b0f19 !important;
        color: #f1f5f9 !important;
    }

    /* Cabeceras de bloques principales con estilo industrial */
    .block-card {
        background: #111827;
        border: 1px solid #1f293d;
        border-radius: 8px;
        padding: 14px;
        margin-bottom: 12px;
        box-shadow: 0 4px 15px rgba(0, 0, 0, 0.4);
    }
    
    .block-header {
        background: linear-gradient(90deg, #1e293b, #0f172a);
        border-left: 4px solid #38bdf8;
        padding: 8px 14px;
        font-size: 0.95rem;
        font-weight: 700;
        letter-spacing: 1px;
        text-transform: uppercase;
        color: #f8fafc;
        border-radius: 4px;
        margin-bottom: 12px;
        display: flex;
        align-items: center;
        gap: 8px;
    }

    /* Display digital para números principales */
    .display-box {
        background: #070a12;
        border: 1px solid #1e293b;
        border-radius: 6px;
        padding: 10px;
        text-align: center;
        margin-top: -10px;
        margin-bottom: 12px;
        box-shadow: inset 0 2px 8px rgba(0,0,0,0.6);
    }
    .display-label {
        font-size: 0.75rem;
        color: #94a3b8;
        letter-spacing: 1.5px;
        text-transform: uppercase;
        font-weight: 600;
        margin-bottom: 2px;
    }
    .display-value {
        font-size: 2.2rem;
        font-weight: 800;
        line-height: 1.1;
        letter-spacing: -0.5px;
    }
    .display-badge {
        display: inline-block;
        padding: 3px 10px;
        border-radius: 20px;
        font-size: 0.78rem;
        font-weight: 700;
        letter-spacing: 0.5px;
        margin-top: 5px;
        text-transform: uppercase;
    }

    /* Tablas normativas técnicas */
    .tech-table {
        width: 100%;
        border-collapse: separate;
        border-spacing: 0;
        font-size: 0.8rem;
        margin-top: 5px;
    }
    .tech-table th {
        background-color: #1e293b;
        color: #94a3b8;
        font-weight: 700;
        text-transform: uppercase;
        font-size: 0.7rem;
        letter-spacing: 0.5px;
        padding: 5px 8px;
        border-bottom: 1px solid #334155;
    }
    .tech-table td {
        padding: 6px 8px;
        border-bottom: 1px solid #1e293b;
        color: #cbd5e1;
    }
    /* Resaltado dinámico de la fila activa */
    .active-row {
        background-color: rgba(56, 189, 248, 0.18) !important;
        border-left: 3px solid #38bdf8 !important;
        font-weight: bold !important;
        color: #ffffff !important;
    }

    /* Personalización de inputs en Streamlit */
    div[data-baseweb="input"], div[data-baseweb="select"] {
        background-color: #1e293b !important;
        border-color: #334155 !important;
        color: #ffffff !important;
    }
</style>
""", unsafe_allow_html=True)

# --- MODELOS TERMOFISIOLÓGICOS RIGUROSOS ---

def calcular_tmrt(tg: float, ta: float, va: float, d: float = 0.15, emissivity: float = 0.95) -> float:
    """ISO 7726: Balance radiativo y convectivo del termómetro de globo"""
    tg_k, ta_k = tg + 273.15, ta + 273.15
    h_forced = 1.1e8 * (va ** 0.6) / (emissivity * (d ** 0.4))
    h_natural = (0.25e8 / emissivity) * (abs(tg - ta) / d) ** 0.25
    h_conv = max(h_forced, h_natural)
    val = (tg_k ** 4) + h_conv * (tg - ta)
    return round((val ** 0.25) - 273.15, 2) if val > 0 else tg

def calcular_pmv_ppd(ta: float, tr: float, va: float, rh: float, met: float, clo: float):
    """ISO 7730: Ecuación analítica de confort de Fanger"""
    pa = (rh / 100.0) * 10.0 * 0.1 * math.exp(16.6536 - 4030.183 / (ta + 235.0))
    m = met * 58.15
    w = 0.0
    icl = 0.155 * clo
    fcl = 1.05 + 0.645 * icl if icl > 0.078 else 1.0 + 1.29 * icl
    
    tcl = ta + (35.5 - ta) / (3.5 * (6.45 * icl + 0.1))
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
    
    pmv = (0.303 * math.exp(-0.036 * m) + 0.028) * (m - w - (hl1 + hl2 + hl3 + hl4 + hl5 + hl6))
    ppd = 100.0 - 95.0 * math.exp(-0.03353 * (pmv ** 4) - 0.2179 * (pmv ** 2))
    return round(pmv, 2), round(min(ppd, 100.0), 1)

def calcular_utci(ta: float, tr: float, va: float, rh: float) -> float:
    """COST Action 730 / WMO: Índice bioclimático universal"""
    try:
        from pythermalcomfort.models import utci as ptc_utci
        v10 = max(0.5, va * (math.log(10.0 / 0.01) / math.log(1.2 / 0.01)))
        res = ptc_utci(tdb=ta, tr=tr, v=v10, rh=rh, limit_inputs=False)
        return round(float(res.utci if hasattr(res, 'utci') else res), 1)
    except:
        # Modelo aproximado polinomial
        d_tr = tr - ta
        pa_kpa = (rh / 100.0) * 0.61078 * math.exp((17.27 * ta) / (ta + 237.3))
        v10 = max(0.5, va * 1.44)
        val = ta + (0.6075 * d_tr) - (0.0288 * d_tr * v10) - (0.0036 * ta * v10) + (0.09 * (pa_kpa - 1.0))
        return round(val, 1)

# ==============================================================================
# BLOQUE 1: PANEL DE ENTRADA DE DATOS DE ALTO RENDIMIENTO
# ==============================================================================
st.markdown('<div class="block-header">🎛️ 1. VARIABLES DE ENTRADA | Monitorización Experimental y Parámetros Operativos</div>', unsafe_allow_html=True)

with st.container():
    # Fila rápida de presets de prueba
    col_preset, col_info = st.columns([1.5, 3.5])
    with col_preset:
        preset = st.selectbox(
            "⚡ Cargar Escenario Predefinido:",
            ("Personalizado", "Oficina Confortable (Verano)", "Nave Industrial Severa", "Exterior Carga Solar Extrema")
        )

    # Carga de valores por defecto según preset
    if preset == "Oficina Confortable (Verano)":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 24.0, 24.5, 17.5, 50.0, 0.15, 1.2, 0.5, False
    elif preset == "Nave Industrial Severa":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 32.0, 39.0, 25.0, 60.0, 0.30, 2.0, 0.8, False
    elif preset == "Exterior Carga Solar Extrema":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 36.5, 48.0, 28.0, 45.0, 1.80, 2.4, 0.6, True
    else:
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 28.0, 34.5, 23.0, 55.0, 0.40, 1.2, 0.6, False

    # Distribución en 5 columnas con bordes y fondo oscuro
    c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1.1, 1.2])

    with c1:
        st.caption("🌡️ **TEMPERATURAS (°C)**")
        t_a = st.number_input("Aire / Seco (Ta):", value=d_ta, step=0.1, format="%.2f")
        t_g = st.number_input("Globo (Tg):", value=d_tg, step=0.1, format="%.2f")

    with c2:
        st.caption("💧 **HUMEDADES**")
        t_nw = st.number_input("Bulbo Húmedo (Tnw):", value=d_tnw, step=0.1, format="%.2f")
        rh = st.slider("Humedad Relativa (RH %):", 5.0, 100.0, d_rh, 1.0)

    with c3:
        st.caption("💨 **DINÁMICA & SOL**")
        v_a = st.number_input("Velocidad Aire (va, m/s):", value=d_va, min_value=0.05, max_value=25.0, step=0.05)
        con_sol = st.checkbox("Radiación Solar Directa", value=d_sol)

    with c4:
        st.caption("👤 **FACTORES HUMANOS (PMV)**")
        met = st.number_input("Metabolismo (met):", value=d_met, min_value=0.7, max_value=4.0, step=0.1)
        clo = st.number_input("Aislamiento Ropa (clo):", value=d_clo, min_value=0.0, max_value=2.5, step=0.05)

    with c5:
        st.caption("📋 **NORMATIVA (WBGT ISO 7243)**")
        clase_iso = st.selectbox("Actividad metabólica:", [
            "Clase 1: Ligera (115-180 W)",
            "Clase 2: Moderada (180-300 W)",
            "Clase 3: Pesada (300-390 W)",
            "Clase 4: Muy pesada (>390 W)"
        ], index=1)
        aclimatado = st.toggle("Operario Aclimatado", value=True)

# --- PROCESAMIENTO MATEMÁTICO ---
t_mrt = calcular_tmrt(t_g, t_a, v_a)
pmv, ppd = calcular_pmv_ppd(t_a, t_mrt, v_a, rh, met, clo)

# WBGT
wbgt = round(0.7 * t_nw + 0.2 * t_g + 0.1 * t_a, 2) if con_sol else round(0.7 * t_nw + 0.3 * t_g, 2)
limites_wbgt = {
    "Clase 1: Ligera (115-180 W)": 30.0 if aclimatado else 29.0,
    "Clase 2: Moderada (180-300 W)": 28.0 if aclimatado else 26.0,
    "Clase 3: Pesada (300-390 W)": 26.0 if aclimatado else 23.0,
    "Clase 4: Muy pesada (>390 W)": 25.0 if aclimatado else 20.0
}
limite_wbgt = limites_wbgt[clase_iso]
delta_wbgt = round(wbgt - limite_wbgt, 2)

# UTCI
utci = calcular_utci(t_a, t_mrt, v_a, rh)

st.write("") # Espaciador vertical

# ==============================================================================
# FILA DE RESULTADOS: 3 COLUMNAS INDEPENDIENTES (BLOQUES 2, 3 Y 4)
# ==============================================================================
col_b2, col_b3, col_b4 = st.columns(3)

# ------------------------------------------------------------------------------
# BLOQUE 2: PMV / PPD (ISO 7730)
# ------------------------------------------------------------------------------
with col_b2:
    st.markdown('<div class="block-header">🏢 2. CONFORT TÉRMICO | PMV & PPD (ISO 7730)</div>', unsafe_allow_html=True)
    
    # Semicírculo Plotly
    fig_pmv = go.Figure(go.Indicator(
        mode="gauge",
        value=pmv,
        gauge={
            'axis': {'range': [-3.0, 3.0], 'tickwidth': 1.5, 'tickcolor': "#94a3b8", 'tickvals': [-3, -2, -1, 0, 1, 2, 3]},
            'bar': {'color': "#ffffff", 'thickness': 0.32},
            'bgcolor': "#070a12",
            'steps': [
                {'range': [-3.0, -2.0], 'color': '#1d4ed8'},  # Azul intenso
                {'range': [-2.0, -0.7], 'color': '#0284c7'},  # Azul cielo
                {'range': [-0.7, 0.7],  'color': '#10b981'},  # Verde confort ISO
                {'range': [0.7, 2.0],   'color': '#f59e0b'},  # Naranja
                {'range': [2.0, 3.0],   'color': '#ef4444'}   # Rojo
            ],
            'threshold': {'line': {'color': '#38bdf8', 'width': 4}, 'thickness': 0.8, 'value': pmv}
        }
    ))
    fig_pmv.update_layout(height=170, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_pmv, use_container_width=True)

    # Color del badge y diagnóstico
    if abs(pmv) <= 0.5:
        pmv_color = "#10b981"; pmv_bg = "rgba(16,185,129,0.2)"; pmv_txt = "CONFORT ÓPTIMO (CAT. A/B)"
    elif abs(pmv) <= 0.7:
        pmv_color = "#38bdf8"; pmv_bg = "rgba(56,189,248,0.2)"; pmv_txt = "CONFORT ACEPTABLE (CAT. C)"
    elif pmv > 0.7:
        pmv_color = "#ef4444"; pmv_bg = "rgba(239,68,68,0.2)"; pmv_txt = "SOBRECALENTAMIENTO"
    else:
        pmv_color = "#0284c7"; pmv_bg = "rgba(2,132,199,0.2)"; pmv_txt = "SUBENFRIAMIENTO"

    # Display digital LCD
    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Voto Medio Estimado / Insatisfechos</div>
        <div class="display-value" style="color: {pmv_color};">PMV: {pmv:+.2f}</div>
        <div style="font-size:0.95rem; color:#f8fafc; font-weight:600; margin-top:2px;">PPD: {ppd:.1f}% insatisfechos</div>
        <div class="display-badge" style="background:{pmv_bg}; color:{pmv_color};">{pmv_txt}</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:4px;">Tmrt calculada (ISO 7726): <b>{t_mrt:.1f} °C</b></div>
    </div>
    """, unsafe_allow_html=True)

    # Escala con fila activa resaltada
    r1 = "active-row" if pmv > 2.0 else ""
    r2 = "active-row" if 0.7 < pmv <= 2.0 else ""
    r3 = "active-row" if -0.5 <= pmv <= 0.5 else ""
    r4 = "active-row" if -2.0 <= pmv < -0.7 else ""
    r5 = "active-row" if pmv < -2.0 else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>PMV</th><th>Sensación</th><th>PPD</th><th>Norma</th></tr>
        <tr class="{r1}"><td>+2 a +3</td><td style="color:#ef4444;">Calor severo</td><td>&gt; 75%</td><td>Fuera norma</td></tr>
        <tr class="{r2}"><td>+0.7 a +2</td><td style="color:#f59e0b;">Lig. Cálido</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr class="{r3}"><td style="color:#10b981;">-0.5 a +0.5</td><td style="color:#10b981;">Neutro / Confort</td><td>&lt; 10%</td><td>Cat. B (Estándar)</td></tr>
        <tr class="{r4}"><td>-2 a -0.7</td><td style="color:#0284c7;">Lig. Frío</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr class="{r5}"><td>-3 a -2</td><td style="color:#1d4ed8;">Frío severo</td><td>&gt; 75%</td><td>Fuera norma</td></tr>
    </table>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 3: WBGT (ISO 7243)
# ------------------------------------------------------------------------------
with col_b3:
    st.markdown('<div class="block-header">🏭 3. ESTRÉS TÉRMICO LABORAL | WBGT (ISO 7243)</div>', unsafe_allow_html=True)

    # Semicírculo Plotly
    fig_wbgt = go.Figure(go.Indicator(
        mode="gauge",
        value=wbgt,
        gauge={
            'axis': {'range': [10.0, 42.0], 'tickwidth': 1.5, 'tickcolor': "#94a3b8", 'tickvals': [10, 20, limite_wbgt, 40]},
            'bar': {'color': "#ffffff", 'thickness': 0.32},
            'bgcolor': "#070a12",
            'steps': [
                {'range': [10.0, limite_wbgt - 2.0], 'color': '#10b981'},  # Zona Segura
                {'range': [limite_wbgt - 2.0, limite_wbgt], 'color': '#f59e0b'},  # Zona Alerta
                {'range': [limite_wbgt, 42.0], 'color': '#ef4444'}   # Peligro
            ],
            'threshold': {'line': {'color': '#f8fafc', 'width': 4}, 'thickness': 0.8, 'value': limite_wbgt}
        }
    ))
    fig_wbgt.update_layout(height=170, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_wbgt, use_container_width=True)

    # Color del badge y diagnóstico
    if wbgt <= (limite_wbgt - 2.0):
        wbgt_color = "#10b981"; wbgt_bg = "rgba(16,185,129,0.2)"; wbgt_txt = "ZONA SEGURA CONFORME"
        diff_txt = f"{abs(delta_wbgt):.2f} °C por debajo del límite seguro"
    elif wbgt <= limite_wbgt:
        wbgt_color = "#f59e0b"; wbgt_bg = "rgba(245,158,11,0.2)"; wbgt_txt = "ALERTA PREVENTIVA"
        diff_txt = f"A solo {abs(delta_wbgt):.2f} °C de superar el límite"
    else:
        wbgt_color = "#ef4444"; wbgt_bg = "rgba(239,68,68,0.2)"; wbgt_txt = "EXCEDE LÍMITE LABORAL"
        diff_txt = f"+{delta_wbgt:.2f} °C de sobrecarga térmica"

    # Display digital LCD
    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Índice WBGT Ocupacional Medido</div>
        <div class="display-value" style="color: {wbgt_color};">{wbgt:.2f} °C</div>
        <div style="font-size:0.95rem; color:#f8fafc; font-weight:600; margin-top:2px;">Límite Normativo: {limite_wbgt:.1f} °C</div>
        <div class="display-badge" style="background:{wbgt_bg}; color:{wbgt_color};">{wbgt_txt}</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:4px;">{diff_txt}</div>
    </div>
    """, unsafe_allow_html=True)

    # Escala con fila activa resaltada
    w_r1 = "active-row" if wbgt <= (limite_wbgt - 2.0) else ""
    w_r2 = "active-row" if (limite_wbgt - 2.0) < wbgt <= limite_wbgt else ""
    w_r3 = "active-row" if wbgt > limite_wbgt else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>Rango WBGT</th><th>Evaluación</th><th>Acción Preventiva</th></tr>
        <tr class="{w_r1}"><td style="color:#10b981;">&lt; {limite_wbgt - 2:.1f} °C</td><td style="color:#10b981;">Zona Conforme</td><td>Jornada continua / Normal</td></tr>
        <tr class="{w_r2}"><td style="color:#f59e0b;">{limite_wbgt - 2:.1f} a {limite_wbgt:.1f} °C</td><td style="color:#f59e0b;">Alerta Térmica</td><td>Pausas 15 min / Hidratación</td></tr>
        <tr class="{w_r3}"><td style="color:#ef4444;">&gt; {limite_wbgt:.1f} °C</td><td style="color:#ef4444;">Sobrecarga Severa</td><td>Pausas 50% / Suspender</td></tr>
    </table>
    <div style="font-size:0.7rem; color:#64748b; margin-top:5px;">*Límite dinámico según tasa metabólica seleccionada.</div>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 4: UTCI (COST ACTION 730 / WMO)
# ------------------------------------------------------------------------------
with col_b4:
    st.markdown('<div class="block-header">🌍 4. BIOCLIMA EXTERIOR | UTCI (COST 730 / WMO)</div>', unsafe_allow_html=True)

    # Semicírculo Plotly
    fig_utci = go.Figure(go.Indicator(
        mode="gauge",
        value=utci,
        gauge={
            'axis': {'range': [-15.0, 52.0], 'tickwidth': 1.5, 'tickcolor': "#94a3b8", 'tickvals': [-13, 0, 9, 26, 32, 38, 46]},
            'bar': {'color': "#ffffff", 'thickness': 0.32},
            'bgcolor': "#070a12",
            'steps': [
                {'range': [-15.0, 0.0], 'color': '#1d4ed8'},  # Azul frío
                {'range': [0.0, 9.0],   'color': '#0284c7'},  # Frío ligero
                {'range': [9.0, 26.0],  'color': '#10b981'},  # Confort universal
                {'range': [26.0, 32.0], 'color': '#eab308'},  # Calor moderado
                {'range': [32.0, 38.0], 'color': '#f97316'},  # Calor fuerte
                {'range': [38.0, 46.0], 'color': '#ef4444'},  # Muy fuerte
                {'range': [46.0, 52.0], 'color': '#a855f7'}   # Extremo
            ],
            'threshold': {'line': {'color': '#f8fafc', 'width': 4}, 'thickness': 0.8, 'value': utci}
        }
    ))
    fig_utci.update_layout(height=170, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_utci, use_container_width=True)

    # Diagnóstico UTCI
    if utci > 46:
        u_col = "#a855f7"; u_bg = "rgba(168,85,247,0.2)"; u_diag = "ESTRÉS EXTREMO CALOR"
    elif utci > 38:
        u_col = "#ef4444"; u_bg = "rgba(239,68,68,0.2)"; u_diag = "ESTRÉS MUY FUERTE CALOR"
    elif utci > 32:
        u_col = "#f97316"; u_bg = "rgba(249,115,22,0.2)"; u_diag = "ESTRÉS FUERTE CALOR"
    elif utci > 26:
        u_col = "#eab308"; u_bg = "rgba(234,179,8,0.2)"; u_diag = "ESTRÉS MODERADO CALOR"
    elif utci >= 9:
        u_col = "#10b981"; u_bg = "rgba(16,185,129,0.2)"; u_diag = "ZONA DE CONFORT TÉRMICO"
    elif utci >= 0:
        u_col = "#0284c7"; u_bg = "rgba(2,132,199,0.2)"; u_diag = "ESTRÉS LIGERO FRÍO"
    else:
        u_col = "#1d4ed8"; u_bg = "rgba(29,78,216,0.2)"; u_diag = "ESTRÉS FUERTE FRÍO"

    # Display digital LCD
    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Temperatura Climatológica Equivalente</div>
        <div class="display-value" style="color: {u_col};">{utci:+.1f} °C</div>
        <div style="font-size:0.95rem; color:#f8fafc; font-weight:600; margin-top:2px;">{u_diag}</div>
        <div class="display-badge" style="background:{u_bg}; color:{u_col};">ESCALA INTERNACIONAL WMO</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:4px;">Viento a 10m normalizado: <b>{max(0.5, v_a * 1.44):.1f} m/s</b></div>
    </div>
    """, unsafe_allow_html=True)

    # Escala con fila activa resaltada
    u1 = "active-row" if utci > 46 else ""
    u2 = "active-row" if 38 < utci <= 46 else ""
    u3 = "active-row" if 32 < utci <= 38 else ""
    u4 = "active-row" if 26 < utci <= 32 else ""
    u5 = "active-row" if 9 <= utci <= 26 else ""
    u6 = "active-row" if utci < 9 else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>Rango UTCI</th><th>Sensación Fisiológica</th><th>Respuesta Humana</th></tr>
        <tr class="{u1}"><td style="color:#a855f7;">&gt; +46 °C</td><td style="color:#a855f7;">Calor Extremo</td><td>Peligro térmico inminente</td></tr>
        <tr class="{u2}"><td style="color:#ef4444;">+38 a +46 °C</td><td style="color:#ef4444;">Calor Muy Fuerte</td><td>Sudoración masiva / Taquicardia</td></tr>
        <tr class="{u3}"><td style="color:#f97316;">+32 a +38 °C</td><td style="color:#f97316;">Calor Fuerte</td><td>Sobrecarga circulatoria</td></tr>
        <tr class="{u4}"><td style="color:#eab308;">+26 a +32 °C</td><td style="color:#eab308;">Calor Moderado</td><td>Bochorno / Sudoración</td></tr>
        <tr class="{u5}"><td style="color:#10b981;">+9 a +26 °C</td><td style="color:#10b981;">Confort Térmico</td><td>Equilibrio óptimo</td></tr>
        <tr class="{u6}"><td style="color:#0284c7;">&lt; +9 °C</td><td style="color:#0284c7;">Estrés por Frío</td><td>Vasoconstricción periférica</td></tr>
    </table>
    """, unsafe_allow_html=True)
