import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math

# Configuración de página
st.set_page_config(
    page_title="Monitor Bioclimático de Alta Precisión | Clean Tech",
    page_icon="🌡️",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# --- CSS DE ALTO CONTRASTE (TEMA CLARO / LABORATORIO) ---
st.markdown("""
<style>
    /* Tipografía técnica Arial Narrow / Sans-Serif condensada */
    * {
        font-family: 'Arial Narrow', 'Franklin Gothic Medium', 'Segoe UI', sans-serif !important;
    }
    
    /* Fondo general claro de alto contraste */
    .stApp {
        background-color: #f1f5f9 !important;
        color: #0f172a !important;
    }

    /* Tarjetas y bloques principales */
    .block-card {
        background: #ffffff;
        border: 1px solid #cbd5e1;
        border-radius: 8px;
        padding: 14px;
        margin-bottom: 12px;
        box-shadow: 0 2px 6px rgba(15, 23, 42, 0.06);
    }
    
    /* Encabezados de bloques estilo instrumentación */
    .block-header {
        background: #0f172a;
        border-left: 5px solid #0284c7;
        padding: 8px 14px;
        font-size: 0.92rem;
        font-weight: 700;
        letter-spacing: 0.8px;
        text-transform: uppercase;
        color: #ffffff !important;
        border-radius: 4px;
        margin-bottom: 10px;
    }

    /* Display digital para números principales */
    .display-box {
        background: #ffffff;
        border: 1.5px solid #e2e8f0;
        border-radius: 6px;
        padding: 10px;
        text-align: center;
        margin-top: -10px;
        margin-bottom: 12px;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.04);
    }
    .display-label {
        font-size: 0.76rem;
        color: #475569;
        letter-spacing: 1.2px;
        text-transform: uppercase;
        font-weight: 700;
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
        padding: 3px 12px;
        border-radius: 16px;
        font-size: 0.78rem;
        font-weight: 800;
        letter-spacing: 0.5px;
        margin-top: 6px;
        text-transform: uppercase;
    }

    /* Tablas técnicas de alta legibilidad */
    .tech-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.82rem;
        margin-top: 5px;
        background: #ffffff;
    }
    .tech-table th {
        background-color: #e2e8f0;
        color: #1e293b;
        font-weight: 800;
        text-transform: uppercase;
        font-size: 0.72rem;
        letter-spacing: 0.5px;
        padding: 6px 8px;
        border-bottom: 2px solid #cbd5e1;
        text-align: left;
    }
    .tech-table td {
        padding: 6px 8px;
        border-bottom: 1px solid #f1f5f9;
        color: #1e293b;
    }
    
    /* Resaltado dinámico e inequívoco de la fila activa */
    .active-row {
        background-color: #e0f2fe !important;
        border-left: 4px solid #0284c7 !important;
        font-weight: 800 !important;
        color: #0369a1 !important;
    }

    /* Etiquetas de entradas */
    .input-label {
        color: #0f172a;
        font-weight: 700;
        font-size: 0.85rem;
        text-transform: uppercase;
        margin-bottom: 4px;
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
        d_tr = tr - ta
        pa_kpa = (rh / 100.0) * 0.61078 * math.exp((17.27 * ta) / (ta + 237.3))
        v10 = max(0.5, va * 1.44)
        val = ta + (0.6075 * d_tr) - (0.0288 * d_tr * v10) - (0.0036 * ta * v10) + (0.09 * (pa_kpa - 1.0))
        return round(val, 1)

# ==============================================================================
# BLOQUE 1: PANEL DE ENTRADA DE DATOS (FONDO CLARO Y ALTO CONTRASTE)
# ==============================================================================
st.markdown('<div class="block-header">🎛️ 1. VARIABLES DE ENTRADA | Monitorización Experimental y Parámetros Operativos</div>', unsafe_allow_html=True)

with st.container():
    col_preset, col_pad = st.columns([1.6, 3.4])
    with col_preset:
        preset = st.selectbox(
            "Cargar Escenario Predefinido:",
            ("Personalizado", "Oficina Confortable (Verano)", "Nave Industrial Severa", "Exterior Carga Solar Extrema")
        )

    if preset == "Oficina Confortable (Verano)":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 24.0, 24.5, 17.5, 50.0, 0.15, 1.2, 0.5, False
    elif preset == "Nave Industrial Severa":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 32.0, 39.0, 25.0, 60.0, 0.30, 2.0, 0.8, False
    elif preset == "Exterior Carga Solar Extrema":
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 36.5, 48.0, 28.0, 45.0, 1.80, 2.4, 0.6, True
    else:
        d_ta, d_tg, d_tnw, d_rh, d_va, d_met, d_clo, d_sol = 28.0, 34.5, 23.0, 55.0, 0.40, 1.2, 0.6, False

    c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1.1, 1.2])

    with c1:
        st.markdown('<div class="input-label">🌡️ Temperaturas (°C)</div>', unsafe_allow_html=True)
        t_a = st.number_input("Aire / Seco (Ta):", value=d_ta, step=0.1, format="%.2f")
        t_g = st.number_input("Globo (Tg):", value=d_tg, step=0.1, format="%.2f")

    with c2:
        st.markdown('<div class="input-label">💧 Humedades</div>', unsafe_allow_html=True)
        t_nw = st.number_input("Bulbo Húmedo (Tnw) [°C]:", value=d_tnw, step=0.1, format="%.2f")
        rh = st.slider("Humedad Relativa (RH %):", 5.0, 100.0, d_rh, 1.0)

    with c3:
        st.markdown('<div class="input-label">💨 Dinámica & Sol</div>', unsafe_allow_html=True)
        v_a = st.number_input("Velocidad Aire (va, m/s):", value=d_va, min_value=0.05, max_value=25.0, step=0.05)
        con_sol = st.checkbox("Radiación Solar Directa", value=d_sol)

    with c4:
        st.markdown('<div class="input-label">👤 Factores Humanos</div>', unsafe_allow_html=True)
        met = st.number_input("Metabolismo (met):", value=d_met, min_value=0.7, max_value=4.0, step=0.1)
        clo = st.number_input("Aislamiento Ropa (clo):", value=d_clo, min_value=0.0, max_value=2.5, step=0.05)

    with c5:
        st.markdown('<div class="input-label">📋 Normativa WBGT</div>', unsafe_allow_html=True)
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
# FILA DE RESULTADOS: BLOQUE 2, BLOQUE 3 Y BLOQUE 4
# ==============================================================================
col_b2, col_b3, col_b4 = st.columns(3)

# ------------------------------------------------------------------------------
# BLOQUE 2: PMV / PPD (ISO 7730)
# ------------------------------------------------------------------------------
with col_b2:
    st.markdown('<div class="block-header">🏢 2. CONFORT TÉRMICO | PMV & PPD (ISO 7730)</div>', unsafe_allow_html=True)
    
    # Semicírculo Plotly con aguja negro grafito de alto contraste
    fig_pmv = go.Figure(go.Indicator(
        mode="gauge",
        value=pmv,
        gauge={
            'axis': {'range': [-3.0, 3.0], 'tickwidth': 1.5, 'tickcolor': "#334155", 'tickfont': {'color': '#0f172a', 'size': 11, 'family': 'Arial Narrow'}, 'tickvals': [-3, -2, -1, 0, 1, 2, 3]},
            'bar': {'color': "#0f172a", 'thickness': 0.34}, # Aguja oscura contrastada
            'bgcolor': "#e2e8f0",
            'steps': [
                {'range': [-3.0, -2.0], 'color': '#2563eb'},  # Azul frío
                {'range': [-2.0, -0.7], 'color': '#38bdf8'},  # Azul celeste
                {'range': [-0.7, 0.7],  'color': '#16a34a'},  # Verde confort ISO
                {'range': [0.7, 2.0],   'color': '#f59e0b'},  # Naranja
                {'range': [2.0, 3.0],   'color': '#dc2626'}   # Rojo
            ],
            'threshold': {'line': {'color': '#0f172a', 'width': 3}, 'thickness': 0.8, 'value': pmv}
        }
    ))
    fig_pmv.update_layout(height=165, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_pmv, use_container_width=True)

    # Diagnóstico y paleta
    if abs(pmv) <= 0.5:
        pmv_color = "#15803d"; pmv_bg = "#dcfce7"; pmv_border = "#86efac"; pmv_txt = "🟢 CONFORT ÓPTIMO (CAT. A/B)"
    elif abs(pmv) <= 0.7:
        pmv_color = "#0369a1"; pmv_bg = "#e0f2fe"; pmv_border = "#7dd3fc"; pmv_txt = "🔵 CONFORT ACEPTABLE (CAT. C)"
    elif pmv > 0.7:
        pmv_color = "#b91c1c"; pmv_bg = "#fee2e2"; pmv_border = "#fca5a5"; pmv_txt = "🔴 SOBRECALENTAMIENTO"
    else:
        pmv_color = "#1d4ed8"; pmv_bg = "#dbeafe"; pmv_border = "#93c5fd"; pmv_txt = "🔵 SUBENFRIAMIENTO"

    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Voto Medio Estimado / Insatisfechos</div>
        <div class="display-value" style="color: {pmv_color};">PMV: {pmv:+.2f}</div>
        <div style="font-size:0.95rem; color:#0f172a; font-weight:700; margin-top:2px;">PPD: {ppd:.1f}% insatisfechos</div>
        <div class="display-badge" style="background:{pmv_bg}; color:{pmv_color}; border:1px solid {pmv_border};">{pmv_txt}</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:5px;">Tmrt calculada (ISO 7726): <b>{t_mrt:.1f} °C</b></div>
    </div>
    """, unsafe_allow_html=True)

    # Filas activas
    r1 = "active-row" if pmv > 2.0 else ""
    r2 = "active-row" if 0.7 < pmv <= 2.0 else ""
    r3 = "active-row" if -0.5 <= pmv <= 0.5 else ""
    r4 = "active-row" if -2.0 <= pmv < -0.7 else ""
    r5 = "active-row" if pmv < -2.0 else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>PMV</th><th>Sensación</th><th>PPD</th><th>Norma</th></tr>
        <tr class="{r1}"><td>+2 a +3</td><td style="color:#dc2626; font-weight:700;">Calor severo</td><td>&gt; 75%</td><td>Fuera norma</td></tr>
        <tr class="{r2}"><td>+0.7 a +2</td><td style="color:#d97706; font-weight:700;">Lig. Cálido</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr class="{r3}"><td style="color:#16a34a; font-weight:700;">-0.5 a +0.5</td><td style="color:#16a34a; font-weight:700;">Neutro / Confort</td><td>&lt; 10%</td><td>Cat. B (Estándar)</td></tr>
        <tr class="{r4}"><td>-2 a -0.7</td><td style="color:#0284c7; font-weight:700;">Lig. Frío</td><td>15 - 75%</td><td>Cat. C</td></tr>
        <tr class="{r5}"><td>-3 a -2</td><td style="color:#2563eb; font-weight:700;">Frío severo</td><td>&gt; 75%</td><td>Fuera norma</td></tr>
    </table>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 3: WBGT (ISO 7243)
# ------------------------------------------------------------------------------
with col_b3:
    st.markdown('<div class="block-header">🏭 3. ESTRÉS TÉRMICO LABORAL | WBGT (ISO 7243)</div>', unsafe_allow_html=True)

    # Semicírculo Plotly con aguja oscuro grafito
    fig_wbgt = go.Figure(go.Indicator(
        mode="gauge",
        value=wbgt,
        gauge={
            'axis': {'range': [10.0, 42.0], 'tickwidth': 1.5, 'tickcolor': "#334155", 'tickfont': {'color': '#0f172a', 'size': 11, 'family': 'Arial Narrow'}, 'tickvals': [10, 20, limite_wbgt, 40]},
            'bar': {'color': "#0f172a", 'thickness': 0.34},
            'bgcolor': "#e2e8f0",
            'steps': [
                {'range': [10.0, limite_wbgt - 2.0], 'color': '#16a34a'},  # Zona Segura
                {'range': [limite_wbgt - 2.0, limite_wbgt], 'color': '#f59e0b'},  # Alerta
                {'range': [limite_wbgt, 42.0], 'color': '#dc2626'}   # Peligro
            ],
            'threshold': {'line': {'color': '#0f172a', 'width': 3}, 'thickness': 0.8, 'value': limite_wbgt}
        }
    ))
    fig_wbgt.update_layout(height=165, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_wbgt, use_container_width=True)

    if wbgt <= (limite_wbgt - 2.0):
        wbgt_color = "#15803d"; wbgt_bg = "#dcfce7"; wbgt_border = "#86efac"; wbgt_txt = "🟢 ZONA SEGURA CONFORME"
        diff_txt = f"{abs(delta_wbgt):.2f} °C bajo el umbral permisible"
    elif wbgt <= limite_wbgt:
        wbgt_color = "#b45309"; wbgt_bg = "#fef3c7"; wbgt_border = "#fde68a"; wbgt_txt = "🟡 ALERTA PREVENTIVA"
        diff_txt = f"A {abs(delta_wbgt):.2f} °C de superar el límite normativo"
    else:
        wbgt_color = "#b91c1c"; wbgt_bg = "#fee2e2"; wbgt_border = "#fca5a5"; wbgt_txt = "🔴 EXCEDE LÍMITE LABORAL"
        diff_txt = f"+{delta_wbgt:.2f} °C de sobrecarga térmica laboral"

    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Índice WBGT Ocupacional Medido</div>
        <div class="display-value" style="color: {wbgt_color};">{wbgt:.2f} °C</div>
        <div style="font-size:0.95rem; color:#0f172a; font-weight:700; margin-top:2px;">Límite Normativo: {limite_wbgt:.1f} °C</div>
        <div class="display-badge" style="background:{wbgt_bg}; color:{wbgt_color}; border:1px solid {wbgt_border};">{wbgt_txt}</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:5px;">{diff_txt}</div>
    </div>
    """, unsafe_allow_html=True)

    w_r1 = "active-row" if wbgt <= (limite_wbgt - 2.0) else ""
    w_r2 = "active-row" if (limite_wbgt - 2.0) < wbgt <= limite_wbgt else ""
    w_r3 = "active-row" if wbgt > limite_wbgt else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>Rango WBGT</th><th>Evaluación</th><th>Acción Preventiva</th></tr>
        <tr class="{w_r1}"><td style="color:#16a34a; font-weight:700;">&lt; {limite_wbgt - 2:.1f} °C</td><td style="color:#16a34a; font-weight:700;">Zona Conforme</td><td>Jornada continua / Normal</td></tr>
        <tr class="{w_r2}"><td style="color:#d97706; font-weight:700;">{limite_wbgt - 2:.1f} a {limite_wbgt:.1f} °C</td><td style="color:#d97706; font-weight:700;">Alerta Térmica</td><td>Pausas 15 min / Hidratación</td></tr>
        <tr class="{w_r3}"><td style="color:#dc2626; font-weight:700;">&gt; {limite_wbgt:.1f} °C</td><td style="color:#dc2626; font-weight:700;">Sobrecarga Severa</td><td>Pausas 50% / Suspender</td></tr>
    </table>
    <div style="font-size:0.72rem; color:#64748b; margin-top:6px;">*Límite dinámico según tasa metabólica seleccionada en panel.</div>
    """, unsafe_allow_html=True)

# ------------------------------------------------------------------------------
# BLOQUE 4: UTCI (COST ACTION 730 / WMO)
# ------------------------------------------------------------------------------
with col_b4:
    st.markdown('<div class="block-header">🌍 4. BIOCLIMA EXTERIOR | UTCI (COST 730 / WMO)</div>', unsafe_allow_html=True)

    # Semicírculo Plotly con aguja oscuro grafito
    fig_utci = go.Figure(go.Indicator(
        mode="gauge",
        value=utci,
        gauge={
            'axis': {'range': [-15.0, 52.0], 'tickwidth': 1.5, 'tickcolor': "#334155", 'tickfont': {'color': '#0f172a', 'size': 11, 'family': 'Arial Narrow'}, 'tickvals': [-13, 0, 9, 26, 32, 38, 46]},
            'bar': {'color': "#0f172a", 'thickness': 0.34},
            'bgcolor': "#e2e8f0",
            'steps': [
                {'range': [-15.0, 0.0], 'color': '#2563eb'},  # Azul frío
                {'range': [0.0, 9.0],   'color': '#38bdf8'},  # Frío ligero
                {'range': [9.0, 26.0],  'color': '#16a34a'},  # Confort universal
                {'range': [26.0, 32.0], 'color': '#eab308'},  # Calor moderado
                {'range': [32.0, 38.0], 'color': '#ea580c'},  # Calor fuerte
                {'range': [38.0, 46.0], 'color': '#dc2626'},  # Muy fuerte
                {'range': [46.0, 52.0], 'color': '#9333ea'}   # Extremo
            ],
            'threshold': {'line': {'color': '#0f172a', 'width': 3}, 'thickness': 0.8, 'value': utci}
        }
    ))
    fig_utci.update_layout(height=165, margin=dict(l=15, r=15, t=15, b=5), paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig_utci, use_container_width=True)

    if utci > 46:
        u_col = "#7e22ce"; u_bg = "#f3e8ff"; u_border = "#d8b4fe"; u_diag = "🟣 ESTRÉS EXTREMO CALOR"
    elif utci > 38:
        u_col = "#b91c1c"; u_bg = "#fee2e2"; u_border = "#fca5a5"; u_diag = "🔴 ESTRÉS MUY FUERTE CALOR"
    elif utci > 32:
        u_col = "#c2410c"; u_bg = "#ffedd5"; u_border = "#fed7aa"; u_diag = "🟠 ESTRÉS FUERTE CALOR"
    elif utci > 26:
        u_col = "#a16207"; u_bg = "#fef9c3"; u_border = "#fde047"; u_diag = "🟡 ESTRÉS MODERADO CALOR"
    elif utci >= 9:
        u_col = "#15803d"; u_bg = "#dcfce7"; u_border = "#86efac"; u_diag = "🟢 ZONA DE CONFORT TÉRMICO"
    elif utci >= 0:
        u_col = "#0369a1"; u_bg = "#e0f2fe"; u_border = "#7dd3fc"; u_diag = "🔵 ESTRÉS LIGERO FRÍO"
    else:
        u_col = "#1d4ed8"; u_bg = "#dbeafe"; u_border = "#93c5fd"; u_diag = "🔷 ESTRÉS FUERTE FRÍO"

    st.markdown(f"""
    <div class="display-box">
        <div class="display-label">Temperatura Climatológica Equivalente</div>
        <div class="display-value" style="color: {u_col};">{utci:+.1f} °C</div>
        <div style="font-size:0.95rem; color:#0f172a; font-weight:700; margin-top:2px;">{u_diag}</div>
        <div class="display-badge" style="background:{u_bg}; color:{u_col}; border:1px solid {u_border};">ESCALA OFICIAL WMO</div>
        <div style="font-size:0.75rem; color:#64748b; margin-top:5px;">Viento a 10m proyectado: <b>{max(0.5, v_a * 1.44):.1f} m/s</b></div>
    </div>
    """, unsafe_allow_html=True)

    u1 = "active-row" if utci > 46 else ""
    u2 = "active-row" if 38 < utci <= 46 else ""
    u3 = "active-row" if 32 < utci <= 38 else ""
    u4 = "active-row" if 26 < utci <= 32 else ""
    u5 = "active-row" if 9 <= utci <= 26 else ""
    u6 = "active-row" if utci < 9 else ""

    st.markdown(f"""
    <table class="tech-table">
        <tr><th>Rango UTCI</th><th>Sensación Fisiológica</th><th>Respuesta Humana</th></tr>
        <tr class="{u1}"><td style="color:#9333ea; font-weight:700;">&gt; +46 °C</td><td style="color:#9333ea; font-weight:700;">Calor Extremo</td><td>Peligro térmico inminente</td></tr>
        <tr class="{u2}"><td style="color:#dc2626; font-weight:700;">+38 a +46 °C</td><td style="color:#dc2626; font-weight:700;">Calor Muy Fuerte</td><td>Sudoración masiva / Taquicardia</td></tr>
        <tr class="{u3}"><td style="color:#ea580c; font-weight:700;">+32 a +38 °C</td><td style="color:#ea580c; font-weight:700;">Calor Fuerte</td><td>Sobrecarga circulatoria</td></tr>
        <tr class="{u4}"><td style="color:#ca8a04; font-weight:700;">+26 a +32 °C</td><td style="color:#ca8a04; font-weight:700;">Calor Moderado</td><td>Bochorno / Sudoración</td></tr>
        <tr class="{u5}"><td style="color:#16a34a; font-weight:700;">+9 a +26 °C</td><td style="color:#16a34a; font-weight:700;">Confort Térmico</td><td>Equilibrio óptimo</td></tr>
        <tr class="{u6}"><td style="color:#0284c7; font-weight:700;">&lt; +9 °C</td><td style="color:#0284c7; font-weight:700;">Estrés por Frío</td><td>Vasoconstricción periférica</td></tr>
    </table>
    """, unsafe_allow_html=True)
