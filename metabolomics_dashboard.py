"""
ADHD Metabolomics Insights Dashboard
Visualizes ALLCompounds_ADHD_2026_insight.xlsx with group-level comparisons
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import openpyxl
import io
import os
import re
from collections import Counter
from typing import Optional
import warnings
warnings.filterwarnings("ignore")

_HERE = os.path.dirname(os.path.abspath(__file__))

# ─── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Metabolomics Dashboard",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Custom CSS ───────────────────────────────────────────────────────────────
st.markdown("""
<style>
    div[data-testid="stMarkdownContainer"] p { color: #1a1a1a; }
    .stText, label { color: #1a1a1a; }

    .main-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        padding: 2rem; border-radius: 12px; margin-bottom: 1.5rem;
        text-align: center;
    }
    .main-header h1 { font-size: 2.4rem; margin: 0; font-weight: 700; color: #ffffff !important; }
    div[data-testid="stMarkdownContainer"] .main-header p,
    .main-header p  { margin: 0.4rem 0 0; font-size: 1.1rem; color: #e8eaf6 !important; }

    .metric-card {
        background: #f0f2f6; border-radius: 10px; padding: 1.2rem;
        border-left: 5px solid #0f3460; text-align: center;
    }
    .metric-card .value { font-size: 2rem; font-weight: 700; color: #0f3460; }
    .metric-card .label { font-size: 0.9rem; font-weight: 600; color: #333; margin-top: 0.3rem; }

    .section-header {
        background: #0f3460; color: white; padding: 0.7rem 1.2rem;
        border-radius: 8px; margin: 1.5rem 0 0.8rem;
        font-weight: 700; font-size: 1rem; letter-spacing: 0.02em;
    }
    .badge-adhd { background:#c0392b; color:white; padding:3px 10px; border-radius:12px;
                  font-size:0.85rem; font-weight:600; }
    .badge-ctrl { background:#1a6ea8; color:white; padding:3px 10px; border-radius:12px;
                  font-size:0.85rem; font-weight:600; }

    .dl-card {
        background: #f8f9fa; border: 1px solid #dee2e6;
        border-radius: 10px; padding: 1rem 1.2rem; margin-bottom: 0.5rem;
    }
    .dl-card-title { font-weight: 700; font-size: 0.95rem; color: #1a1a2e; margin-bottom: 0.2rem; }
    .dl-card-desc  { font-size: 0.82rem; color: #666; margin-bottom: 0.6rem; }
    .dl-card-count { font-size: 0.82rem; color: #0f3460; font-weight: 600; }

    .export-box {
        background: linear-gradient(135deg, #e8f4f8, #f0f7ff);
        border: 2px solid #1a6ea8; border-radius: 12px;
        padding: 1.2rem 1.5rem; margin-bottom: 1rem;
    }
    .export-box-title { font-weight: 700; font-size: 1rem; color: #0f3460; margin-bottom: 0.3rem; }
    .export-box-desc  { font-size: 0.85rem; color: #444; margin-bottom: 0.8rem; }

    div[data-testid="stExpander"] { border: 1px solid #bbb; border-radius: 8px; }
    div[data-testid="stExpander"] summary { font-size: 1rem; font-weight: 600; color: #111; }
    .css-1d391kg, .css-qrbaxs { color: #111 !important; font-weight: 500; }
    button[data-baseweb="tab"] { font-size: 0.95rem; font-weight: 600; }
</style>
""", unsafe_allow_html=True)

# ─── Plot template ─────────────────────────────────────────────────────────────
PLOT_TEMPLATE = dict(
    layout=dict(
        font=dict(family="Arial, sans-serif", size=13, color="#111"),
        title=dict(font=dict(size=15, color="#111"), x=0.02),
        xaxis=dict(tickfont=dict(size=12, color="#222"), title_font=dict(size=13, color="#222"),
                   linecolor="#888", gridcolor="#e0e0e0"),
        yaxis=dict(tickfont=dict(size=12, color="#222"), title_font=dict(size=13, color="#222"),
                   linecolor="#888", gridcolor="#e0e0e0"),
        paper_bgcolor="white", plot_bgcolor="#fafafa",
        legend=dict(font=dict(size=12, color="#111")),
    )
)

# ─── Constants ────────────────────────────────────────────────────────────────
CLINICAL_FILE   = os.path.join(_HERE, "data", "ADHD_all_data.xlsx")
DEFAULT_INSIGHT = os.path.join(_HERE, "data", "ALLCompounds_ADHD_2026_insight.xlsx")

COLORS = {"ADHD": "#c0392b", "Control": "#1a6ea8", "neutral": "#7f8c8d"}

CAT_COLORS = {
    "BIOLOGICAL":   "#27ae60",
    "CLINICAL":     "#2980b9",
    "FUNCTIONAL":   "#e67e22",
    "SENSORY":      "#8e44ad",
    "UNCLASSIFIED": "#7f8c8d",
}

# Standard export columns (no intensity columns)
EXPORT_BASE = ['Original_Name', 'Formula_Final', 'MW_Final', 'HMDB_ID', 'KEGG_ID_Final',
               'Main_Categories', 'LM_Category', 'LM_Main_Class', 'LM_Sub_Class',
               'OC_Ratio', 'NOSC', 'KEGG_Pathways', 'Neuro_Trap',
               'Clinical_Tags', 'Structural_Tags']

# ─── Download Helpers ──────────────────────────────────────────────────────────

def to_csv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode('utf-8')


def to_excel(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine='openpyxl') as w:
        df.to_excel(w, index=False)
    return buf.getvalue()


def export_cols(df: pd.DataFrame, name_col: str, fc_col: str, pval_col: str, direction: bool = True) -> list:
    """Return available export columns from EXPORT_BASE + stats, present in df."""
    stat_cols = [fc_col, pval_col, '_direction'] if direction else []
    cols = [name_col] + stat_cols + [c for c in EXPORT_BASE if c in df.columns and c != name_col]
    return [c for c in cols if c in df.columns]


def build_metaboanalyst(df: pd.DataFrame, sample_cols: list, grp_norm: dict):
    """Build MetaboAnalyst peak table and metadata DataFrames."""
    def best_id(row):
        pubchem = str(row.get('PubChem_CID', '') or '').strip()
        hmdb    = str(row.get('HMDB_ID',     '') or '').strip()
        kegg    = str(row.get('KEGG_ID_Final','') or '').strip()
        if pubchem and pubchem not in ('nan', '-', ''):
            return pubchem
        if hmdb.startswith('HMDB') and len(hmdb) > 4:
            return hmdb
        if kegg and kegg not in ('-', 'nan', ''):
            return kegg
        return None

    df2 = df.copy()
    df2['_ma_id'] = df2.apply(best_id, axis=1)
    included = df2[df2['_ma_id'].notna()].copy()
    n_excluded = len(df2) - len(included)

    # Peak table: rows = compounds, cols = samples
    peak_rows = []
    for _, row in included.iterrows():
        r = {'Label': row['_ma_id']}
        for s in sample_cols:
            val = pd.to_numeric(row.get(s, np.nan), errors='coerce')
            r[s] = 0 if pd.isna(val) else val
        peak_rows.append(r)
    peak_table = pd.DataFrame(peak_rows)

    # Metadata
    meta = pd.DataFrame([
        {'Sample': s, 'Group': grp_norm.get(s, 'Unknown')}
        for s in sample_cols
    ])

    return peak_table, meta, len(included), n_excluded


def dl_row(label: str, desc: str, df: pd.DataFrame, key: str,
           count: Optional[int] = None, name_col_rename: Optional[str] = None):
    """Render a compact download card with CSV + Excel buttons."""
    count_str = f"&nbsp;·&nbsp;<span class='dl-card-count'>{count:,} compounds</span>" if count is not None else ""
    st.markdown(
        f"<div class='dl-card'>"
        f"<div class='dl-card-title'>{label}</div>"
        f"<div class='dl-card-desc'>{desc}{count_str}</div>"
        f"</div>",
        unsafe_allow_html=True
    )
    export_df = df.rename(columns={name_col_rename: 'Compound'}) if name_col_rename else df
    c1, c2, _ = st.columns([1, 1, 4])
    c1.download_button("⬇ CSV",   to_csv(export_df),   file_name=f"{key}.csv",
                       mime='text/csv',             key=f"csv_{key}", use_container_width=True)
    c2.download_button("⬇ Excel", to_excel(export_df), file_name=f"{key}.xlsx",
                       mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                       key=f"xl_{key}", use_container_width=True)


# ─── Data loaders ─────────────────────────────────────────────────────────────

def code_to_col(clinical_code: str) -> str:
    m = re.match(r'^(\d+)([A-Za-z]+)$', str(clinical_code))
    return (m.group(2).upper() + m.group(1)) if m else str(clinical_code)


@st.cache_data(show_spinner="Loading clinical group data…")
def load_group_map():
    wb = openpyxl.load_workbook(CLINICAL_FILE, read_only=True)
    ws = wb['ALL_data']
    mapping = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        code, group = row[5], row[6]
        if code and group:
            mapping[code_to_col(code)] = group.strip().upper()
    return mapping


@st.cache_data(show_spinner="Parsing insights file…")
def load_insights(path: str):
    wb = openpyxl.load_workbook(path, read_only=True)
    ws = wb['Sheet1']
    headers = list(next(ws.iter_rows(min_row=1, max_row=1, values_only=True)))
    rows    = list(ws.iter_rows(min_row=2, values_only=True))
    return pd.DataFrame(rows, columns=headers)


def detect_sample_cols(df):
    return [c for c in df.columns if isinstance(c, str) and re.match(r'^[AB]\d+$', c)]


def explode_tags(series: pd.Series) -> pd.Series:
    return (series.dropna()
                  .apply(lambda x: [t.strip() for t in str(x).split(';')
                                    if t.strip() and t.strip() != '-'])
                  .explode()
                  .reset_index(drop=True))


def apply_template(fig):
    fig.update_layout(**PLOT_TEMPLATE['layout'])
    return fig


# ─── Main app ─────────────────────────────────────────────────────────────────

def main():
    st.markdown("""
    <div class="main-header">
        <h1>⚗️ Metabolomics Dashboard</h1>
        <p>Untargeted Metabolomics Analysis · Tel-Hai Nutrition &amp; Bioinformatics Lab</p>
    </div>""", unsafe_allow_html=True)

    # ── Sidebar ──────────────────────────────────────────────────────────────
    with st.sidebar:
        st.header("📂 Data Source")
        uploaded = st.file_uploader("Upload Insights XLSX", type=["xlsx"])
        if uploaded:
            import tempfile
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx")
            tmp.write(uploaded.read()); tmp.close()
            insights_path = tmp.name
        else:
            insights_path = DEFAULT_INSIGHT
            st.info("Default: ALLCompounds_ADHD_2026_insight.xlsx")

        st.divider()
        st.header("🔬 Statistical Filters")
        fc_thresh   = st.slider("|Log₂ FC| threshold", 0.0, 4.0, 1.0, 0.25)
        pval_thresh = st.selectbox("P-value threshold",
                                   [1.0, 0.5, 0.1, 0.05, 0.01, 0.001], index=3,
                                   format_func=lambda v: "Any (no filter)" if v == 1.0 else str(v))
        show_adj    = st.checkbox("Use adjusted p-value (FDR)", value=False)
        st.divider()
        st.caption("Tel-Hai College · Nutrition Lab · 2026")

    # ── Load data ────────────────────────────────────────────────────────────
    group_map   = load_group_map()
    df          = load_insights(insights_path)
    sample_cols = detect_sample_cols(df)

    grp_norm     = {k: ("ADHD" if v == "ADHD" else "Control") for k, v in group_map.items()}
    adhd_cols    = [c for c in sample_cols if grp_norm.get(c) == "ADHD"]
    control_cols = [c for c in sample_cols if grp_norm.get(c) == "Control"]

    fc_col   = 'Log2 Fold Change: (Control) / (ADHD)'
    pval_col = 'P-value: (Control) / (ADHD)'
    adjp_col = 'Adj. P-value: (Control) / (ADHD)'
    p_use    = adjp_col if show_adj else pval_col
    name_col = 'Name' if 'Name' in df.columns else 'Original_Name'

    for col in [fc_col, pval_col, adjp_col, 'OC_Ratio', 'NOSC']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    df['_direction'] = df.apply(
        lambda r: ('↑ ADHD' if r.get(fc_col, np.nan) < -fc_thresh
                   else ('↑ Control' if r.get(fc_col, np.nan) > fc_thresh else 'n/s'))
        if (pd.notna(r.get(fc_col)) and pd.notna(r.get(p_use))
            and abs(r.get(fc_col, 0)) >= fc_thresh
            and r.get(p_use, 1) < pval_thresh)
        else 'n/s', axis=1)

    # ── Project banner — group names extracted from column name ───────────────
    _m = re.search(r'\(([^)]+)\)\s*/\s*\(([^)]+)\)', fc_col)
    group_b_name = _m.group(1) if _m else "Group B"   # Control
    group_a_name = _m.group(2) if _m else "Group A"   # ADHD
    st.markdown(
        f"<div style='background:#f0f7ff;border-left:4px solid #0f3460;"
        f"padding:0.6rem 1.2rem;border-radius:6px;margin-bottom:1rem;"
        f"font-size:0.95rem;color:#1a1a2e'>"
        f"<b>Study:</b> {group_a_name} vs {group_b_name} &nbsp;·&nbsp; "
        f"<b>{len(df):,}</b> compounds &nbsp;·&nbsp; "
        f"<b>{len(sample_cols)}</b> samples "
        f"({len(adhd_cols)} {group_a_name} · {len(control_cols)} {group_b_name})"
        f"</div>",
        unsafe_allow_html=True
    )

    # ── TABS ─────────────────────────────────────────────────────────────────
    tabs = st.tabs([
        "📊 Overview",
        "🧪 Lipids",
        "⚡ Oxidation",
        "🗺 KEGG Pathways",
        "🧠 Neurochemistry",
        f"📈 {group_a_name} vs {group_b_name}",
        "🔍 Compound Explorer",
        "📥 Downloads & Exports",
    ])

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 1 — OVERVIEW
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[0]:
        n_sig_adhd = int((df['_direction'] == '↑ ADHD').sum())
        n_sig_ctrl = int((df['_direction'] == '↑ Control').sum())

        metrics = [
            (len(df),                                          "Total Compounds"),
            (len(adhd_cols),                                   "ADHD Samples"),
            (len(control_cols),                                "Control Samples"),
            (int(df['LM_Category'].notna().sum()),             "Classified Lipids"),
            (int((df['KEGG_Pathways'].notna() & (df['KEGG_Pathways'] != '-')).sum()), "In KEGG Pathways"),
            (int((df['Neuro_Trap'].notna() & (df['Neuro_Trap'] != '-')).sum()),        "Neuro-Active"),
        ]
        cols = st.columns(6)
        for col, (val, label) in zip(cols, metrics):
            col.markdown(f"""<div class="metric-card">
                <div class="value">{val:,}</div>
                <div class="label">{label}</div></div>""", unsafe_allow_html=True)

        st.markdown("<div class='section-header'>Functional Classification of All Compounds</div>",
                    unsafe_allow_html=True)
        cl, cr = st.columns(2)

        cat_ctr = Counter()
        for v in df['Main_Categories'].dropna():
            for c in str(v).split(';'):
                c = c.strip()
                if c: cat_ctr[c] += 1
        cat_df = pd.DataFrame(cat_ctr.items(), columns=['Category','Count'])
        fig_pie = px.pie(cat_df, names='Category', values='Count',
                         color='Category', color_discrete_map=CAT_COLORS,
                         title='Compound Categories (multi-label)', hole=0.42)
        fig_pie.update_traces(textinfo='percent+label', textfont_size=13, pull=[0.04]*len(cat_df))
        fig_pie.update_layout(showlegend=False, margin=dict(t=50,b=10))
        apply_template(fig_pie)
        cl.plotly_chart(fig_pie, use_container_width=True)

        stag = explode_tags(df['Structural_Tags'])
        stag_df = stag.value_counts().head(15).reset_index()
        stag_df.columns = ['Tag','Count']
        stag_df['Tag'] = stag_df['Tag'].str.split(':').str[-1].str.replace('_',' ')
        fig_stag = px.bar(stag_df, x='Count', y='Tag', orientation='h',
                          color='Count', color_continuous_scale='Blues',
                          title='Top 15 Structural / Biological Tags')
        fig_stag.update_layout(yaxis={'categoryorder':'total ascending'},
                                coloraxis_showscale=False, margin=dict(t=50,b=10))
        apply_template(fig_stag)
        cr.plotly_chart(fig_stag, use_container_width=True)

        st.markdown("<div class='section-header'>Significance Summary (current filter thresholds)</div>",
                    unsafe_allow_html=True)
        sc1, sc2, sc3 = st.columns(3)
        sc1.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#c0392b">{n_sig_adhd}</div>
            <div class="label">↑ Elevated in ADHD</div></div>""", unsafe_allow_html=True)
        sc2.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#1a6ea8">{n_sig_ctrl}</div>
            <div class="label">↑ Elevated in Control</div></div>""", unsafe_allow_html=True)
        sc3.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#555">{len(df)-n_sig_adhd-n_sig_ctrl:,}</div>
            <div class="label">Not Significant</div></div>""", unsafe_allow_html=True)

        st.markdown("<div class='section-header'>Sample Group Map</div>", unsafe_allow_html=True)
        smap = pd.DataFrame([{'Sample': s, 'Group': grp_norm.get(s, 'Unknown')} for s in sample_cols])
        fig_map = px.strip(smap, x='Sample', color='Group', color_discrete_map=COLORS,
                           category_orders={'Group': ['ADHD','Control']},
                           title=f'{len(sample_cols)} samples — {len(adhd_cols)} ADHD · {len(control_cols)} Control')
        fig_map.update_traces(marker_size=14, jitter=0)
        fig_map.update_layout(xaxis_tickangle=45, margin=dict(t=50,b=90), xaxis_tickfont_size=11)
        apply_template(fig_map)
        st.plotly_chart(fig_map, use_container_width=True)

        with st.expander("📋 Full sample list by group"):
            ea, ec = st.columns(2)
            ea.markdown("**ADHD samples**");  ea.write(sorted(adhd_cols))
            ec.markdown("**Control samples**"); ec.write(sorted(control_cols))

        # ── Quick Download ──────────────────────────────────────────────────
        st.markdown("<div class='section-header'>Quick Download — Full Dataset</div>",
                    unsafe_allow_html=True)
        ecols = export_cols(df, name_col, fc_col, pval_col)
        dl_row("Full Annotated Dataset",
               "All compounds with annotation, statistics, and classification",
               df[ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                          fc_col:'Log2FC', pval_col:'p-value'}),
               key="overview_full", count=len(df), name_col_rename=None)

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 2 — LIPIDS
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[1]:
        lipid_df = df[df['LM_Category'].notna()].copy()
        st.markdown(f"### {len(lipid_df):,} classified lipids (LIPID MAPS annotation)")

        cl, cr = st.columns(2)
        sb_df = lipid_df[['LM_Category','LM_Main_Class','LM_Sub_Class']].fillna('Unknown')
        sb_df = sb_df[sb_df['LM_Category'] != 'Unknown'].copy()
        sb_df['n'] = 1
        fig_sun = px.sunburst(sb_df, path=['LM_Category','LM_Main_Class','LM_Sub_Class'],
                              values='n', title='Lipid Class Hierarchy — Sunburst',
                              color_discrete_sequence=px.colors.qualitative.Set2)
        fig_sun.update_layout(margin=dict(t=50,b=10), title_font_size=15)
        cl.plotly_chart(fig_sun, use_container_width=True)

        lm_mc = lipid_df['LM_Main_Class'].value_counts().head(18).reset_index()
        lm_mc.columns = ['Class','Count']
        lm_mc['Class'] = lm_mc['Class'].str.replace(r'\s*\[.*?\]', '', regex=True)
        fig_lm = px.bar(lm_mc, x='Count', y='Class', orientation='h',
                        color='Count', color_continuous_scale='Teal',
                        title='Lipid Main Classes (top 18)')
        fig_lm.update_layout(yaxis={'categoryorder':'total ascending'},
                              coloraxis_showscale=False, margin=dict(t=50,b=10))
        apply_template(fig_lm)
        cr.plotly_chart(fig_lm, use_container_width=True)

        # Download all lipids
        lipid_ecols = export_cols(lipid_df, name_col, fc_col, pval_col)
        dl_row("All Classified Lipids",
               "All LIPID MAPS annotated compounds with statistics",
               lipid_df[lipid_ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                       fc_col:'Log2FC', pval_col:'p-value'}),
               key="lipids_all", count=len(lipid_df))

        st.markdown("<div class='section-header'>Lipid Sub-Classes — Expand for Group Comparison</div>",
                    unsafe_allow_html=True)
        lm_sub_top = lipid_df['LM_Sub_Class'].value_counts()
        for _lm_idx, (sub_class, count) in enumerate(lm_sub_top.items()):
            if sub_class in ('-','Unknown') or not sub_class: continue
            sub_df = lipid_df[lipid_df['LM_Sub_Class'] == sub_class]
            short  = re.sub(r'\s*\[.*?\]','', sub_class)
            n_adhd_up = int((sub_df['_direction'] == '↑ ADHD').sum())
            n_ctrl_up = int((sub_df['_direction'] == '↑ Control').sum())
            n_ns      = count - n_adhd_up - n_ctrl_up
            med_fc    = sub_df[fc_col].median() if fc_col in sub_df.columns else np.nan
            with st.expander(f"**{short}** — {count} compounds"):
                g1, g2 = st.columns(2)
                with g1:
                    st.markdown("**Directional compound counts:**")
                    st.markdown(f"<span class='badge-adhd'>↑ ADHD: {n_adhd_up}</span>&nbsp;&nbsp;"
                                f"<span class='badge-ctrl'>↑ Control: {n_ctrl_up}</span>&nbsp;&nbsp;"
                                f"<span style='background:#888;color:white;padding:3px 8px;border-radius:12px;"
                                f"font-size:0.85rem;font-weight:600'>n/s: {n_ns}</span>",
                                unsafe_allow_html=True)
                    if np.isfinite(med_fc):
                        direction = "↑ Control" if med_fc > 0 else "↑ ADHD"
                        st.metric("Median Log₂FC (Control/ADHD)", f"{med_fc:+.2f}",
                                  help=f"Positive = higher in Control. {direction}")
                    if n_adhd_up > 0 or n_ctrl_up > 0:
                        mini = pd.DataFrame({
                            'Group': ['↑ ADHD','↑ Control'],
                            'Count': [-n_adhd_up, n_ctrl_up],
                        })
                        fig_mini = px.bar(mini, x='Count', y='Group', orientation='h',
                                          color='Group',
                                          color_discrete_map={'↑ ADHD':COLORS['ADHD'],
                                                              '↑ Control':COLORS['Control']},
                                          title='Significant compounds')
                        fig_mini.update_layout(showlegend=False, height=120,
                                               margin=dict(t=30,b=5,l=5,r=5),
                                               xaxis=dict(tickformat='d', zeroline=True,
                                                          zerolinecolor='#333', zerolinewidth=2),
                                               font=dict(size=12,color='#111'),
                                               paper_bgcolor='white', plot_bgcolor='#fafafa')
                        st.plotly_chart(fig_mini, use_container_width=True,
                                        key=f"lipid_mini_{_lm_idx}")
                    # Download this sub-class
                    sub_ecols = export_cols(sub_df, name_col, fc_col, pval_col)
                    sub_export = sub_df[sub_ecols].rename(
                        columns={name_col:'Compound', '_direction':'Direction',
                                 fc_col:'Log2FC', pval_col:'p-value'})
                    sub_key = re.sub(r'[^\w]', '_', short)[:30]
                    c1, c2, _ = st.columns([1,1,4])
                    c1.download_button("⬇ CSV", to_csv(sub_export),
                                       file_name=f"lipid_{sub_key}.csv", mime='text/csv',
                                       key=f"csv_lip_{_lm_idx}", use_container_width=True)
                    c2.download_button("⬇ Excel", to_excel(sub_export),
                                       file_name=f"lipid_{sub_key}.xlsx",
                                       mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                       key=f"xl_lip_{_lm_idx}", use_container_width=True)
                with g2:
                    names = sub_df[name_col].dropna().unique()
                    st.markdown(f"**{len(names)} compounds:**")
                    st.write(list(names[:30]))

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 3 — OXIDATION
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[2]:
        st.markdown("### Oxidative Stress Markers")
        oc_df = df[df['OC_Ratio'].notna()].copy()
        oc_df['Oxidation_Level'] = pd.cut(oc_df['OC_Ratio'],
            bins=[-np.inf, 0.2, 0.4, 0.6, 0.8, np.inf],
            labels=['Very Low','Low','Medium','High','Very High'])

        cl, cr = st.columns(2)
        fig_oc = px.histogram(oc_df, x='OC_Ratio', nbins=60,
                              color_discrete_sequence=['#c0392b'],
                              title='O/C Ratio Distribution — Oxidation Index',
                              labels={'OC_Ratio':'O/C Ratio', 'count':'# Compounds'})
        fig_oc.add_vline(x=0.6, line_dash='dash', line_color='#1a1a8c',
                         annotation_text='  High oxidation threshold (0.6)',
                         annotation_font_size=12)
        fig_oc.update_layout(margin=dict(t=50,b=10))
        apply_template(fig_oc)
        cl.plotly_chart(fig_oc, use_container_width=True)

        ox_c = oc_df['Oxidation_Level'].value_counts().reset_index()
        ox_c.columns = ['Level','Count']
        fig_ox_pie = px.pie(ox_c, names='Level', values='Count', color='Level',
                            color_discrete_sequence=['#2ecc71','#f1c40f','#e67e22','#e74c3c','#8e44ad'],
                            title='Compounds by Oxidation Level', hole=0.4)
        fig_ox_pie.update_traces(textinfo='percent+label', textfont_size=12,
                                  textposition='auto',
                                  insidetextorientation='horizontal')
        fig_ox_pie.update_layout(
            margin=dict(t=50, b=80, l=40, r=40),
            uniformtext_minsize=9,
            uniformtext_mode='hide',
            showlegend=True,
            legend=dict(orientation='v', x=1.02, y=0.5,
                        font=dict(size=12), bgcolor='rgba(0,0,0,0)'),
        )
        cr.plotly_chart(fig_ox_pie, use_container_width=True)

        st.markdown("<div class='section-header'>Clinical Oxidative Stress Categories</div>",
                    unsafe_allow_html=True)
        ox_kw = ['Oxylipin','Epoxide','Oxidative','Reactive','Peroxide']
        ox_dict = Counter()
        for v in df['Clinical_Tags'].dropna():
            for tag in str(v).split(';'):
                tag = tag.strip()
                if any(k in tag for k in ox_kw):
                    ox_dict[tag.split(':')[-1].replace('_',' ')] += 1
        if ox_dict:
            ox_df2 = pd.DataFrame(ox_dict.most_common(12), columns=['Type','Count'])
            fig_oxt = px.bar(ox_df2, x='Type', y='Count',
                             color='Count', color_continuous_scale='Reds',
                             title='Oxidative Stress Compound Types')
            fig_oxt.update_layout(coloraxis_showscale=False,
                                   margin=dict(t=50,b=90), xaxis_tickangle=30)
            apply_template(fig_oxt)
            st.plotly_chart(fig_oxt, use_container_width=True)

        st.markdown("<div class='section-header'>High-Oxidation Compounds (O/C > 0.6) — ADHD vs Control</div>",
                    unsafe_allow_html=True)
        high_ox = oc_df[oc_df['OC_Ratio'] > 0.6].copy()
        n_hox_adhd = int((high_ox['_direction'] == '↑ ADHD').sum())
        n_hox_ctrl = int((high_ox['_direction'] == '↑ Control').sum())
        n_hox_ns   = len(high_ox) - n_hox_adhd - n_hox_ctrl

        hox_c1, hox_c2, hox_c3 = st.columns(3)
        hox_c1.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#c0392b">{n_hox_adhd}</div>
            <div class="label">↑ ADHD (sig.)</div></div>""", unsafe_allow_html=True)
        hox_c2.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#1a6ea8">{n_hox_ctrl}</div>
            <div class="label">↑ Control (sig.)</div></div>""", unsafe_allow_html=True)
        hox_c3.markdown(f"""<div class="metric-card">
            <div class="value" style="color:#555">{n_hox_ns}</div>
            <div class="label">Not significant</div></div>""", unsafe_allow_html=True)

        hox_fc = high_ox[high_ox[fc_col].notna() & high_ox[pval_col].notna()].copy()
        hox_fc = hox_fc.drop_duplicates(subset=[name_col]).nlargest(30, 'OC_Ratio')
        hox_fc['label_color'] = hox_fc['_direction'].map(
            {'↑ ADHD': COLORS['ADHD'], '↑ Control': COLORS['Control'], 'n/s': '#aaa'})
        hox_fc_sorted = hox_fc.sort_values(fc_col)
        fig_hox = go.Figure(go.Bar(
            x=hox_fc_sorted[fc_col], y=hox_fc_sorted[name_col], orientation='h',
            marker_color=hox_fc_sorted['label_color'].tolist(),
            hovertemplate='<b>%{y}</b><br>Log₂FC: %{x:.2f}<extra></extra>',
        ))
        fig_hox.add_vline(x=0, line_color='#333', line_width=1.5)
        fig_hox.add_vline(x= fc_thresh, line_dash='dash', line_color='#888', line_width=1)
        fig_hox.add_vline(x=-fc_thresh, line_dash='dash', line_color='#888', line_width=1)
        fig_hox.update_layout(
            title=dict(text='Top 30 High-Oxidation Compounds (O/C > 0.6) — Log₂ Fold Change',
                       font_size=14),
            xaxis=dict(title='Log₂ FC (Control / ADHD)', zeroline=True, zerolinecolor='#333'),
            height=max(400, len(hox_fc_sorted) * 20),
            margin=dict(t=80, b=30, l=10, r=20),
            paper_bgcolor='white', plot_bgcolor='#fafafa', font=dict(size=12, color='#111'),
        )
        st.plotly_chart(fig_hox, use_container_width=True)

        hox_ecols = export_cols(high_ox, name_col, fc_col, pval_col)
        with st.expander(f"📋 All {len(high_ox)} highly-oxidised compounds (O/C > 0.6)"):
            show = [name_col,'OC_Ratio','NOSC', fc_col, pval_col, '_direction','Main_Categories','Clinical_Tags']
            show = [c for c in show if c in high_ox.columns]
            st.dataframe(high_ox[show].rename(columns={fc_col:'Log2FC', pval_col:'p-value',
                                                        '_direction':'Direction'})
                         .reset_index(drop=True), use_container_width=True)

        dl_row("High-Oxidation Compounds (O/C > 0.6)",
               "All compounds above the high oxidation threshold",
               high_ox[hox_ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                    fc_col:'Log2FC', pval_col:'p-value'}),
               key="high_ox", count=len(high_ox))

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 4 — KEGG PATHWAYS
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[3]:
        kegg_df = df[df['KEGG_Pathways'].notna() & (df['KEGG_Pathways'] != '-')].copy()
        st.markdown(f"### {len(kegg_df):,} compounds mapped to KEGG pathways")

        path_all  = Counter()
        path_adhd = Counter()
        path_ctrl = Counter()
        for _, row in kegg_df.iterrows():
            paths = [p.strip() for p in str(row['KEGG_Pathways']).split(';') if p.strip()]
            fc = row.get(fc_col, np.nan); pv = row.get(p_use, np.nan)
            is_adhd_up = pd.notna(fc) and pd.notna(pv) and fc < -fc_thresh and pv < pval_thresh
            is_ctrl_up = pd.notna(fc) and pd.notna(pv) and fc >  fc_thresh and pv < pval_thresh
            for p in paths:
                path_all[p]  += 1
                if is_adhd_up: path_adhd[p] += 1
                if is_ctrl_up: path_ctrl[p]  += 1

        all_pathway_names = [p for p, _ in path_all.most_common()]
        selected_paths = st.multiselect(
            "Pathways to display (add or remove freely):",
            options=all_pathway_names, default=all_pathway_names[:25], key='kegg_multiselect')

        if not selected_paths:
            st.info("Select at least one pathway to display the chart.")
        else:
            pw_df = pd.DataFrame({
                'Pathway':   selected_paths,
                'ADHD_n':    [path_adhd[p] for p in selected_paths],
                'Control_n': [path_ctrl[p] for p in selected_paths],
                'Total':     [path_all[p]  for p in selected_paths],
            }).sort_values('ADHD_n')

            bar_h = max(28, min(48, 600 // max(len(pw_df), 1)))
            fig_h = max(420, len(pw_df) * bar_h + 160)
            pv_label = "Any" if pval_thresh >= 1.0 else f"p<{pval_thresh}"
            fc_label = f"|FC|≥{fc_thresh}" if fc_thresh > 0 else "any FC"

            fig_mirror = go.Figure()
            fig_mirror.add_trace(go.Bar(
                name='↑ ADHD', x=[-n for n in pw_df['ADHD_n']], y=pw_df['Pathway'],
                orientation='h', marker=dict(color=COLORS['ADHD'], line=dict(color='white', width=0.8)),
                customdata=pw_df['ADHD_n'].values,
                hovertemplate='<b>%{y}</b><br>↑ ADHD: <b>%{customdata}</b><extra></extra>',
                text=[str(n) if n > 0 else '' for n in pw_df['ADHD_n']],
                textposition='inside', textfont=dict(size=13, color='white'), insidetextanchor='middle',
            ))
            fig_mirror.add_trace(go.Bar(
                name='↑ Control', x=pw_df['Control_n'], y=pw_df['Pathway'],
                orientation='h', marker=dict(color=COLORS['Control'], line=dict(color='white', width=0.8)),
                customdata=pw_df['Control_n'].values,
                hovertemplate='<b>%{y}</b><br>↑ Control: <b>%{customdata}</b><extra></extra>',
                text=[str(n) if n > 0 else '' for n in pw_df['Control_n']],
                textposition='inside', textfont=dict(size=13, color='white'), insidetextanchor='middle',
            ))
            for i, row_pw in pw_df.iterrows():
                if row_pw['Total'] > 0:
                    fig_mirror.add_annotation(
                        x=max(pw_df['Control_n'].max(), 1) * 1.15, y=row_pw['Pathway'],
                        text=f"n={row_pw['Total']}", showarrow=False,
                        font=dict(size=10, color='#666'), xanchor='left')

            x_max = max(pw_df['ADHD_n'].max(), pw_df['Control_n'].max(), 1)
            x_pad = x_max * 0.25
            fig_mirror.update_layout(
                barmode='relative',
                title=dict(text=(f'<b>Metabolites per KEGG Pathway</b>'
                                 f'<br><span style="font-size:12px;color:#555">'
                                 f'Left ← ↑ ADHD | ↑ Control → Right | {fc_label}, {pv_label}</span>'),
                           font=dict(size=16, color='#111'), x=0.5, xanchor='center'),
                xaxis=dict(title=dict(text='Number of metabolites', font=dict(size=13)),
                           zeroline=True, zerolinecolor='#222', zerolinewidth=2,
                           gridcolor='#e8e8e8', range=[-(x_max+x_pad), x_max+x_pad*3],
                           tickmode='array',
                           tickvals=list(range(-(x_max+1), x_max+2)),
                           ticktext=[str(abs(v)) for v in range(-(x_max+1), x_max+2)]),
                yaxis=dict(tickfont=dict(size=12, color='#111'), automargin=True),
                legend=dict(orientation='h', x=0.5, y=-0.08, xanchor='center',
                            font=dict(size=13), bgcolor='rgba(0,0,0,0)'),
                height=fig_h, margin=dict(t=90, b=70, l=20, r=80),
                paper_bgcolor='white', plot_bgcolor='white',
                font=dict(family='Arial, sans-serif', size=13, color='#111'), bargap=0.35,
            )
            for i, pathway in enumerate(pw_df['Pathway']):
                if i % 2 == 0:
                    fig_mirror.add_hrect(y0=i-0.5, y1=i+0.5,
                                          fillcolor='rgba(0,0,0,0.025)', line_width=0)
            st.plotly_chart(fig_mirror, use_container_width=True, key='kegg_mirror')

            st.markdown("<div class='section-header'>Pathway Detail — Compounds &amp; Download</div>",
                        unsafe_allow_html=True)
            for _pw_idx, path in enumerate(sorted(selected_paths, key=lambda p: -path_all[p])):
                path_cpds = kegg_df[kegg_df['KEGG_Pathways'].str.contains(re.escape(path), na=False)]
                adhd_up   = path_cpds[path_cpds['_direction'] == '↑ ADHD']
                ctrl_up   = path_cpds[path_cpds['_direction'] == '↑ Control']
                label = (f"**{path}** — {path_all[path]} total "
                         f"| ↑ADHD: {path_adhd[path]} | ↑Control: {path_ctrl[path]}")
                with st.expander(label):
                    ec1, ec2 = st.columns(2)
                    with ec1:
                        st.markdown(f"<span class='badge-adhd'>↑ ADHD ({len(adhd_up)})</span>",
                                    unsafe_allow_html=True)
                        st.write(list(adhd_up[name_col].dropna().unique()))
                    with ec2:
                        st.markdown(f"<span class='badge-ctrl'>↑ Control ({len(ctrl_up)})</span>",
                                    unsafe_allow_html=True)
                        st.write(list(ctrl_up[name_col].dropna().unique()))
                    st.markdown("**All compounds in this pathway:**")
                    st.write(list(path_cpds[name_col].dropna().unique()))

                    # Download this pathway
                    pw_ecols = export_cols(path_cpds, name_col, fc_col, pval_col)
                    pw_export = path_cpds[pw_ecols].rename(
                        columns={name_col:'Compound', '_direction':'Direction',
                                 fc_col:'Log2FC', pval_col:'p-value'})
                    pw_key = re.sub(r'[^\w]', '_', path)[:30]
                    c1, c2, _ = st.columns([1,1,4])
                    c1.download_button("⬇ CSV", to_csv(pw_export),
                                       file_name=f"kegg_{pw_key}.csv", mime='text/csv',
                                       key=f"csv_kegg_{_pw_idx}", use_container_width=True)
                    c2.download_button("⬇ Excel", to_excel(pw_export),
                                       file_name=f"kegg_{pw_key}.xlsx",
                                       mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                       key=f"xl_kegg_{_pw_idx}", use_container_width=True)

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 5 — NEUROCHEMISTRY
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[4]:
        neuro_df = df[df['Neuro_Trap'].notna() & (df['Neuro_Trap'] != '-')].copy()
        st.markdown(f"### {len(neuro_df):,} neuro-active compounds detected")

        neuro_ctr = Counter()
        for v in neuro_df['Neuro_Trap']:
            for t in str(v).split(';'):
                t = t.strip()
                if t: neuro_ctr[t] += 1

        neuro_top = pd.DataFrame(neuro_ctr.most_common(20), columns=['Axis','Count'])
        cl, cr = st.columns([3,2])
        fig_neuro = px.bar(neuro_top, x='Count', y='Axis', orientation='h',
                           color='Count', color_continuous_scale='Purples',
                           title='Neuro-active Compound Axes')
        fig_neuro.update_layout(yaxis={'categoryorder':'total ascending'},
                                 coloraxis_showscale=False, margin=dict(t=50,b=10))
        apply_template(fig_neuro)
        cl.plotly_chart(fig_neuro, use_container_width=True)

        major_axes = [a for a in neuro_ctr if 'axis' in a.lower() or '[Pathway' in a][:8]
        if major_axes:
            rv = [neuro_ctr[a] for a in major_axes]
            fig_radar = go.Figure(go.Scatterpolar(
                r=rv + [rv[0]], theta=major_axes + [major_axes[0]],
                fill='toself', line_color='#6c3483', fillcolor='rgba(108,52,131,0.25)',
            ))
            fig_radar.update_layout(
                polar=dict(radialaxis=dict(visible=True, tickfont_size=11)),
                title=dict(text='NT Axis Coverage', font_size=14),
                margin=dict(t=60,b=20), font=dict(size=12, color='#111'),
            )
            cr.plotly_chart(fig_radar, use_container_width=True)

        # Download all neuro-active
        neuro_ecols = export_cols(neuro_df, name_col, fc_col, pval_col)
        dl_row("All Neuro-Active Compounds",
               "Compounds with Neuro_Trap annotation (all axes)",
               neuro_df[neuro_ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                       fc_col:'Log2FC', pval_col:'p-value'}),
               key="neuro_all", count=len(neuro_df))

        st.markdown("<div class='section-header'>Neuro-Axis Detail — Compounds &amp; Download</div>",
                    unsafe_allow_html=True)
        for _na_idx, (axis, cnt) in enumerate(neuro_ctr.most_common(12)):
            axis_df = neuro_df[neuro_df['Neuro_Trap'].str.contains(re.escape(axis), na=False)]
            adhd_up = axis_df[axis_df['_direction'] == '↑ ADHD']
            ctrl_up = axis_df[axis_df['_direction'] == '↑ Control']
            n_ns    = cnt - len(adhd_up) - len(ctrl_up)
            med_fc  = axis_df[fc_col].median() if fc_col in axis_df.columns else np.nan
            with st.expander(f"**{axis}** — {cnt} compounds"):
                n1, n2, n3, n4 = st.columns(4)
                n1.markdown(f"<span class='badge-adhd'>↑ ADHD</span><br>"
                            f"<b style='font-size:1.4rem;color:#c0392b'>{len(adhd_up)}</b>",
                            unsafe_allow_html=True)
                n2.markdown(f"<span class='badge-ctrl'>↑ Control</span><br>"
                            f"<b style='font-size:1.4rem;color:#1a6ea8'>{len(ctrl_up)}</b>",
                            unsafe_allow_html=True)
                n3.markdown(f"<span style='background:#888;color:white;padding:3px 8px;"
                            f"border-radius:12px;font-size:0.85rem;font-weight:600'>n/s</span>"
                            f"<br><b style='font-size:1.4rem;color:#555'>{n_ns}</b>",
                            unsafe_allow_html=True)
                if np.isfinite(med_fc):
                    n4.metric("Median Log₂FC", f"{med_fc:+.2f}", help="Positive = higher in Control")

                ea, ec = st.columns(2)
                ea.markdown("<span class='badge-adhd'>↑ ADHD compounds</span>", unsafe_allow_html=True)
                ea.write(list(adhd_up[name_col].dropna().unique()))
                ec.markdown("<span class='badge-ctrl'>↑ Control compounds</span>", unsafe_allow_html=True)
                ec.write(list(ctrl_up[name_col].dropna().unique()))

                ax_ecols = export_cols(axis_df, name_col, fc_col, pval_col)
                ax_export = axis_df[ax_ecols].rename(
                    columns={name_col:'Compound', '_direction':'Direction',
                             fc_col:'Log2FC', pval_col:'p-value'})
                ax_key = re.sub(r'[^\w]', '_', axis)[:30]
                c1, c2, _ = st.columns([1,1,4])
                c1.download_button("⬇ CSV", to_csv(ax_export),
                                   file_name=f"neuro_{ax_key}.csv", mime='text/csv',
                                   key=f"csv_neuro_{_na_idx}", use_container_width=True)
                c2.download_button("⬇ Excel", to_excel(ax_export),
                                   file_name=f"neuro_{ax_key}.xlsx",
                                   mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                                   key=f"xl_neuro_{_na_idx}", use_container_width=True)

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 6 — ADHD vs CONTROL
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[5]:
        st.markdown("### Differential Analysis — ADHD vs Control")

        vol_df = df[[fc_col, pval_col, adjp_col, 'Original_Name', name_col]].copy()
        vol_df = vol_df.dropna(subset=[fc_col, pval_col])
        vol_df['-log10p'] = -np.log10(vol_df[pval_col].clip(lower=1e-20))
        p_s = (vol_df[adjp_col] if show_adj else vol_df[pval_col])
        vol_df['Significance'] = np.where(
            (vol_df[fc_col].abs() >= fc_thresh) & (p_s < pval_thresh),
            np.where(vol_df[fc_col] > 0, '↑ Control', '↑ ADHD'), 'n/s')

        color_map  = {'↑ Control': COLORS['Control'], '↑ ADHD': COLORS['ADHD'], 'n/s': '#cccccc'}
        sig_counts = vol_df['Significance'].value_counts()

        m1, m2, m3 = st.columns(3)
        m1.metric("↑ in ADHD",      sig_counts.get('↑ ADHD',    0))
        m2.metric("↑ in Control",   sig_counts.get('↑ Control', 0))
        m3.metric("Not significant", sig_counts.get('n/s',       0))

        fig_vol = px.scatter(
            vol_df, x=fc_col, y='-log10p',
            color='Significance', color_discrete_map=color_map,
            hover_data={name_col: True, fc_col: ':.2f',
                        '-log10p': ':.2f', pval_col: ':.4f', 'Significance': False},
            opacity=0.75,
            title=f'Volcano Plot  |  |FC| ≥ {fc_thresh}  |  p < {pval_thresh}'
                  f'{"  (FDR adjusted)" if show_adj else ""}',
            labels={fc_col: 'Log₂ Fold Change (Control/ADHD)', '-log10p': '–log₁₀(p)'},
        )
        fig_vol.add_vline(x= fc_thresh, line_dash='dash', line_color='#555', line_width=1)
        fig_vol.add_vline(x=-fc_thresh, line_dash='dash', line_color='#555', line_width=1)
        fig_vol.add_hline(y=-np.log10(pval_thresh), line_dash='dash', line_color='#555', line_width=1)
        xmax = vol_df[fc_col].abs().max() * 0.9
        fig_vol.add_annotation(x=-xmax, y=vol_df['-log10p'].max()*0.95,
                                text="← Elevated in ADHD", showarrow=False,
                                font=dict(size=13, color=COLORS['ADHD']))
        fig_vol.add_annotation(x= xmax, y=vol_df['-log10p'].max()*0.95,
                                text="Elevated in Control →", showarrow=False,
                                font=dict(size=13, color=COLORS['Control']))
        fig_vol.update_traces(marker_size=5)
        fig_vol.update_layout(height=540, margin=dict(t=60,b=20),
                               paper_bgcolor='white', plot_bgcolor='#fafafa',
                               font=dict(size=13, color='#111'))
        st.plotly_chart(fig_vol, use_container_width=True)

        sig_df = (vol_df[vol_df['Significance'] != 'n/s']
                  .merge(df[['Original_Name','Main_Categories','LM_Main_Class',
                              'KEGG_Pathways','Neuro_Trap']], on='Original_Name', how='left'))

        st.markdown("<div class='section-header'>Significant Compounds — Tables</div>",
                    unsafe_allow_html=True)
        tl, tr = st.columns(2)
        with tl:
            st.markdown("#### ↑ Higher in ADHD")
            up_adhd = sig_df[sig_df['Significance'] == '↑ ADHD'].sort_values(fc_col, ascending=True)
            n_show_a = st.slider("Show top N (ADHD-elevated)", 5, max(5,len(up_adhd)),
                                  min(25, len(up_adhd)), 5, key='n_adhd')
            min_fc_a = st.slider("|FC| minimum (ADHD)", 0.0, 5.0, 0.0, 0.25, key='fc_adhd')
            shown_a  = up_adhd[up_adhd[fc_col].abs() >= min_fc_a].head(n_show_a)
            st.dataframe(shown_a[[name_col, fc_col, pval_col, 'LM_Main_Class','KEGG_Pathways']]
                         .rename(columns={name_col:'Compound', fc_col:'Log2FC', pval_col:'p-value'})
                         .reset_index(drop=True), use_container_width=True, height=350)
        with tr:
            st.markdown("#### ↑ Higher in Control")
            up_ctrl = sig_df[sig_df['Significance'] == '↑ Control'].sort_values(fc_col, ascending=False)
            n_show_c = st.slider("Show top N (Control-elevated)", 5, max(5,len(up_ctrl)),
                                  min(25, len(up_ctrl)), 5, key='n_ctrl')
            min_fc_c = st.slider("|FC| minimum (Control)", 0.0, 5.0, 0.0, 0.25, key='fc_ctrl')
            shown_c  = up_ctrl[up_ctrl[fc_col].abs() >= min_fc_c].head(n_show_c)
            st.dataframe(shown_c[[name_col, fc_col, pval_col, 'LM_Main_Class','KEGG_Pathways']]
                         .rename(columns={name_col:'Compound', fc_col:'Log2FC', pval_col:'p-value'})
                         .reset_index(drop=True), use_container_width=True, height=350)

        # Downloads for significant compounds
        st.markdown("<div class='section-header'>Download Significant Compounds</div>",
                    unsafe_allow_html=True)
        all_sig = sig_df.rename(columns={name_col:'Compound', fc_col:'Log2FC',
                                          pval_col:'p-value', 'Significance':'Direction'})
        dl_col1, dl_col2 = st.columns(2)
        with dl_col1:
            adhd_export = up_adhd.rename(columns={name_col:'Compound', fc_col:'Log2FC', pval_col:'p-value'})
            st.download_button("⬇ ↑ADHD Compounds (CSV)", to_csv(adhd_export),
                               file_name="significant_ADHD_up.csv", mime='text/csv',
                               key='dl_adhd_sig', use_container_width=True)
        with dl_col2:
            ctrl_export = up_ctrl.rename(columns={name_col:'Compound', fc_col:'Log2FC', pval_col:'p-value'})
            st.download_button("⬇ ↑Control Compounds (CSV)", to_csv(ctrl_export),
                               file_name="significant_Control_up.csv", mime='text/csv',
                               key='dl_ctrl_sig', use_container_width=True)

        # Heatmap
        st.markdown("<div class='section-header'>Intensity Heatmap — Top Significant Compounds (log₁₀)</div>",
                    unsafe_allow_html=True)
        n_heat = st.slider("Compounds in heatmap", 10, 80, 40, 5)
        top_sig_names = (vol_df[vol_df['Significance'] != 'n/s']
                         .nlargest(n_heat, '-log10p')['Original_Name'].tolist())
        heat_df = df[df['Original_Name'].isin(top_sig_names)].drop_duplicates('Original_Name').copy()
        if not heat_df.empty:
            sorted_samp = adhd_cols + control_cols
            sorted_idx  = [sample_cols.index(s) for s in sorted_samp if s in sample_cols]
            heat_vals   = heat_df[sample_cols].apply(pd.to_numeric, errors='coerce').values
            heat_log    = np.log10(np.clip(heat_vals, 1, None))[:,sorted_idx]
            labels_x    = [f"{s} ({grp_norm.get(s,'?')})" for s in sorted_samp if s in sample_cols]
            labels_y    = heat_df[name_col].fillna(heat_df['Original_Name']).tolist()
            fig_heat    = go.Figure(go.Heatmap(
                z=heat_log, x=labels_x, y=labels_y, colorscale='RdBu_r',
                zmid=float(np.nanmedian(heat_log)),
                hovertemplate='%{y}<br>%{x}<br>log₁₀ intensity: %{z:.2f}<extra></extra>',
                colorbar=dict(title=dict(text='log₁₀(intensity)', font=dict(size=12))),
            ))
            fig_heat.add_vrect(x0=-0.5, x1=len(adhd_cols)-0.5,
                               fillcolor="rgba(192,57,43,0.06)", line_width=0)
            fig_heat.add_vrect(x0=len(adhd_cols)-0.5,
                               x1=len(adhd_cols)+len(control_cols)-0.5,
                               fillcolor="rgba(26,110,168,0.06)", line_width=0)
            fig_heat.update_layout(
                height=max(500, n_heat * 17), margin=dict(t=50,b=50),
                xaxis=dict(tickangle=45, tickfont_size=10), yaxis=dict(tickfont_size=11),
                title=dict(text=f'Top {n_heat} significant compounds | ADHD → Control', font_size=14),
                font=dict(size=12, color='#111'), paper_bgcolor='white',
            )
            st.plotly_chart(fig_heat, use_container_width=True)

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 7 — COMPOUND EXPLORER
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[6]:
        st.markdown("### Interactive Compound Explorer")
        cf1, cf2, cf3 = st.columns(3)
        search_term = cf1.text_input("🔍 Search compound name", placeholder="e.g. Arachidonic")
        cat_filter  = cf2.multiselect("Category", sorted(df['Main_Categories'].dropna().unique()))
        sig_filter  = cf3.selectbox("Direction filter",
                                     ['All','Significant only','↑ ADHD','↑ Control'])

        ex_df = df.copy()
        if search_term:
            mask = (ex_df[name_col].str.contains(search_term, case=False, na=False) |
                    ex_df['Original_Name'].str.contains(search_term, case=False, na=False))
            ex_df = ex_df[mask]
        if cat_filter:
            ex_df = ex_df[ex_df['Main_Categories'].apply(
                lambda x: any(c in str(x) for c in cat_filter) if pd.notna(x) else False)]
        if sig_filter == 'Significant only':
            ex_df = ex_df[ex_df['_direction'] != 'n/s']
        elif sig_filter in ('↑ ADHD','↑ Control'):
            ex_df = ex_df[ex_df['_direction'] == sig_filter]

        st.markdown(f"**{len(ex_df):,} compounds** match filters")

        disp_cols = [name_col,'Main_Categories','LM_Category','LM_Main_Class',
                     fc_col, pval_col,'OC_Ratio','KEGG_Pathways','Neuro_Trap','_direction']
        disp_cols = [c for c in disp_cols if c in ex_df.columns]
        show_df   = (ex_df[disp_cols]
                     .rename(columns={name_col:'Compound', fc_col:'Log2FC',
                                      pval_col:'p-value','_direction':'Direction'})
                     .reset_index(drop=True))
        st.dataframe(
            show_df, use_container_width=True, height=400,
            column_config={
                'Log2FC':   st.column_config.NumberColumn(format="%.2f"),
                'p-value':  st.column_config.NumberColumn(format="%.4f"),
                'OC_Ratio': st.column_config.ProgressColumn(min_value=0, max_value=1, format="%.2f"),
            }
        )

        # Download filtered results
        ex_ecols = export_cols(ex_df, name_col, fc_col, pval_col)
        ex_export = ex_df[ex_ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                      fc_col:'Log2FC', pval_col:'p-value'})
        c1, c2, _ = st.columns([1,1,4])
        c1.download_button("⬇ Download filtered (CSV)", to_csv(ex_export),
                           file_name="explorer_filtered.csv", mime='text/csv',
                           key='dl_explorer_csv', use_container_width=True)
        c2.download_button("⬇ Download filtered (Excel)", to_excel(ex_export),
                           file_name="explorer_filtered.xlsx",
                           mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                           key='dl_explorer_xl', use_container_width=True)

        # Compound Detail
        st.markdown("<div class='section-header'>Compound Detail View</div>", unsafe_allow_html=True)
        all_names = ex_df[name_col].dropna().unique()
        if len(all_names) == 0:
            st.info("No compounds match the current filters.")
        else:
            det_search = st.text_input("Filter compound list", placeholder="type to narrow...",
                                        key='det_search')
            if det_search:
                all_names = [n for n in all_names if det_search.lower() in n.lower()]
            if len(all_names) > 0:
                sel_name = st.selectbox(
                    f"Select compound ({len(all_names):,} available)",
                    options=all_names[:500],
                    format_func=lambda n: (n[:80] + '…') if len(n) > 80 else n,
                )
                row = ex_df[ex_df[name_col] == sel_name].iloc[0]
                d1, d2, d3 = st.columns(3)
                for col, pairs in [
                    (d1, [('Formula','Formula_Final'),('MW','MW_Final'),
                           ('KEGG ID','KEGG_ID_Final'),('HMDB','HMDB_ID')]),
                    (d2, [('Category','Main_Categories'),('LM Class','LM_Main_Class'),
                           ('O/C Ratio','OC_Ratio'),('NOSC','NOSC')]),
                    (d3, [('Log2FC',fc_col),('p-value',pval_col),
                           ('Neuro','Neuro_Trap'),('Clinical Tags','Clinical_Tags')]),
                ]:
                    for label, key in pairs:
                        val = row.get(key, 'N/A')
                        if isinstance(val, float):
                            val = f"{val:.4f}" if label in ('p-value','Log2FC','O/C Ratio') else f"{val:.3f}"
                        col.markdown(f"**{label}:** {val}")

                if sample_cols:
                    int_vals = pd.to_numeric(row[sample_cols], errors='coerce')
                    int_df = pd.DataFrame({
                        'Sample': sample_cols, 'Intensity': int_vals.values,
                        'Group': [grp_norm.get(s,'Unknown') for s in sample_cols],
                    }).dropna(subset=['Intensity'])
                    int_df = int_df[int_df['Intensity'] > 0].sort_values('Group')
                    fig_int = px.bar(int_df, x='Sample', y='Intensity', color='Group',
                                     color_discrete_map=COLORS, log_y=True,
                                     title=f"Per-sample intensity — {sel_name[:60]}",
                                     labels={'Intensity':'Intensity (log scale)'})
                    fig_int.update_layout(xaxis_tickangle=45, margin=dict(t=60,b=90))
                    apply_template(fig_int)
                    st.plotly_chart(fig_int, use_container_width=True)

                kegg_val = row.get('KEGG_Pathways','')
                if kegg_val and str(kegg_val) != '-':
                    st.markdown(f"**KEGG Pathways:** {kegg_val}")

    # ═════════════════════════════════════════════════════════════════════════
    # TAB 8 — DOWNLOADS & EXPORTS
    # ═════════════════════════════════════════════════════════════════════════
    with tabs[7]:
        st.markdown("### 📥 Downloads & Exports")
        st.markdown("All files are generated from the current dataset and statistical filters.")

        # ── Section 1: Full & Processed Files ────────────────────────────────
        st.markdown("<div class='section-header'>Full & Processed Datasets</div>",
                    unsafe_allow_html=True)

        full_ecols = export_cols(df, name_col, fc_col, pval_col)
        full_export = df[full_ecols].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                      fc_col:'Log2FC', pval_col:'p-value'})

        sig_mask = df['_direction'] != 'n/s'
        sig_export = df[sig_mask][full_ecols].rename(
            columns={name_col:'Compound', '_direction':'Direction',
                     fc_col:'Log2FC', pval_col:'p-value'})

        adhd_up_export = df[df['_direction'] == '↑ ADHD'][full_ecols].rename(
            columns={name_col:'Compound', '_direction':'Direction', fc_col:'Log2FC', pval_col:'p-value'})
        ctrl_up_export = df[df['_direction'] == '↑ Control'][full_ecols].rename(
            columns={name_col:'Compound', '_direction':'Direction', fc_col:'Log2FC', pval_col:'p-value'})

        r1c1, r1c2 = st.columns(2)
        with r1c1:
            st.markdown("<div class='dl-card'><div class='dl-card-title'>Full Annotated Dataset</div>"
                        f"<div class='dl-card-desc'>All compounds · annotation · statistics · classification"
                        f"<br><span class='dl-card-count'>{len(full_export):,} compounds</span></div></div>",
                        unsafe_allow_html=True)
            b1, b2 = st.columns(2)
            b1.download_button("⬇ CSV",   to_csv(full_export),   "full_dataset.csv",
                               mime='text/csv', key='dl_full_csv', use_container_width=True)
            b2.download_button("⬇ Excel", to_excel(full_export), "full_dataset.xlsx",
                               mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                               key='dl_full_xl', use_container_width=True)
        with r1c2:
            st.markdown("<div class='dl-card'><div class='dl-card-title'>Significant Compounds Only</div>"
                        f"<div class='dl-card-desc'>Compounds passing current FC & p-value thresholds"
                        f"<br><span class='dl-card-count'>{len(sig_export):,} compounds</span></div></div>",
                        unsafe_allow_html=True)
            b1, b2 = st.columns(2)
            b1.download_button("⬇ CSV",   to_csv(sig_export),   "significant_compounds.csv",
                               mime='text/csv', key='dl_sig_csv', use_container_width=True)
            b2.download_button("⬇ Excel", to_excel(sig_export), "significant_compounds.xlsx",
                               mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                               key='dl_sig_xl', use_container_width=True)

        r2c1, r2c2 = st.columns(2)
        with r2c1:
            st.markdown(f"<div class='dl-card'><div class='dl-card-title'>↑ Elevated in ADHD</div>"
                        f"<div class='dl-card-desc'>Compounds significantly higher in ADHD group"
                        f"<br><span class='dl-card-count'>{len(adhd_up_export):,} compounds</span></div></div>",
                        unsafe_allow_html=True)
            b1, b2 = st.columns(2)
            b1.download_button("⬇ CSV",   to_csv(adhd_up_export),   "elevated_ADHD.csv",
                               mime='text/csv', key='dl_adhd_up_csv', use_container_width=True)
            b2.download_button("⬇ Excel", to_excel(adhd_up_export), "elevated_ADHD.xlsx",
                               mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                               key='dl_adhd_up_xl', use_container_width=True)
        with r2c2:
            st.markdown(f"<div class='dl-card'><div class='dl-card-title'>↑ Elevated in Control</div>"
                        f"<div class='dl-card-desc'>Compounds significantly higher in Control group"
                        f"<br><span class='dl-card-count'>{len(ctrl_up_export):,} compounds</span></div></div>",
                        unsafe_allow_html=True)
            b1, b2 = st.columns(2)
            b1.download_button("⬇ CSV",   to_csv(ctrl_up_export),   "elevated_Control.csv",
                               mime='text/csv', key='dl_ctrl_up_csv', use_container_width=True)
            b2.download_button("⬇ Excel", to_excel(ctrl_up_export), "elevated_Control.xlsx",
                               mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                               key='dl_ctrl_up_xl', use_container_width=True)

        # ── Section 2: Thematic Slices ────────────────────────────────────────
        st.markdown("<div class='section-header'>Thematic Datasets</div>", unsafe_allow_html=True)

        lipid_slice  = df[df['LM_Category'].notna()]
        kegg_slice   = df[df['KEGG_Pathways'].notna() & (df['KEGG_Pathways'] != '-')]
        neuro_slice  = df[df['Neuro_Trap'].notna() & (df['Neuro_Trap'] != '-')]
        highox_slice = df[df['OC_Ratio'].notna() & (df['OC_Ratio'] > 0.6)]

        for label, desc, subset, key in [
            ("Classified Lipids (LIPID MAPS)",
             "All lipids with LM category, main class, and sub-class",
             lipid_slice, "slice_lipids"),
            ("KEGG-Mapped Compounds",
             "Compounds with at least one KEGG pathway annotation",
             kegg_slice, "slice_kegg"),
            ("Neuro-Active Compounds",
             "Compounds with Neuro_Trap annotation (any neurochemical axis)",
             neuro_slice, "slice_neuro"),
            ("High-Oxidation Compounds (O/C > 0.6)",
             "Compounds above the high oxidation threshold",
             highox_slice, "slice_highox"),
        ]:
            ecols_s = export_cols(subset, name_col, fc_col, pval_col)
            exp_s = subset[ecols_s].rename(columns={name_col:'Compound', '_direction':'Direction',
                                                     fc_col:'Log2FC', pval_col:'p-value'})
            st.markdown(f"<div class='dl-card'><div class='dl-card-title'>{label}</div>"
                        f"<div class='dl-card-desc'>{desc}"
                        f"<br><span class='dl-card-count'>{len(exp_s):,} compounds</span></div></div>",
                        unsafe_allow_html=True)
            b1, b2, _ = st.columns([1,1,5])
            b1.download_button("⬇ CSV",   to_csv(exp_s),   f"{key}.csv",
                               mime='text/csv', key=f'dl_{key}_csv', use_container_width=True)
            b2.download_button("⬇ Excel", to_excel(exp_s), f"{key}.xlsx",
                               mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                               key=f'dl_{key}_xl', use_container_width=True)

        # ── Section 3: MetaboAnalyst Export ──────────────────────────────────
        st.markdown("<div class='section-header'>MetaboAnalyst Export</div>", unsafe_allow_html=True)

        peak_table, meta_df, n_included, n_excluded = build_metaboanalyst(df, sample_cols, grp_norm)

        st.markdown(f"""
        <div class='export-box'>
            <div class='export-box-title'>MetaboAnalyst Ready Files</div>
            <div class='export-box-desc'>
                Two files are required for MetaboAnalyst statistical analysis:<br>
                <b>Peak Table</b> — rows = compounds (with PubChem/HMDB/KEGG ID), columns = samples<br>
                <b>Metadata</b> — sample names with group labels
            </div>
            <span style='color:#27ae60;font-weight:600'>✓ {n_included:,} compounds included</span>
            &nbsp;&nbsp;
            <span style='color:#c0392b;font-weight:600'>✗ {n_excluded:,} excluded (no recognized ID)</span>
            <br><small style='color:#777'>ID priority: PubChem CID → HMDB → KEGG. Compounds without any of these are excluded.</small>
        </div>
        """, unsafe_allow_html=True)

        ma1, ma2 = st.columns(2)
        with ma1:
            st.markdown("**File 1 — Peak Table**")
            st.caption(f"Rows: {n_included:,} compounds · Columns: {len(sample_cols)} samples")
            st.download_button(
                "⬇ Download Peak Table (CSV)",
                to_csv(peak_table),
                file_name="metaboanalyst_peak_table.csv",
                mime='text/csv',
                key='dl_ma_peak', use_container_width=True)
            with st.expander("Preview (first 5 rows)"):
                st.dataframe(peak_table.head(), use_container_width=True)

        with ma2:
            st.markdown("**File 2 — Sample Metadata**")
            st.caption(f"{len(meta_df)} samples · 2 columns: Sample, Group")
            st.download_button(
                "⬇ Download Metadata (CSV)",
                to_csv(meta_df),
                file_name="metaboanalyst_metadata.csv",
                mime='text/csv',
                key='dl_ma_meta', use_container_width=True)
            with st.expander("Preview"):
                st.dataframe(meta_df, use_container_width=True)

        # ── Section 4: With Intensity Values ─────────────────────────────────
        st.markdown("<div class='section-header'>Full Data With Sample Intensities</div>",
                    unsafe_allow_html=True)
        st.markdown("Includes all annotation columns **plus** raw intensity per sample. "
                    "Large file — use for downstream R/Python analysis.")

        intensity_export = df[full_ecols + sample_cols].rename(
            columns={name_col:'Compound', '_direction':'Direction',
                     fc_col:'Log2FC', pval_col:'p-value'})
        b1, b2, _ = st.columns([1,1,5])
        b1.download_button("⬇ CSV (with intensities)", to_csv(intensity_export),
                           file_name="full_with_intensities.csv", mime='text/csv',
                           key='dl_int_csv', use_container_width=True)
        b2.download_button("⬇ Excel (with intensities)", to_excel(intensity_export),
                           file_name="full_with_intensities.xlsx",
                           mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                           key='dl_int_xl', use_container_width=True)


# ─── Entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
