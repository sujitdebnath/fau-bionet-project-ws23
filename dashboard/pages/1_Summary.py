import os
import streamlit as st

BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')

st.set_page_config(
    page_title="BioNets - Summary",
    page_icon=":dna:",
    layout='wide',
    initial_sidebar_state='expanded'
)

# ----- Hide the streamlit header and footer -----
hide_steamlit_style = """
<style>
#MainMenu {visibility: hidden;}
.stDeployButton {display: none;}
footer {visibility: hidden;}
</style>
"""
st.markdown(hide_steamlit_style, unsafe_allow_html=True)
# ----- Hide the streamlit header and footer -----

# ---------- Title ----------
st.markdown('<h1 style="text-align: center; color: black;">🧬 Project BioNets Dashboard</h1>', unsafe_allow_html=True)

# ---------- Metrices ----------
st.markdown("### Metrics")

st.markdown(
    """
    <style>
        /* Card */
        div.metric-card {
            background-color: #FFFFFF;
            border: 1px solid #CCCCCC;
            padding: 5% 5% 5% 10%;
            border-radius: 5px;
            border-left: 0.5rem solid Grey;
            box-shadow: 0 0.15rem 1.75rem 0 rgba(58, 59, 69, 0.15);
            margin-bottom: 1rem;
        }

        /* Metric Label */
        label.metric-label {
            color: Grey;
            font-weight: 700;
            text-transform: uppercase;
        }
    </style>
    """,
    unsafe_allow_html=True
)

col1, col2, col3, col4 = st.columns(4)

disease_ids = [dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))]
dataset_ids = [dataset_id for disease_id in disease_ids for dataset_id in os.listdir(os.path.join(BASE_RES_DIR, disease_id))]

# Metrics in card-like design
with col1:
    st.markdown('<div class="metric-card"><label class="metric-label">Diseases</label><br><span>' + str(len(disease_ids)) + '</span></div>', unsafe_allow_html=True)
with col2:
    st.markdown('<div class="metric-card"><label class="metric-label">Datasets</label><br><span>' + str(len(dataset_ids)) + '</span></div>', unsafe_allow_html=True)
with col3:
    st.markdown('<div class="metric-card"><label class="metric-label">Annotation Methods</label><br><span>2</span></div>', unsafe_allow_html=True)
with col4:
    st.markdown('<div class="metric-card"><label class="metric-label">DGE Methods</label><br><span>4</span></div>', unsafe_allow_html=True)
# ---------- Metrices ----------

# ---------- Summary of Cell Annotation ----------
st.markdown("### Summary of Cell type Annotation")
# ---------- Summary of Cell Annotation ----------

# ---------- Summary of DGE Analysis ----------
st.markdown("### Summary of DGE Analysis")
# ---------- Summary of DGE Analysis ----------