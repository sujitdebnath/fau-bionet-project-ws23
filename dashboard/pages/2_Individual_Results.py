import os
import pandas as pd
import streamlit as st
from PIL import Image

BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')

def resize_image(img_path, ratio):
    img     = Image.open(img_path)
    new_img = img.resize((int(img.size[0]*ratio), int(img.size[1]*ratio)))

    return new_img

st.set_page_config(
    page_title="BioNets - Results",
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

# ---------- Drop-down ----------
col1, col2 = st.columns(2)
with col1:
    # Select Disease
    disease_ids = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))])
    disease_id = st.selectbox('Select Disease', disease_ids)
with col2:
    # Select Dataset
    disease_dir = os.path.join(BASE_RES_DIR, disease_id)
    dataset_ids = sorted([dir for dir in os.listdir(disease_dir)])
    dataset_id = st.selectbox('Select Dataset', dataset_ids)
# ---------- Drop-down ----------

# ---------- Submit Button ----------
_, _, _, _, col5 = st.columns(5)
with col5:
    click_submit = st.button('Submit', use_container_width=True)
# ---------- Submit Button ----------

# ---------- Submit Button Mechanism ----------
if click_submit:
    results_dir = os.path.join(BASE_RES_DIR, disease_id, dataset_id)
    st.markdown(f'<h4 style="text-align: center; color: black;">Results for {dataset_id.capitalize()} of {disease_id.upper()} Disease</h1>', unsafe_allow_html=True)

    st.markdown("#### 1. Donor (Case vs Control) and Clustering")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir, 'donor_cells.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir, 'leiden.png'), ratio=0.8))

    st.markdown("#### 2. Cell Annotation - SCSA")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir, 'scsa_cellmarker.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir, 'scsa_panglaodb.png'), ratio=0.6))
    
    st.markdown("#### 3. Cell Annotation - MetaTiME and Ro/e Viz (if available)")
    st.image(resize_image(img_path=os.path.join(results_dir, 'metatime_minor.png'), ratio=0.7))

    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir, 'metatime_major.png'), ratio=0.8))
    
    roe_img = os.path.join(results_dir, 'roe.png')
    if os.path.exists(roe_img):
        with col2:
            st.image(resize_image(img_path=roe_img, ratio=0.75))
    
    dge_results = pd.read_csv(os.path.join(results_dir, 'dge_analysis_result.csv'))
    st.markdown("#### 4. Differential Gene Expression Analysis Results")
    st.dataframe(dge_results, height=1010, use_container_width=True, hide_index=True)
# ---------- Submit Button Mechanism ----------