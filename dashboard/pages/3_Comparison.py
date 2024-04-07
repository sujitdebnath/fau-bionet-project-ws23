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
    page_title="BioNets - Comparison",
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

# ---------- First Dataset: Drop-down ----------
st.markdown("##### Choose disease to compare:")
# Select Disease
disease_ids = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))])
disease_id  = st.selectbox('Select Disease', disease_ids)
# ---------- First Dataset: Drop-down ----------

# ---------- Second Dataset: Drop-down ----------
disease_dir = os.path.join(BASE_RES_DIR, disease_id)
dataset_ids = sorted([dir for dir in os.listdir(disease_dir)])

st.markdown("##### Choose two datasets to compare:")
col1, col2 = st.columns(2)
with col1:
    # Select Dataset1
    dataset_id1  = st.selectbox('Select Dataset1', dataset_ids)
with col2:
    # Select Dataset2
    dataset_ids.remove(dataset_id1)
    dataset_id2 = st.selectbox('Select Dataset2', dataset_ids)
# ---------- Second Dataset: Drop-down ----------

# ---------- Submit Button ----------
_, _, _, _, col5 = st.columns(5)
with col5:
    click_submit = st.button('Submit', use_container_width=True)
# ---------- Submit Button ----------

if click_submit:
    results_dir1 = os.path.join(BASE_RES_DIR, disease_id, dataset_id1)
    results_dir2 = os.path.join(BASE_RES_DIR, disease_id, dataset_id2)
    st.markdown(f'<h4 style="text-align: center; color: black;">Comparison between {dataset_id1.capitalize()} and {dataset_id2.capitalize()} of {disease_id.upper()} Disease (Left: {dataset_id1.capitalize()}, Right: {dataset_id2.capitalize()})</h1>', unsafe_allow_html=True)

    # ---------- Comparison (Row 1) ----------
    st.markdown("#### 1. Donor (Case vs Control)")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir1, 'donor_cells.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir2, 'donor_cells.png'), ratio=0.8))
    # ---------- Comparison (Row 1) ----------

    # ---------- Comparison (Row 2) ----------
    st.markdown("#### 2. Leiden Clustering")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir1, 'leiden.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir2, 'leiden.png'), ratio=0.8))
    # ---------- Comparison (Row 2) ----------

    # ---------- Comparison (Row 3) ----------
    st.markdown("#### 3. Cell Annotation - SCSA")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir1, 'scsa_cellmarker.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir2, 'scsa_cellmarker.png'), ratio=0.8))
    # ---------- Comparison (Row 3) ----------

    # ---------- Comparison (Row 4) ----------
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir1, 'scsa_panglaodb.png'), ratio=0.6))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir2, 'scsa_panglaodb.png'), ratio=0.6))
    # ---------- Comparison (Row 4) ----------

    # ---------- Comparison (Row 5) ----------
    st.markdown("#### 4. Cell Annotation - MetaTiME")
    col1, col2 = st.columns(2)
    with col1:
        st.image(resize_image(img_path=os.path.join(results_dir1, 'metatime_major.png'), ratio=0.8))
    with col2:
        st.image(resize_image(img_path=os.path.join(results_dir2, 'metatime_major.png'), ratio=0.8))
    # ---------- Comparison (Row 5) ----------

    # ---------- Comparison (Row 6) ----------
    roe_img1 = os.path.join(results_dir1, 'roe.png')
    roe_img2 = os.path.join(results_dir1, 'roe.png')
    roe_comp = False

    if os.path.exists(roe_img1) and os.path.exists(roe_img2):
        st.markdown("#### 5. Ratio of Observed to Expected cell numbers (Ro/e)")
        col1, col2 = st.columns(2)
        with col1:
            st.image(resize_image(img_path=roe_img1, ratio=0.75))
        with col2:
            st.image(resize_image(img_path=roe_img2, ratio=0.75))
        
        roe_comp = True
    # ---------- Comparison (Row 6) ----------

    # ---------- Comparison (Row 7) ----------
    dge_results1 = pd.read_csv(os.path.join(results_dir1, 'dge_analysis_result.csv'))
    dge_results2 = pd.read_csv(os.path.join(results_dir2, 'dge_analysis_result.csv'))

    if roe_comp:
        st.markdown("#### 6. Differential Gene Expression Analysis Results")
    else:
        st.markdown("#### 5. Differential Gene Expression Analysis Results")
    
    col1, col2 = st.columns(2)
    with col1:
        st.dataframe(dge_results1, height=1010, use_container_width=True, hide_index=True)
    with col2:
        st.dataframe(dge_results2, height=1010, use_container_width=True, hide_index=True)
    # ---------- Comparison (Row 7) ----------