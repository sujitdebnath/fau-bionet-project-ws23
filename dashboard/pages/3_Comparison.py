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
st.markdown("#### Choose 1st dataset to compare:")
col1, col2 = st.columns(2)
with col1:
    # Select Disease
    disease_ids1 = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))])
    disease_id1 = st.selectbox('Select Disease', disease_ids1, key='1')
with col2:
    # Select Dataset
    disease_dir1 = os.path.join(BASE_RES_DIR, disease_id1)
    dataset_ids1 = sorted([dir for dir in os.listdir(disease_dir1)])
    dataset_id1 = st.selectbox('Select Dataset', dataset_ids1, key='2')
# ---------- First Dataset: Drop-down ----------

# ---------- Second Dataset: Drop-down ----------
st.markdown("#### Choose 2nd dataset to compare:")
col1, col2 = st.columns(2)
with col1:
    # Select Disease
    disease_ids2 = sorted([dir for dir in os.listdir(BASE_RES_DIR) if os.path.isdir(os.path.join(BASE_RES_DIR, dir))])
    disease_id2 = st.selectbox('Select Disease', disease_ids2, key='3')
with col2:
    # Select Dataset
    disease_dir2 = os.path.join(BASE_RES_DIR, disease_id2)
    dataset_ids2 = sorted([dir for dir in os.listdir(disease_dir2)])

    if disease_id1 == disease_id2:
        dataset_ids2.remove(dataset_id1)
    
    dataset_id2 = st.selectbox('Select Dataset', dataset_ids2, key='4')
# ---------- Second Dataset: Drop-down ----------

# ---------- Submit Button ----------
_, _, _, _, col5 = st.columns(5)
with col5:
    click_submit = st.button('Submit', use_container_width=True)
# ---------- Submit Button ----------