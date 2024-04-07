import os
import base64
import streamlit as st

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

st.set_page_config(
    page_title="BioNets - Home",
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

# ---------- Banner ----------
banner_file     = open(os.path.join(BASE_DIR, 'img', 'banner2.gif'), "rb")
banner_contents = banner_file.read()
data_url        = base64.b64encode(banner_contents).decode("utf-8")
banner_file.close()

banner_style = """ 
<style>
  .container {
    display: flex;
    justify-content: center;
  }
  .image {
    max-width: 100%;
    height: auto;
  }
</style>
"""

banner_html = banner_style + f"""
<div class="container">
  <div style="width: 20%;"></div>
  <div style="width: 70%;">
    <img class="image" src="data:image/gif;base64,{data_url}" alt="banner" width="900" height="450">
  </div>
  <div style="width: 20%;"></div>
</div>
"""

st.markdown(banner_html, unsafe_allow_html=True)
# ---------- Banner ----------

# ---------- Content ----------
st.markdown('''
Welcome to the Biomedical Network Science Project for the Winter'23/24 semester at [Friedrich-Alexander University Erlangen-Nürnberg](https://www.fau.eu/)!

Our project called **"Large-scale Differential Gene Expression Analysis in scRNA-seq Data"**, proposed by Biomedical Network Science ([BIONETS](https://www.bionets.tf.fau.de/)) lab, supervised by [Prof. Dr. David B. Blumenthal](https://www.bionets.tf.fau.de/person/david-b-blumenthal/), and [Dr. Anne Hartebrodt](https://www.bionets.tf.fau.de/person/anne-hartebrodt/) at FAU Erlangen-Nürnberg.
''')

st.subheader("Project Goals")
st.markdown("The project aims to retrieve suitable scRNA-seq data for a specific disease, implement pipelines for standard preprocssing, clustering, automatic cell type annotation and DEG identification, with plans to create an interactive dashboard and extend to multiple diseases.")

st.subheader("Project Description")
st.markdown('''
1. **Dataset:** Two data sources of two diseases, such as [type II Diabetes Mellitus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515) and [Myeloproliferative Neoplasm (MPN)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244589) have been obtained from the [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/).
2. **Pipelines:** The overall pipeline consists of three key stages: _preprocessing adata_, _automatic cell-type annotation_, and _differential gene expression (DGE) analysis_.
    - **Preprocessing adata:** Standard preprocessing has been done, including quantity control, normalization, dimensionality reduction using PCA, and etc.
    - **Automatic Cell-type Annotation:** In this stage, automatic cell type annotation has been performed using two methods: [SCSA](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00490/full) and [MetaTiME](https://www.nature.com/articles/s41467-023-38333-8).
    - **DGE Analysis:** The differential gene expression (DGE) analysis pipeline performs gene expression analysis using various methods such as _**t-test, wilcoxon rank-sum, logistic regression, and t-test with overestimated variance**_.
3. **Dashboard:** Built an interactive dashboard using [Streamlit](https://streamlit.io). And you are here! :crossed_fingers:
''')

st.subheader("Contributors")
st.markdown('''
- [Farzam Taghipour](https://www.linkedin.com/in/farzamtaghipour/), Graduate Student in Artificial Intelligence at FAU Erlangen-Nürnberg
- [Sujit Debnath](https://www.linkedin.com/in/sujit-debnath/), Graduate Student in Artificial Intelligence at FAU Erlangen-Nürnberg
''')
# ---------- Content ----------