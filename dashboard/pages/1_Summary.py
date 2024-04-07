import os
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go

BASE_DIR      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BASE_RES_DIR  = os.path.join(BASE_DIR, 'results')

def read_csv(csv_fpath):
    return pd.read_csv(csv_fpath)

def map_cell_anno_method_to_col_name(cell_anno_method):
    if cell_anno_method == "SCSA - cellmarker":
        cell_anno_col_name = 'scsa_celltype_cellmarker'
    elif cell_anno_method == "SCSA - panglaodb":
        cell_anno_col_name = 'scsa_celltype_panglaodb'
    elif cell_anno_method == "MetaTiME":
        cell_anno_col_name = 'Major_MetaTiME'
    
    return cell_anno_col_name

def cell_anno_box_plot(df, cell_anno_method):
    cell_anno_col_name = map_cell_anno_method_to_col_name(cell_anno_method)
    
    # Aggregate counts of cell types for each disease and dataset
    agg_df = df.groupby(['disease_id', 'dataset_id', cell_anno_col_name]).size().reset_index(name='count')

    # Filter out rows with 'Unknown' cell types
    agg_df = agg_df[agg_df[cell_anno_col_name] != 'Unknown']

    # Calculate total count of each cell type across all datasets and diseases
    total_counts = agg_df.groupby(cell_anno_col_name)['count'].sum().reset_index()
    total_counts = total_counts.sort_values(by='count', ascending=False)  # Sort by count

    # Sort agg_df based on total count order
    agg_df[cell_anno_col_name] = pd.Categorical(agg_df[cell_anno_col_name], categories=total_counts[cell_anno_col_name], ordered=True)
    agg_df = agg_df.sort_values(by=cell_anno_col_name)

    # Create a Plotly box plot
    fig = px.box(
        agg_df,
        x=cell_anno_col_name,
        y='count',
        color='disease_id',
        category_orders={cell_anno_col_name: list(total_counts[cell_anno_col_name].values),
                         'disease_id': sorted(df['disease_id'].unique())},
        title=f'Distribution of Cell Types Across Datasets for Each Disease ({cell_anno_method})',
        labels={'count': 'Count', cell_anno_col_name: 'Cell Type'},
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_traces(boxmean='sd')  # Show mean and standard deviation on the box plot

    return fig

def cell_anno_box_plot2(df, cell_anno_method):
    cell_anno_col_name = map_cell_anno_method_to_col_name(cell_anno_method)
    
    # Set the style
    sns.set_style("whitegrid")

    # Aggregate counts of cell types for each disease and dataset
    agg_df = df.groupby(['disease_id', 'dataset_id', cell_anno_col_name]).size().reset_index(name='count')
    agg_df = agg_df[agg_df[cell_anno_col_name] != 'Unknown'] # exclude 'Unknown' labels

    # Calculate total count of each cell type across all datasets and diseases
    total_counts = agg_df.groupby(cell_anno_col_name)['count'].sum().reset_index()
    total_counts = total_counts.sort_values(by='count', ascending=False)  # Sort by count

    # Sort agg_df based on total count order
    agg_df[cell_anno_col_name] = pd.Categorical(agg_df[cell_anno_col_name], categories=total_counts[cell_anno_col_name], ordered=True)
    agg_df = agg_df.sort_values(by=cell_anno_col_name)

    # Plotting
    fig = plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=agg_df,
        x=cell_anno_col_name,
        y='count',
        hue='disease_id',
        hue_order=sorted(df['disease_id'].unique()),
        palette='Set2'
    )
    plt.title(f'Distribution of Cell Types Across Datasets for Each Disease ({cell_anno_method})')
    plt.xlabel('Cell Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    return fig

import plotly.graph_objects as go

def cell_anno_heatmap(df, disease_id, anno_method1, anno_method2):
    cell_anno_col1 = map_cell_anno_method_to_col_name(anno_method1)
    cell_anno_col2 = map_cell_anno_method_to_col_name(anno_method2)

    # Filter data by disease_id
    df_filtered = df[df['disease_id'] == disease_id]

    # Compute the co-occurrence matrix
    co_occurrence_matrix = df_filtered.groupby([cell_anno_col1, cell_anno_col2]).size().unstack(fill_value=0)

    # Create a Plotly heatmap
    fig = go.Figure(data=go.Heatmap(
        z=co_occurrence_matrix.values,
        x=co_occurrence_matrix.columns,
        y=co_occurrence_matrix.index,
        colorscale='Viridis',  # Choose a color scale
        hoverongaps=False  # Disable hover on gaps
    ))

    # Customize layout
    fig.update_layout(
        title=f"Cell Type Co-occurrence Heatmap for {disease_id.upper()} Disease",
        xaxis_title=anno_method2,
        yaxis_title=anno_method1,
        xaxis=dict(side="top"),  # Display x-axis on top
        height=800,
        width=800
    )

    return fig

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
st.markdown("#### Metrics")

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
st.markdown("#### 1. Summary of Cell type Annotation")

# ---------- Summary 1 ----------
st.markdown("##### 1.1. Distribution of Cell Types (unknown labels excluded)")
cell_anno_methods = ['SCSA - cellmarker', 'SCSA - panglaodb', 'MetaTiME']
cell_anno_methods_box = st.selectbox('Select Cell-type Annotation Method', cell_anno_methods, key='1')

cell_anno_df = read_csv(csv_fpath=os.path.join(BASE_DIR, 'dashboard', 'cell_anno_res.csv'))
# st.pyplot(cell_anno_box_plot2(df=cell_anno_df, cell_anno_method=cell_anno_methods_box))
st.plotly_chart(cell_anno_box_plot(df=cell_anno_df, cell_anno_method=cell_anno_methods_box), use_container_width=True)
# ---------- Summary 1 ----------

# ---------- Summary 2 ----------
st.markdown("##### 1.2. Heatmap of Cell Type Co-occurrence")

col1, col2, col3 = st.columns(3)
with col1:
    disease_id = st.selectbox('Select Disease', sorted(disease_ids))
with col2:
    cell_anno_methods_box1 = st.selectbox('Select 1st Annotation Method', cell_anno_methods, key='2')
with col3:
    cell_anno_methods2     = [anno_method for anno_method in cell_anno_methods if anno_method != cell_anno_methods_box1]
    cell_anno_methods_box2 = st.selectbox('Select 2nd Annotation Method', cell_anno_methods2, key='3')

st.plotly_chart(cell_anno_heatmap(
                    df=cell_anno_df,
                    disease_id=disease_id,
                    anno_method1=cell_anno_methods_box1,
                    anno_method2=cell_anno_methods_box2),
                use_container_width=True)
# ---------- Summary 2 ----------
# ---------- Summary of Cell Annotation ----------

# ---------- Summary of DGE Analysis ----------
# st.markdown("#### 2. Summary of DGE Analysis")
# ---------- Summary of DGE Analysis ----------