#!/usr/bin/env python
# coding: utf-8

# This method was first introduced by [Farlik, Halbritter, et al.](https://doi.org/10.1016/j.stem.2016.10.019) to computationally reconstruct the human hematopoietic lineage from **DNA methylation** profiles. It was later adapted by [Traxler, Reichl, et al.](https://doi.org/10.1016/j.cels.2025.101346) to infer functional relationships between gene knockouts from **scCRISPR-seq** perturbation signatures, successfully validated by identifying members of known protein complexes.

####libraries
import os
import pathlib
import pandas as pd
import numpy as np
# for classification
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
# for visualization
import networkx as nx
import matplotlib.pyplot as plt
# pydot is used as the interface to Graphviz
import pydot

#### configs

# input
data_path = snakemake.input["data"]
metadata_path = snakemake.input["metadata"]
feature_annotation_path = snakemake.input["feature_annotation"]

# output
adjacency_matrix_path = snakemake.output["adjacency_matrix"]
top_features_path =  snakemake.output["top_features"]
graph_path = snakemake.output["graph"]

# params
group_var = snakemake.params["group_var"]
top_features_n = snakemake.params["top_features_n"]
prune_th = snakemake.params["prune_th"]
feature_annotation_var = snakemake.params["feature_annotation_var"]

### Load & prepare data and metadata
data = pd.read_csv(data_path, index_col=0, header=0)
metadata = pd.read_csv(metadata_path, index_col=0, header=0)
feature_annotation = pd.read_csv(feature_annotation_path, index_col=0, header=0)

# sort them the same
data = data.loc[:,metadata.index]

# Prepare data for training
X = np.array(data.T)
X = X-X.mean(axis=0) # center data
y = metadata[group_var]

#### Prepare Classifier (Logistic Regression w/ mostly default parameters)
clf = LogisticRegression(
    penalty="elasticnet",
    solver="saga",
    multi_class="multinomial",
    max_iter=100, # default 100
    n_jobs=-1,
    random_state=42,
    verbose=1,
    l1_ratio=0.5
)

#################################################################
#### Get connectivity matrix from LOO-CV strategy prediction probabilities
#################################################################

classnames, groups = np.unique(y, return_inverse=True)
cv = sklearn.model_selection.LeaveOneGroupOut()
pred = sklearn.model_selection.cross_val_predict(estimator=clf, X=X, y=y, groups=groups, cv=cv, n_jobs=-1, method='predict_proba')

# sanity check if approach works as intended -> pred. prob for correct label has to be 0!
result=pd.DataFrame(pred, columns=classnames, index=y)
for index, row in result.iterrows():
    if row[index]!=0:
        print('Something went wrong. The prediction probabilities of the correct class are not 0.')

# aggregate the prediction probabilities using arithmetic mean
conn_norm = pd.DataFrame(columns=classnames, index=classnames)
for col in classnames:
    for row in classnames:
        conn_norm.loc[row,col]=result.loc[row,col].mean()

# sanity check if approach works as intended -> diagonal has to be 0
if sum(np.diag(conn_norm)!=0)>0:
    print('Something went wrong. The average prediction probabilities of the correct class are not 0.')

# save normalized probability connectivity matrix
conn_norm.to_csv(adjacency_matrix_path)

#################################################################
#### Train Complete Model and Extract Top Features per Class
#################################################################

# Train the classifier on the entire dataset
clf.fit(X, y)

# Get top features (i.e., coefficients by importances) for each class  
top_features_per_class = {}

# Iterate over each class from the fitted model
for i, class_name in enumerate(clf.classes_):
    # Get coefficients for the current class. For multinomial logistic regression, coef_ has shape (n_classes, n_features)
    class_coeffs = clf.coef_[i]

    # Get feature names of the {top_features_n} features based on coefficient value
    # Use np.argsort to get indices that would sort the array, then take the last {top_features_n} i.e., largest
    top_features_per_class[class_name] = data.index[np.argsort(class_coeffs)[-top_features_n:]][::-1] # Reverse to show most important first

# Create a mapping dictionary
feature_name_map = feature_annotation[feature_annotation_var].to_dict()
# Iterate through the original dictionary and map features to annotation
top_features_mapped = {}
for class_name, features in top_features_per_class.items():
    # For each list of feature IDs, create a new list with the corresponding alternative names.
    # The .get(feature, feature) method ensures that if an ID is not found in the map, the original ID is used instead.
    top_features_mapped[class_name] = [feature_name_map.get(feature, feature) for feature in features]

# save top features
pd.DataFrame(top_features_mapped).to_csv(top_features_path)

#################################################################
#### HIERARCHICAL VISUALIZATION WITH GRAPHVIZ ####
#################################################################

# --- Pruning the Adjacency Matrix ---
# conn_viz = conn_norm.copy()
conn_viz = pd.read_csv(adjacency_matrix_path, index_col=0, header=0)
conn_viz[conn_viz < prune_th] = 0

# Create a directed graph from the PRUNED adjacency matrix
G = nx.from_pandas_adjacency(conn_viz, create_using=nx.DiGraph())

# --- Visualization Setup ---
# A wider figure is better for hierarchical layouts
plt.figure(figsize=(4, 4))

# Use graphviz_layout to create a hierarchical ("dot") layout.
# - 'prog="dot"' specifies the hierarchical layout engine.
# - The engine automatically places nodes with an in-degree of 0 (like HSC) at the top.
# - It also handles disconnected components by plotting them separately.
pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='dot')

# --- Drawing the Graph ---
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos, node_size=300, node_color='skyblue')
nx.draw_networkx_edges(
    G,
    pos,
    width=[w * 4 for w in edge_weights],
    edge_color='grey',
    node_size=300,
    arrowstyle='->',
    arrowsize=15,
    connectionstyle='arc3,rad=0.2' # Use slight curves for clarity
)

# Draw labels
nx.draw_networkx_labels(
    G,
    pos,
    font_size=6,
    font_family='sans-serif'
)

# save plot
plt.axis('off')
plt.savefig(graph_path, dpi=300, bbox_inches='tight')
# plt.show()