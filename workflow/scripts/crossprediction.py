#!/usr/bin/env python
# coding: utf-8

####libraries
import os
import pathlib
import pandas as pd
import numpy as np
# for classification
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV

#### configs

# input
data_path = snakemake.input["data"]
metadata_path = snakemake.input["metadata"]

# output
adjacency_matrix_path = snakemake.output["adjacency_matrix"]
top_features_path =  snakemake.output["top_features"]

# params
group_var = snakemake.params["group_var"]
top_features_n = 5

### Load & prepare data and metadata
data = pd.read_csv(data_path, index_col=0, header=0)
metadata = pd.read_csv(metadata_path, index_col=0, header=0)

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

#### get connectivity matrix from LOO-CV strategy prediction probabilities
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

#### Train Complete Model and Extract Top Features per Class

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

# Print results
print("Top {} features per class:".format(top_features_n))
for class_name, features in top_features_per_class.items():
    print(f"\nClass: {class_name}")
    for j, feature in enumerate(features):
        print(f"  {j+1}. {feature}")

# save top features
pd.DataFrame(top_features_per_class).to_csv(top_features_path)
