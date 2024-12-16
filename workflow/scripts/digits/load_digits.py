#### laod the digits dataset from sklearn ####

#### libraries
from sklearn.datasets import load_digits
import pandas as pd

# outputs
data_path = snakemake.output["data"]
labels_path = snakemake.output["labels"]

# load the digits dataset
digits = load_digits()

# save to CSV with indices as the first column
pd.DataFrame(digits.data, columns=digits.feature_names).to_csv(data_path, index=True)
pd.DataFrame(digits.target, columns=['target']).to_csv(labels_path, index=True)
