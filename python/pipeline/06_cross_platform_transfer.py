import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import root_mean_squared_error, r2_score

# 1. Load the merged dataset
file_path = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
df = pd.read_csv(file_path)

hania_paper_id = 'Hanna_2025'  # Exact string for Hania's dataset

# 2. Convert Bacteria to Relative Abundance (Row-wise normalization)
genus_cols = [col for col in df.columns if col.startswith('g__')]
df[genus_cols] = df[genus_cols].fillna(0)

# Divide each cell by its sample's total genus reads
row_sums = df[genus_cols].sum(axis=1)
# Prevent division by zero if a sample has 0 total reads
row_sums = row_sums.replace(0, 1) 
df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

# 3. Hard Split: Train on Illumina (Master), Test on Nanopore (Hania)
train_df = df[df['Paper_ID'] != hania_paper_id].copy()
test_df = df[df['Paper_ID'] == hania_paper_id].copy()

# 4. Define Features and Target
TARGET = 'Caproate'
op_and_substrates = ['PH', 'TEMP', 'HRT', 'Lactate', 'Acetate', 'Ethanol']

# Intersect to keep only features present in the dataframe
features = genus_cols + [col for col in op_and_substrates if col in df.columns]

X_train = train_df[features].apply(pd.to_numeric, errors='coerce').fillna(0)
y_train = pd.to_numeric(train_df[TARGET], errors='coerce').fillna(0)

X_test = test_df[features].apply(pd.to_numeric, errors='coerce').fillna(0)
y_test = pd.to_numeric(test_df[TARGET], errors='coerce').fillna(0)

# 5. Train on Short-Read, Test on Long-Read
model = xgb.XGBRegressor(
    n_estimators=150, 
    learning_rate=0.05, 
    max_depth=3, 
    random_state=42
).fit(X_train, y_train)

y_pred = model.predict(X_test)

# 6. Print Results
print("--- RELATIVE ABUNDANCE ILLUMINA -> NANOPORE TRANSFER ---")
print(f"Training Samples (Illumina): {X_train.shape[0]}")
print(f"Testing Samples (Nanopore): {X_test.shape[0]}")
print(f"Transfer R2:   {r2_score(y_test, y_pred):.3f}")
print(f"Transfer RMSE: {root_mean_squared_error(y_test, y_pred):.3f} mM C")



