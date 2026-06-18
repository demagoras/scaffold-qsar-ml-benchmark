from pathlib import Path
import pandas as pd

CURRENT_DIR = Path(".").resolve()

DATA_PATH = CURRENT_DIR / "data"
DATA_PATH.mkdir(parents=True, exist_ok=True)
csv_path = DATA_PATH / "lipophilicity.csv"

FIGURES_DIR = CURRENT_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

DATA_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv"

if csv_path.exists():
    print(f"File already exists at {csv_path} . Skipping download.")
else:
    print("Downloading Lipophilicity dataset from MoleculeNet...")
    print(f"[SUCCESS]: Dataset saved to {csv_path}")

df = pd.read_csv(DATA_URL)
df.to_csv(csv_path, index=False)
print(f"Dataset shape: {df.shape}")
print(f"First 5 rows of the dataset:\n{df.head()}")
print("\nInfo about the dataset:")
print(df.info())
print("\nMissing values count:")
print(df.isnull().sum())
print("\nDuplicated SMILES count:")
print(df.duplicated(subset=['smiles']).sum())

