from rinpy import RINProcess
import pandas as pd
import os

proc = RINProcess(output_path="rinpy_out", input_path="rinpy_in", num_workers=1)
proc.start_process()

for root, dirs, files in os.walk("rinpy_out"):
    for file in files:
        if file.endswith('.csv') or file.endswith('.tsv') or file.endswith('.txt'):
            print(f"--- {file} ---")
            try:
                df = pd.read_csv(os.path.join(root, file))
                print(df.head(2))
            except Exception as e:
                with open(os.path.join(root, file)) as f:
                    print(f.read()[:200])
