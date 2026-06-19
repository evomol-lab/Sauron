import pandas as pd
df = pd.DataFrame({'grp': ['A', 'A', 'B', 'B'], 'val': ['3.500', '10.200', 'nan', '2.100']})
print(df.groupby('grp')['val'].idxmin())
