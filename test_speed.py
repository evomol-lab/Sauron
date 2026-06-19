import pandas as pd
import time
import numpy as np

# Create 100k edges
data = {
    'NodeId1': np.random.randint(0, 1000, 100000).astype(str),
    'NodeId2': np.random.randint(0, 1000, 100000).astype(str),
    'Interaction': np.random.choice(['HBOND:MC_MC', 'VDW:SC_SC', 'IONIC:SC_SC'], 100000),
    'Distance': np.random.uniform(2.0, 5.0, 100000).astype(str)
}
df = pd.DataFrame(data)

start = time.time()
df['Int_Prefix'] = df['Interaction'].apply(lambda x: x.split(':')[0] if ':' in x else x)
df = df.loc[df.groupby(['NodeId1', 'NodeId2', 'Int_Prefix'])['Distance'].idxmin()]
df = df.drop(columns=['Int_Prefix'])
print("Time taken:", time.time() - start)
