import pandas as pd

data = [[1.1235, 1.9654, 2.6874], [6.5124, 4.2210, 2.2899]]

df = pd.DataFrame(data)
df = df.round(1)
print(df.round(1)) 
print(df) 