import numpy as np
import pandas as pd

first_array = np.array([0, 1, 2, 3])
second_array = np.array(['a', 'b', 'c', 'd'])

first_df = pd.DataFrame(index=first_array, columns=second_array)
filling = np.arange(100, 116)
#for i in range(len(first_df)):
 #   for j in range(len(first_df.columns)):
  #      first_df.at[first_df.index[i], first_df.columns[j]] = filling[len(first_df.columns) + j] 

first_df.loc[:, :] = filling.reshape(first_df.shape)

print(first_df.head(5))