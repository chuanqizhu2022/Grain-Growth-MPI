import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nx = 100
ny = 100
ns = 1000
# step_arr = np.arange(1000, ns*1, ns)
# for step in step_arr:
step = 1000
df = pd.read_csv(f"data/2d{step}.csv", header=None)
arr = df[0].values
mat = arr.reshape(nx, ny)
plt.imshow(mat)
plt.savefig(f"figs/2d{step}.png")
plt.close()
