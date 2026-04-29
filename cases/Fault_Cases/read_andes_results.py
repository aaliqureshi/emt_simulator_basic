import numpy as np

# Load the .npz file
data = np.load("ieee14_fault_barq_out.npz")

# List arrays stored inside
print(data.files)

# Access a specific array
arr0 = data["arr_0"]

print(arr0.shape)
print(arr0.dtype)
