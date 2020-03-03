import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

x = np.linspace(-100, 100, 101)
y = x ** 2
z = x ** 2 + x * 3 + 4
a = x ** 3 + x * 3 + 4

lines = [y, z, a]
labels = ['y', 'z_', '_a']

fig = plt.figure()
for i in range(len(lines)):
    label = labels[i]
    if label.startswith('_'):
        label = label[1:]
    plt.plot(x, lines[i], label=label)

plt.legend()
plt.show()
