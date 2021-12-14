import numpy as np
import matplotlib.pyplot as plt

filename = 'data.txt'
m = np.loadtxt(filename)  # читаем файл в 2D матрицу

plt.plot(m[:, 0], m[:, 1])  # строим график - ось `X` - первый столбец, `Y` - второй
plt.show()
