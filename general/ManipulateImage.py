from scipy import misc
from scipy import ndimage
import numpy as np
f = misc.imread("temp_figure_vikram.png")

import matplotlib.pyplot as plt
plt.close()

print(f.shape)
f_bars = f[:5700,5000:]
plt.imshow(f_bars[120:2000,:, ])
plt.show(block=False)


# # sx=ndimage.sobel(f, axis=0, mode="constant")
# # sy=ndimage.sobel(f, axis=1, mode="constant")

# # sob = np.hypot(sx, sy)

lim_y, lim_x, lim_z = f_bars.shape
bars_l = []
bars_r = []
found = False

line_2 = f_bars[100, :, ]

line_1 = f_bars[0, :, ]

line_2_pls1 = line_1[1:, ]
line_2_min1 = line_1[:-1, ]

line_2_del = line_2_pls1 - line_2_min1

temp = [np.sum(np.float_(i)) for i in f_bars[120,:, ]]
plt.close()

y = temp
x = range(len(temp))
plt.plot(x, y)
plt.show(block=False)
print(line_2_del)
# line_2_del = [np.sum(np.float_(line_2[])**2)]

# return

# for x_i in range(0, lim_x, 100):
# 	profile = [np.sum([])]
# 	for y_i in range(0, 200, 20):
# 		print(f[y_i, x_i])
