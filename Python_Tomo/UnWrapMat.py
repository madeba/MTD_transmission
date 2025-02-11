# -*- coding: utf-8 -*-
"""
@author: Nicolas Verrier
"""

import matplotlib.pyplot as plt
import HoloProcessing as holo
import tifffile as tf
import numpy as np

unwrapname = "ellipse.tiff"
wrapname = "ellipse_wrap.tiff"
Ellipse = tf.imread(unwrapname)
EllipseWrap = tf.imread(wrapname)
plt.imshow(EllipseWrap, cmap="gray")
plt.show()

EllipseUnWrap = holo.unwrapping(EllipseWrap, 5.5e-6, approx=False)
plt.imshow(EllipseUnWrap, cmap="gray")
plt.show()

x = np.linspace(0,1023,1024)
plt.plot(x,EllipseWrap[512,:], label="Wrap")
plt.plot(x,EllipseUnWrap[512,:], label="Unwrap")
plt.plot(x,Ellipse[512,:], label="Groundtruth")
plt.legend()
plt.show()

Diff = Ellipse - EllipseUnWrap
plt.imshow(Diff, cmap="gray")
plt.colorbar()
plt.show()
plt.plot(x,Diff[512,:])
plt.plot()