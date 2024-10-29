# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 10:29:57 2022

@author: bjoerk_c
"""


import matplotlib.pyplot as plt
import numpy as np


class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices // 2

        self.im = ax.imshow(self.X[:, :, self.ind], cmap="gray")
        self.update()

    def onscroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


def plot3d(image):
    fig, ax = plt.subplots(1, 1)
    tracker = IndexTracker(ax, image)
    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    plt.show()


if __name__ == "__main__":
    img = np.array([[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
                     [[0, 0, 0], [1, 1, 1], [0, 0, 0]],
                     [[0, 0, 0], [0, 1, 0], [0, 0, 0]]])

    plot3d(img)
    
