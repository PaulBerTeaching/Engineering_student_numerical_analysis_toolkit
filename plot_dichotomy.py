import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import itertools
plt.rcParams['animation.ffmpeg_path'] = '/Users/paulberraute/Downloads/'

def data_gen():
    for cnt in itertools.count():
        x = cnt * 5 / 100
        yield x, x ** 3 + 4 * x ** 2 - 10


def init():
    ax.set_ylim(-13, 50)
    ax.set_xlim(0, 5)
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    return line,

fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.grid(True, which='both')

ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')
xdata, ydata = [], []


def run(data):
    # update the data
    t, y = data
    xdata.append(t)
    ydata.append(y)
    xmin, xmax = ax.get_xlim()

    if t >= xmax:
        ax.set_xlim(xmin, 2*xmax)
        ax.figure.canvas.draw()
    line.set_data(xdata, ydata)

    return line,

ani = animation.FuncAnimation(fig, run, data_gen, interval=10, init_func=init)
writer = animation.FFMpegWriter(
     fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save("/Users/paulberraute/Downloads/movie.mp4", writer=writer)
plt.show()
