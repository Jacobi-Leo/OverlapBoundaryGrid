import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

point_set = []
f = open("raw.dat")
data = f.readlines()
f.close()
f2 = open("grid.dat")
grid = f2.readlines()[0].split()
grid = map(float, grid)
f2.close()
points = []
for item in data:
    points.append(map(float,item.split()))
del data
#plt.plot(grid, points[200], '-o')
#plt.show()
points2 = []
f2 = open("raw2.dat")
data = f2.readlines()
for item in data:
    points2.append(map(float, item.split()))
del data
point_set.append(points)
point_set.append(points2)

# if len(points[0]) > len(grid):
#     grid.append(1.)
#     grid.insert(0, 0.0)

#fig = plt.figure()
#window = fig.add_subplot(111)
#line, = window.plot(grid, points[0])
#line2, = window.plot(grid, points2[0])

#def update(data):
#    line.set_ydata(data)
#    return line,

#ani = animation.FuncAnimation(fig, update, points, interval=2*10)
#plt.show()

plt.plot(grid, points2[50])

# fig = plt.figure()

# ax = plt.axes(xlim=(0, 1), ylim=(0, 1))

# N = 2
# lines = [plt.plot([], [], 'o-')[0] for _ in range(N)]

# def init():
#     for line in lines:
#         line.set_data([], [])
#     return lines

# def animate(i):
#     for j,line in enumerate(lines):
#         line.set_xdata(grid)
#         line.set_ydata(point_set[j][i])
#     return lines

# # Set up format files for movie
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=20, metadata=dict(artist='Z. Y. Liu'), bitrate=1800)

# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                interval=20, blit=True)

plt.show()

# anim.save("line.mp4", writer = writer)
