import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Graph(object):

    # t: time sequence
    # x: state vector, m*n, n is number of states
    # scale: graph size, [xmin, xmax, ymin, ymax]

    def __init__(self, t, x, scale):
        self.t = t
        self.x = x
        self.scale = scale

    def sys_evaluation(self, t, x):
        return x

    def plot_cartpole_trajectory(self):
        t = self.t
        x = self.x

        theta = x[1][:]
        x_dist = x[0][:]

        boundary = 2
        self.scale[0] -= boundary
        self.scale[1] += boundary
        self.scale[2] = -1
        self.scale[3] = 3

        L = 1.0
        h = 0.2
        w = 0.4
        pend = 0.1

        px = x_dist + L * np.sin(theta)
        py = - L * np.cos(theta)

        for i in range(len(self.t)):

            plt.plot(4*self.scale[0:2],np.array([0,0]), color = 'b', linewidth = 3)
            ax = plt.gca()
            rect1 = patches.Rectangle((x_dist[i]-w/2, -h/2), w, h, edgecolor='k', facecolor='b', linewidth=3)
            ax.add_patch(rect1)
            plt.plot([x_dist[i], x_dist[i]+ L*np.sin(theta[i])], [0, -L*np.cos(theta[i])], color = 'b', linewidth = 3)
            rect2 = patches.Rectangle((x_dist[i]+L*np.sin(theta[i])-pend/2, -L*np.cos(theta[i])-pend/2), pend, pend, edgecolor='k', facecolor='r', linewidth=3)
            ax.add_patch(rect2)
            plt.plot(px[0:i], py[0:i], color = 'g', linewidth=3)
            plt.axis('equal')
            plt.xlim(self.scale[0:2])
            plt.ylim(self.scale[2:4])
            plt.xlabel('y')
            plt.ylabel('z')
            plt.title('Cart-Pole Trajectory, t = {0:.2f}'.format(t[i]))
            plt.show()
