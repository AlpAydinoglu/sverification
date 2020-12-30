import numpy as np
from system import LinearSystem
from vis import Graph

mc = 1
mp = 1
L = 1
g = 9.8100
A = np.array([[0, 0, 1, 0],
      [0, 0, 0, 1],
      [0, g*mp/mc, 0, 0],
      [0, g*(mc+mp)/(L*mc), 0, 0]])
B = np.array([0, 0, 1/mc, 1/(L*mc)])
B.reshape((4,1))


solution = LinearSystem(A, B).simulate_cartpole()
fig = Graph(solution.t, solution.y, np.array([-2, 10, -1, 3])).plot_cartpole_trajectory()