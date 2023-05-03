import matplotlib.pyplot as plt
import numpy as np

def plot_omega(omega, dt):
    t = np.array(range(0, len(omega)))*dt
    grad = np.gradient(omega)*100
    plt.plot(t,omega)
    plt.plot(t,grad, 'r')
    plt.show()
