import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.linalg as linalg
import numpy as np
import time

def plot_all_var(states, n, k, m, N, name, sys_iter, length):
    if length == 1:
        print("States: ", states)
        return

    states = np.squeeze(np.asarray(states[:, :length]))
    iter = states.shape[1]

    fig = plt.figure(figsize=(6, 8))
    gs1 = gridspec.GridSpec(3, 1)
    axs = [fig.add_subplot(ss) for ss in gs1]
    for i in range(n*N):
        axs[0].plot(range(iter), states[i, :], label=r"$x_{}$".format((i % n)+1) + ", t = {}".format((i // n)+1))
        axs[0].grid(True)
        axs[0].legend(loc = 'upper right')
        if np.max(states[i, :]) - np.min(states[i, :]) > 1e2:
            axs[0].set_yscale('log')
    for i in range(n*N, (n+k)*N):
        axs[1].plot(range(iter), states[i, :], label=r"$u_{}$".format(((i-n*N) % k)+1) + ", t = {}".format((i - n*N) // k))
        axs[1].grid(True)
        axs[1].legend()
        if np.max(states[i, :]) - np.min(states[i, :]) > 1e2:
            axs[1].set_yscale('log')
    for i in range((n+k)*N, (n + k + m) * N):
        axs[2].plot(range(iter), states[i, :],
                 label=r"$\lambda_{}$".format(((i - (n+k)*N) % m) + 1) + ", t = {}".format((i - (n+k)*N) // m))
        axs[2].grid(True)
        axs[2].set_xlabel('ADMM iteration k')
        axs[2].legend()
        if np.max(states[i, :]) - np.min(states[i, :]) > 1e2:
            axs[2].set_yscale('log')

    text = "All variables in " + name + ', horizon =' + str(N) + ', sys steps = ' + str(sys_iter)
    fig.suptitle(text)
    gs1.tight_layout(fig, rect=[0, 0.05, 1, 0.95])
    plt.savefig("img/var_in_"+ name +"_step{}".format(sys_iter))
    plt.show()

def plot_all_states(states, n, N, name):
    states = np.squeeze(np.asarray(states))
    iter = states.shape[1]
    for i in range(n):
        plt.plot(range(iter), states[i, :], label=r"$x_{}$".format((i % n)+1) + ", t = {}".format((i // n)+1))
    plt.legend()
    plt.xlabel('system steps t')
    plt.title("All variables in " + name)
    plt.grid(True)
    text = 'horizon =' + str(N)
    plt.figtext(1, 0.01, text, ha='right', va='bottom')
    plt.savefig("img/var_in_"+ name)
    plt.show()

def plot_cost_vs_penalty(miqp_cost, admm_cost, true_cost, penalty_para, N):
    plt.plot(penalty_para, miqp_cost, label="miqp_cost")
    plt.plot(penalty_para, admm_cost, label="admm_cost")
    plt.plot(penalty_para, true_cost, label="true_cost")

    plt.legend()
    plt.xlabel(r"penalty parameter $\rho$")
    plt.ylabel("cost")
    plt.title("cost versus fixed penalty parameter")
    plt.grid(True)
    text = 'horizon =' + str(N)
    plt.figtext(1, 0.01, text, ha='right', va='bottom')
    plt.savefig("img/cost_vs_penalty")
    plt.show()

def plot_residual_vs_iter(s_all, r_all, rho_all, N, k_iter, sys_iter):

    s_all = s_all[:, 1:k_iter]
    r_all = r_all[:, 1:k_iter]
    rho_all = rho_all[:, 1:k_iter]

    eqs_num = s_all.shape[0]
    admm_iter = s_all.shape[1]

    fig = plt.figure(figsize=(6, 8))
    gs1 = gridspec.GridSpec(3, 1)
    axs = [fig.add_subplot(ss) for ss in gs1]

    for i in range(eqs_num):
        axs[0].plot(range(1, admm_iter+1), s_all[i, :], label=r"$||s_{}^k||_2$".format(i+1))
        axs[0].set_yscale('log')
        axs[0].grid(True)
        axs[0].legend(loc = 'upper right')

    for i in range(eqs_num):
        axs[1].plot(range(1, admm_iter+1), r_all[i, :], label=r"$||r_{}^k||_2$".format(i+1))
        axs[1].set_yscale('log')
        axs[1].grid(True)
        axs[1].legend(loc = 'upper right')

    for i in range(eqs_num):
        axs[2].plot(range(1, admm_iter+1), rho_all[i, :], label=r"$||\rho_{}^k||_2$".format(i+1))
        axs[2].set_yscale('log')
        axs[2].grid(True)
        axs[2].legend(loc = 'upper right')
        axs[2].set_xlabel('ADMM iteration k')

    text = 'horizon = ' + str(N) + ', sys steps = ' + str(sys_iter)
    fig.suptitle(text)
    gs1.tight_layout(fig, rect=[0, 0.05, 1, 0.95])
    title = "img/residual_vs_iter_step{}".format(sys_iter)
    plt.savefig(title)
    plt.show()

def plot_constraint_vs_steps(steps, cond, N, penalty_para):

    for i in range(len(penalty_para)):
        text = 'penalty parameter = ' + str(penalty_para[i])
        plt.plot(range(steps), cond[:, i], label = text)
    plt.legend()
    plt.xlabel("ADMM steps")
    plt.ylabel("equality constraint")
    plt.title("Constraint violation vs steps of ADMM")
    plt.grid(True)
    text = 'ADMM steps =' + str(steps) +', horizon =' + str(N)
    plt.figtext(1, 0.01, text, ha='right', va='bottom')
    plt.savefig("img/constraint_vs_steps_ADMM")
    plt.show()

def plot_cost_vs_steps(steps, cost, N, penalty_para):
    for i in range(len(penalty_para)):
        text = 'penalty parameter = ' + str(penalty_para[i])
        plt.plot(range(steps), cost[:, i], label = text)
    plt.legend()

    plt.xlabel("ADMM steps")
    plt.ylabel("cost")
    plt.title("Cost vs steps of ADMM")
    plt.grid(True)
    text = 'ADMM steps =' + str(steps) +', horizon =' + str(N)
    plt.figtext(1, 0.01, text, ha='right', va='bottom')
    plt.savefig("img/cost_vs_steps_ADMM")
    plt.show()

def plot_state_comp(steps, states_solver, states_sim):
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(range(steps), states_solver[0], '-k', label='ADMM states')
    axs[0, 0].plot(range(steps), states_sim[0], '--b', label='sim state')
    axs[0, 0].set_ylabel('Cart Position [m]')
    axs[0, 0].set_xlabel('time step')
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    axs[0, 1].plot(range(steps), states_solver[1], '-k', label='ADMM states')
    axs[0, 1].plot(range(steps), states_sim[1], '--b', label='sim state')
    axs[0, 1].set_title('Pole Angle')
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    axs[1, 0].plot(range(steps), states_solver[2], '-k', label='ADMM states')
    axs[1, 0].plot(range(steps), states_sim[2], '--b', label='sim state')
    axs[1, 0].set_title('Cart Velocity')
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    axs[1, 1].plot(range(steps), states_solver[3], '-k', label='ADMM states')
    axs[1, 1].plot(range(steps), states_sim[3], '--b', label='sim state')
    axs[1, 1].set_title('Pole Angular Velocity')
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    plt.show()

def plot_state(steps, states, N):
    fig = plt.figure(1)
    gs1 = gridspec.GridSpec(2, 2)
    axs = [fig.add_subplot(ss) for ss in gs1]

    axs[0].plot(range(steps), states[0], '-k')
    axs[0].plot(range(steps), np.zeros((steps,1)), '--b', lw = 2, label='steady state')
    axs[0].set_ylabel('Cart Position [m]')
    axs[0].set_xlabel('time step')
    axs[0].grid(True)
    axs[0].legend()

    axs[1].plot(range(steps), states[1], '-k')
    axs[1].plot(range(steps), np.zeros((steps, 1)), '--b', lw=2, label='steady state')
    axs[1].set_ylabel('Pole Angle [rad]')
    axs[1].set_xlabel('time step')
    axs[1].grid(True)
    axs[1].legend()

    axs[2].plot(range(steps), states[2], '-k')
    axs[2].plot(range(steps), np.zeros((steps, 1)), '--b', lw=2, label='steady state')
    axs[2].set_ylabel('Cart Velocity [m/s]')
    axs[2].set_xlabel('time step')
    axs[2].grid(True)
    axs[2].legend()

    axs[3].plot(range(steps), states[3], '-k')
    axs[3].plot(range(steps), np.zeros((steps, 1)), '--b', lw=2, label='steady state')
    axs[3].set_ylabel('Pole Angular Velocity [rad/s]')
    axs[3].set_xlabel('time step')
    axs[3].grid(True)
    axs[3].legend()


    text = 'horizon =' + str(N)
    fig.suptitle(text)
    gs1.tight_layout(fig, rect=[0, 0.05, 1, 0.95])
    fig.savefig('img/states.png')
    plt.show()

def plot_inputs_plots(steps, inputs, N):
    fig, axs = plt.subplots(3)
    axs[0].plot(range(steps), inputs[0], '-g', lw=2)
    axs[0].set_title('Cart Force [N]')
    axs[0].set_xlabel('time step')
    axs[0].grid(True)
    axs[1].plot(range(steps), inputs[1], '-b', lw=2)
    axs[1].set_title('Left Wall Contact Force [N]')
    axs[1].set_xlabel('time step')
    axs[1].grid(True)
    axs[2].plot(range(steps), inputs[2], '-b', lw=2)
    axs[2].set_title('Right Wall Contact Force [N]')
    axs[2].set_xlabel('time step')
    axs[2].grid(True)


    text = 'horizon =' + str(N)
    plt.title(text)
    fig.savefig('img/inputs.png')
    plt.show()

def plot_inputs_lines(steps, inputs, N):
    plt.plot(range(steps), inputs[0], linewidth=2, label = 'Cart Force [N]')
    plt.plot(range(steps), inputs[1], linewidth=2, label = 'Left Wall Contact Force [N]')
    plt.plot(range(steps), inputs[2], linewidth=2, label = 'Right Wall Contact Force [N]')
    plt.xlabel('time step')
    text = 'horizon =' + str(N)
    plt.title(text)
    plt.legend()
    plt.grid(True)
    plt.savefig('img/inputs.png')
    plt.show()

def plot_cost(steps, all_cost_miqp, all_cost_solver, all_cost_real, solver, N):
    plt.plot(range(steps), all_cost_miqp, '-.', lw=2, label='MIQP Cost')
    plt.plot(range(steps), all_cost_solver, '--', lw=2, label = solver + ' Cost')
    plt.plot(range(steps), all_cost_real, ':', lw=2, label ='Real Cost')
    plt.grid(True)
    plt.xlabel('system iteration t')
    plt.yscale("log")
    plt.legend()
    text = 'iter =' + str(steps) +', horizon =' + str(N)
    plt.figtext(1, 0.01, text, ha='right', va='bottom')
    plt.savefig('img/costs.png')
    plt.show()
