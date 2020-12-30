# Explanation

You can run the following files in the given order to recreate the results. Start with `sampling_MIQP.m`.

`sampling_MIQP.m`: Creates 2000 samples using the function `mixed_integer_CS.m`. Saves the samples as `data_training.mat`.

`train_NN.py`: Loads the samples (`data_training.mat`) and trains a NN with ReLU activations. Converts the NN into the LCP format (Section 3.1, Lemma 2). Saves the LCP parameterization of the NN as `LCP_param.mat`.

`complementarity_system.m`: Uses the LCP paramaterization of the NN (`LCP_param.mat`) and creates the complementarity representation of the system (`sys_params.mat`) as in Section 3.2.

`piecewise_quadratic_lyap.m`: Uses the complementarity representation of the closed-loop system (`sys_params.mat`) and finds a Lyapunov function as in Section 4 (Equation 11), then saves the Lyapunov function parameters as `lyap_params.mat`.

`plot_controller.m`: Plots the neural network policy (Figure 3). (requires `LCP_param.mat`)

`plot_level_sets.m`: Plots the sublevel sets of the Lyapunov function (Figure 4) with few trajectories on top. (requires `sys_params.mat` and `lyap_params.mat`)

`run_CS.m`: Plots envelopes for 1000 trajectories with corresponding Lyapunov functions (Figure 5). (requires `sys_params.mat` and `lyap_params.mat`)
