# Explanation

You can run the following files in the given order to recreate the results. Start with `sampling_MIQP.m`.

`example.py`: Creates 5000 samples. Saves the samples as `states_ARE.npy` and `inputs_ARE.npy`.

`train_NN.py`: Loads the samples (`states_ARE.npy` and `inputs_ARE.npy`) and trains a NN with ReLU activations using 4000 samples. Converts the NN into the LCP format (Section 3.1, Lemma 2). Saves the LCP parameterization of the NN as `LCP_param.mat`.

`complementarity_system.m`: Uses the LCP paramaterization of the NN (`LCP_param.mat`) and creates the complementarity representation of the system (`sys_params.mat`) as in Section 3.2.

`piecewise_quadratic_lyap.m`: Uses the complementarity representation of the closed-loop system (`sys_params.mat`) and finds a Lyapunov function as in Section 4 (Equation 11), then saves the Lyapunov function parameters as `lyap_params.mat`.

`run_CS.m`: Plots envelopes for 1000 trajectories with corresponding Lyapunov functions (Figure 7). (requires `sys_params.mat` and `lyap_params.mat`)
