# Explanation

One should run the following files in the given order to recreate the results. Start with `sampling_MIQP.m`.

`sampling_MIQP.m`: Creates 2000 samples using the function `mixed_integer_CS.m`. Saves this data as `data_training.mat`.

`train_NN.py`: Loads the samples (`data_training.mat`) and trains a NN with ReLU activations. Converts the NN into the LCP format (Section 3.1, Lemma 2). Saves the LCP parameterization of the NN as `LCP_param.mat`.

`complementarity_system.m`: Uses the LCP paramaterization of the NN (`LCP_param.mat`) and create 

`five_carts`: Five carts with soft contacts (Example 5.4)
