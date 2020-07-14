# Extended-Target-PMBM-Tracker
MATLAB implementation of the extended target PMBM tracker based on sets of trajectories

This repository contains the Matlab implementations of the Extended target Poisson multi-Bernoulli mixture (PMBM) tracker proposed in 

Yuxuan Xia, Karl Granström, Lennart Svensson, Ángel F. García-Femández, and Jason L. Williams (2019, July). Extended Target Poisson multi-Bernoulli Mixture Trackers Based on Sets of Trajectories. In 2019 22nd International Conference on Information Fusion (FUSION) IEEE.

Full text is available at https://arxiv.org/abs/1911.09025

More details on GGIW-PMBM can be found in the paper

Karl Granström, Maryam Fatemi, and Lennart Svensson. "Poisson multi-Bernoulli mixture conjugate prior for multiple extended target filtering." IEEE Transactions on Aerospace and Electronic Systems (2019).

Full text is available at https://arxiv.org/abs/1605.06311

More details on the point target PMBM tracker can be found in the paper 

Karl Granström, Lennart Svensson, Yuxuan Xia, Jason L. Williams, and Ángel F. García-Femández. "Poisson multi-Bernoulli mixture trackers: continuity through random finite sets of trajectories." In 2018 21st International Conference on Information Fusion (FUSION). IEEE, 2018.

Full text is available at https://arxiv.org/pdf/1812.05131

The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance

A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, Generalized optimal sub-pattern assignment metric, in 20th International
Conference on Information Fusion, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM


- main.m runs the extended target PMBM tracker

- data association algorithm: Stochastic Optimization. Details can be found in the paper

Granström, Karl, Lennart Svensson, Stephan Reuter, Yuxuan Xia, and Maryam Fatemi. "Likelihood-based data association for extended object tracking using sampling methods." IEEE Transactions on intelligent vehicles 3, no. 1 (2017): 30-45.

