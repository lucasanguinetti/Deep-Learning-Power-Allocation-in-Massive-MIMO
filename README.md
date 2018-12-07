# Deep-Learning-Power-Allocation-in-Massive-MIMO

This is the code package related to the follow scientific article:
Luca Sanguinetti, Alessio Zappone, Merouane Debbah 'Deep-Learning-Power-Allocation-in-Massive-MIMO' presented at ASILOMAR 2018.


# Abstract of Article

This work advocates the use of deep learning to perform max-min and max-prod power allocation in the downlink of Massive MIMO networks. More precisely, a deep neural network is trained to learn the map between the positions of user equipments (UEs) and the optimal power allocation policies, and then used to predict the power allocation profiles for a new set of UEs’ positions. The use of deep learning significantly improves the complexity-performance trade-off of power allocation, compared to traditional optimization-oriented methods. Particularly, the proposed approach does not require the computation of any statistical average, which would be instead necessary by using standard methods, and is able to guarantee near-optimal performance.

# Content of Code Package

The package contains a simulation environment, based on MATALB, that allows to produce the data samples that are needed to train the neural network. The max-min and max-prod allocation strategies are simulated with the MR and M-MMSE precoding schemes. The package contains also a simulation environment, based on KERAS, to train and execute the deep learning power allocation algorithm.

# License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
