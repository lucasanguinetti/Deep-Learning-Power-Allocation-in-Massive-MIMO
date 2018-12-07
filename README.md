# Deep-Learning-Power-Allocation-in-Massive-MIMO

This is the code package related to the follow scientific article:
Luca Sanguinetti, Alessio Zappone, Merouane Debbah 'Deep-Learning-Power-Allocation-in-Massive-MIMO' presented at the Asilomar Conference on Signals, Systems, and Computers, 2018. http://www.asilomarsscconf.org

 
# Abstract of Article

This work advocates the use of deep learning to perform max-min and max-prod power allocation in the downlink of Massive MIMO networks. More precisely, a deep neural network is trained to learn the map between the positions of user equipments (UEs) and the optimal power allocation policies, and then used to predict the power allocation profiles for a new set of UEs’ positions. The use of deep learning significantly improves the complexity-performance trade-off of power allocation, compared to traditional optimization-oriented methods. Particularly, the proposed approach does not require the computation of any statistical average, which would be instead necessary by using standard methods, and is able to guarantee near-optimal performance.

# Content of MATLAB Code Package

The package contains a simulation environment, based on MATLAB, that allows to produce the data samples that are needed to train the neural network. The max-min and max-prod allocation strategies are simulated with the MR and M-MMSE precoding schemes. This is the code package used to reproduce Figure 7.2 in the monograph:

Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. 

You can download a free PDF of the manuscript at https://massivemimobook.com 

The latest version of the code package for the entire book is always found in our Github repository:

https://github.com/emilbjornson/massivemimobook

# How to use the MATLAB Code Package

Set the network parameters in main.m. Then, run it. Copy the data samples in MyDataFile_{i}.mat; that is, MyDataFile_0.mat, MyDataFile_1.mat and MyDataFile_2.mat. Unfortunately, we could not provide the data files because they are too big. 

Run transform_data_maxmin.m to put toghetter and reshape all the data samples for the KERAS code package for max-min. The file transform_data_maxprod.m does the same for the max-prod package. 

MyDataFile_0.mat is used from the KERAS code package for testing the neural network after training.
 
# Content of KERAS Code Package

The package contains also a simulation environment, based on KERAS, to train and execute the deep learning power allocation algorithm with the max-min and max-prod power allocation strategies.

# License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
