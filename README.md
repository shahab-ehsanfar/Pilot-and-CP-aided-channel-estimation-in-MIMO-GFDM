# Pilot-and-CP-aided-channel-estimation-in-MIMO-GFDM
Source codes for the paper: 

Ehsanfar, S., Matthé, M., Chafii, M., & Fettweis, G. P. (2018). Pilot- and CP-Aided Channel Estimation in MIMO Non-Orthogonal Multi-Carriers. IEEE Transactions on Wireless Communications, 18(1), 650-664.

Important:

Some functions in the source code have been removed in later releases of MATLAB. You need MATLAB 2018a or older (between 2015a and 2018a) to be able to run the code. 

The modulation and demodulation process in the file main.m are based on the C/C++ based coded modulation library (CML) provided by iterative solutions (www.iterativesolutions.com). You need a Matlab C/C++ compiler to run the BER simulations. For further information on how to install the CML library please refer to http://www.iterativesolutions.com/user/image/readme.pdf

To test if the CML library has been properly installed, type:
edit LTEScenarios

run the file, and then type: 
[sim_param, sim_state] = CmlSimulate('LTEScenarios', [1 3 5 11]);

If the code runs without Matlab crashing, you should be able to run main.m.





