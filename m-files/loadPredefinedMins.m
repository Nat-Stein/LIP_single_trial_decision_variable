% Interim solution to loading Min neurons before code is ready

load('E:\Matalb analyses\sl_allDim_minimal_230522', 'allW', 'dimName')

dim_minC = 31;
dim_minI = 32;

w_minC = allW{dim_minC}{1, sess};
w_minI = allW{dim_minI}{1, sess};
