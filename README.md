# BDStest_JER
This code is used for replicating the illustration section of Kaido &amp; Zhang (2023) JER paper.

Project Title: “Applications of Choquet Expected Utility to Hypothesis Testing with Incompleteness”
Project Description: The code is used to replicate the numerical illustration part of the paper (Figure 1 and Table 1).
Project Authors: Hiroaki Kaido and Yi Zhang (2023)

If you have any questions about code implementation, please send emails to yzhangjnu@outlook.com

Project replication files:
Before starting to replicate the results of the paper, you need to first install an optimization solver called CVX.(For new solver user, please check Reference guide — CVX Users' Guide (cvxr.com)
For installation of Solver, please check Installation—CVX Users' Guide (cvxr.com) for more details.
After finishing installation, please download the zip file and uncompress all files in the same folder.


There are two main files you need to use. The rest of the files are function files for the algorithm we proposed in the paper.
1.	File name: main_tab1fig1rep.m
This file is used to replicate the results of figure 1 and table 1. Follow the file description inside this main file to get all results accordingly.
2.	File name: main_plotQstar.m
This file is used to plot the maximum risk as a function of zeta value and provide root-finding guidance in page 14 of the paper.
