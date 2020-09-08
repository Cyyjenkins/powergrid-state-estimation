# powergrid_state_estimation


The essence of the problem of estimate powergrid state is to make an unbiased estimation of the variables of equations when the number of equations is greater than the number of variables.


<br>For powergrid systems, the variable is the node voltage (i.e.,the state value of the powergrid, which consists of the real part and the imaginary part), and part of the node and branch power (containing errors due to objective conditions) are the quantity measurement, that will be given. The functional structure of the system of equations can be derived from Kirchhoff's theorem.


<br>This code provides two solutions to solve the above problem, namely the least square method and the fast decoupled method. In the least square method, the Newton iteration is adopt to iteratively calculate the node voltage. In the meantime, the fast decoupled method simplified and decoupled the functions of least square method according to the actual characteristic of high-voltage powergrid system.


<br>The input of this code is the modified standard IEEE 30-bus data, and the filename is 'iSE30Bus.txt'. We manually adjusted the input values of node and branch power to simulate the errors contained in the data under real conditions. After running over the code, the estimated powergrid data will be output to 'oStateEstimation.txt'.


<br>All the calculation steps of the two methods are encapsulated in the following file:
```
powergrid_state_estimation.m
```


To quickly understand the project, we present a simple example. Please run
```
example.m
```


Since the project contains a certain amount of computation error, we also provide earlier versions of the code in
```
ealry_version/StateEstimation.m
```
even though the coding of this script is really really really really really awful 🤦‍♂️🤦‍♂️🤦‍♂️


For any problems, please contact us at cyychenyaoyu@163.com

Copyright (c) 2020 Yaoyu Chen
