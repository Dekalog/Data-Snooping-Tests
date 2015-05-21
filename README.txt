Example tests:

Scenario 1: We have many Trading strategies 
and we want to see whether we have acquired a trading strategy
 that has an average return higher than zero.

Data : We have a matrix 'R' of returns where each column stand for a trading strategy
and each row stands for the period (time ). So R(i,j) is Return in period 'i' by model 'j'

Note: A return of -100% is represented by a value of -1 , -10% = -0.10 , etc.


To test with 500 simulations and display 'on' , we type:

[pvalue Vlstar Vl] = WhiteRealityCheck( R ,2, 0 , 500 , 1);



LINKS 4 papers : http://www.ssc.wisc.edu/~bhansen/718/White2000.pdf
http://www.stat.purdue.edu/research/technical_reports/pdfs/1991/tr91-03.pdf



If we want to test vs a benchmark model then we simply need the Returns of the benchmark model and use them as input


[pvalue Vlstar Vl] = WhiteRealityCheck( R ,2, benchmark , 500 , 1);









Scenario 2 : We have 1000 regression models used to make predictions for 'y'
and we want to test whether there is atleast 1 model that outperforms the predictions of a benchmark model ( Random Walk for instance ).

Data : We have a matrix 'e1'; which contains the residuals of the predictions made by our 1000 models.
So ==> e1(3,8) means the residual of the third prediction with model number 8 !

We also have a column vector 'e0'; which contains the residuals of the benchmark model ( Random Walk )

Now if we want to test by Mean Squared error and 300 simulations and display 'off' , we use ==>


[pvalue Vlstar Vl] = WhiteRealityCheck( e1 ,1, e0 , 300 , 0);

Now if we want to test by Mean Absolute error and 1000 simulations and display 'on' , we use ==>

[pvalue Vlstar Vl] = WhiteRealityCheck( e1 ,3, e0 , 1000 , 0);




% This package is FREE of charge and is allowed to be distributed
% It was made because have not found a data snooping test for matlab up till date, so I felt there was a need for it 

 