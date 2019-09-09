# hBReg

Hierarchical Bayesian Regression


This is an R implementation of hierarchical Bayesian Regression.
The file provides a variational approximation of the model as well as Gibbs Sampler.


For a detailed derivation check arxiv: https://arxiv.org/abs/1811.03687


This code is part of the paper:
How to Predict Mood? Delving into Features of Smartphone-Based Data
Presented: Conference: 22nd Americas Conference on Information Systems, August 2016



Hierarchical Beta Regression

I added the Beta Regression counter part.
It utilizes Laplace approximation for the non-conjugate Beta distribution.


Decay Regression 

I added a model with decaying weights.
y= x* beta1*exp(-lambda*(0:M-1)) +beta0

