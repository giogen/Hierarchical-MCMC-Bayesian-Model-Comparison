# Hierarchical-MCMC-Bayesian-Model-Comparison

This repository contains R and Matlab (Python version coming soon) templates to perform Bayesian model comparison using a multidimensional MCMC approach with hierarchical model probabilites.

Coming soon: a version of the templates that relies on pseudoprios, as described for example in "Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and Stan" by John Kruschke.

The templates are based on nested hierarchical logistic regression models, i.e., they compare models with increasing complexity at both the population level and at the level of individual participants or data sources. On top of a standrad implementation of the nested logistic regression models, the templates contain a "model index" parameter that is estimated alongisde the model coefficients with the goal of quantifying the posterior evidence in favor of each model. Priors for the hiearchical model probabilities are defined using a Dirichlet distribution that assigns equal prior probabilites to the different models under investigations.

For comprehensive details on the experiment the templates are based on, or to obtain access to the full dataset itself, please contact the author at gentile.giovanni@gmail.com or ggentile@caltech.edu.
