# TPPE29 - Financial Markets and Instruments

## Option Pricing using the Binomial Model

This repository contains Matlab code for pricing European, American, and exotic options using a binomial tree framework. The project was completed as part of the TPPE29 course in Financial Markets and Instruments at Linköping University.

## Folder Structure
```plaintext
assignments/ 
├── masterfile.m % Main script to run all tasks 
├── task1.m % European call option without dividends 
├── task2.m % European call on OMXS30 incl. implied volatility estimation 
├── task3.m % European & American calls with dividends 
├── task4.m % American call on SAAB B with dividend 
├── task5.m % Up-and-in barrier option pricing
```

## Highlights

- Custom implementation of binomial tree option pricing
- Implied volatility estimation using bisection method
- Analysis of convergence vs. Black-Scholes price
- Support for dividend-adjusted pricing
- Plotting of option value convergence

## Requirements

- Matlab with support for `days252bus()` (for trading day calculations)

Run `masterfile.m` to execute all tasks.