# HTE

Identifying heterogeneous treatment effects (HTEs)—systematically different benefits or risks from a medical therapy among some participants in a study, as compared to the study average—is clinically important for the application of randomized trial data to patient care. HTE estimation equations, derived using either a conventional regression approach or newer machine learning methods, have been increasingly published in high-profile medical journals. We sought to determine how often false identification of HTEs may occur through internal development of risk equations from randomized trials.

We simulated 10,000 randomized trials of 6,000 persons each, for each of four types of trials: trials with and without positive average treatment effects, and with and without HTEs. We systematically varied the correlation among covariates, including covariates that truly affected the outcome, irrelevant correlated covariates, and covariates that determined HTE. We estimated risk/benefit for participants in the trial, and from an independent validation sample, using conventional backwards selection to derive a risk/benefit equation.  We compared the conventional approach to an ensemble of currently-applied machine learning methods (generalized linear modeling with elastic net regularization, gradient boosted trees, random forests, and deep learning).

Sanjay Basu, M.D., Ph.D.,1,2* James Faghmous, Ph.D.,3 Jeremy B. Sussman, M.D., M.S.,4,5 Brian Denton, Ph.D.,6 Joseph Rigdon, Ph.D.,7 Rodney A. Hayward, M.D.4,5

1 Center for Population Health Sciences and Center for Primary Care and Outcomes Research, Departments of Medicine and of Health Research and Policy, Stanford University
2 Center for Primary Care, Harvard Medical School
3 Arnhold Institute for Global Health, Icahn School of Medicine, Mt. Sinai Medical Center 
4 Center for Clinical Management Research, Veterans Affairs Ann Arbor Healthcare System
5 Division of General Medicine, University of Michigan
6 Department of Industrial and Operations Engineering, University of Michigan
7 Quantitative Sciences Unit, Stanford University

*to whom correspondence should be addressed:
Email: basus@stanford.edu
