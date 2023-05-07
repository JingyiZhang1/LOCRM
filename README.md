# Local Continumal Reassessment Methods

R codes to implement local continual reassessment methods for dose finding and optimization in drug-combination trials.

# Description

Due to the limited sample size and large dose exploration space, obtaining a desirable dose combination is a challenging task in the early development of combination treatments for cancer patients. Most existing designs for optimizing the dose combination are model-based, requiring significant efforts to elicit parameters or prior distributions. Model-based designs also rely on intensive model calibration and may yield unstable performance in the case of model misspecification or sparse data. We propose to employ local, underparameterized models for dose exploration to reduce the hurdle of model calibration and enhance the design robustness. Building upon the framework of the partial ordering continual reassessment method (POCRM), we develop local data-based CRM (LOCRM) designs for identifying the maximum tolerated dose combination (MTDC), using toxicity only, and the optimal biological dose combination (OBDC), using both toxicity and efficacy, respectively. The LOCRM designs only model the local data from neighboring dose combinations. Therefore, they are flexible in estimating the local space and circumventing unstable characterization of the entire dose-exploration surface. Our simulation studies show that our approach has competitive performance compared to widely used methods for finding MTDC, and it has advantages over existing model-based methods for optimizing OBDC.

# Functions

The repository includes two functions:

- locrm.R: The R code that includes the function ```locrm()``` to obtain the operating characteristics of the LOCRM design for single agent trials with late-onset toxicities by simulating trials.
  
  ```rscript
  locrm(p.true.tox,target.tox,cutoff.tox,ntrial,nmax,cohortsize)
  ```
  
- locrm12.R: The R code that includes the function ```locrm12()``` to obtain the operating characteristics of the time-to-event model-assisted design for single agent trials with late-onset toxicities by simulating trials.
  
  ```rscipt
  locrm12(p.true.tox,p.true.eff,target.tox,cutoff.tox,target.eff,cutoff.eff,ntrial,nmax,cohortsize)
  ```
  

# Inputs

- `p.true.tox`: A matrix containing the true toxicity probabilities of the investigational dose combinations.
  
- `p.true.eff`: A matrix containing the true efficacy probabilities of the investigational dose combinations.
  
- `target.tox`: The upper limit of toxicity rate, e.g., `target.tox<-0.3`.
  
- `target.eff`: The lower limit of efficacy rate, e.g., `target.eff<-0.2`.
  
- `cutoff.tox`: The cutoff to eliminate an overly toxic dose for safety.
  
- `cutoff.eff`: The cutoff to eliminate an ineffective dose for futility.
  
- `ntrial`: The total number of trial to be simulated.
  
- `nmax`: The maximum sample size.
  
- `cohortsize`: The cohort size.
  

# Outputs

- `locrm()` will return the operating characteristics of the LOCRM design as a data frame, including:

  ```
   (1) selection percentage at each dose combination (sel);  
   (2) the number of toxicities treated at each dose combination (tox);  
   (3) the number of patients treated at each dose combination (pts);  
   (4) the average number of toxicities (ntox);  
   (5) the average number of patients (npts);  
   (6) the percentage of early stopping without selecting the MTDC (p.stop)  
   (7) the average DLT rate (dlt.rate)
  ```
  
- `locrm12()` will return the operating characteristics of the LOCRM12 design as a data frame, including:
  
  ```
   (1) selection percentage at each dose combination (sel);  
   (2) the number of toxicities treated at each dose combination (tox);  
   (3) the number of efficacies treated at each dose combination (eff);  
   (4) the number of patients treated at each dose combination (pts);  
   (5) the average number of toxicities (ntox);  
   (6) the average number of efficacies (neff);  
   (7) the average number of patients (npts);  
   (8) the percentage of early stopping without selecting the OBDC (p.stop)
  ```

# Example

We consider use the LOCRM design as an illustration.

- Suppose the upper limit of toxicity rate is 35%, the lower limit of efficacy rate is 20%. A maximum of 51 patients will be recruited in cohorts of 3. Suppose 5 dose levels of agent A and 3 dose levels of agent B are considered are considered. We can use the following code to simulate the scenario 1 in Table 3.
  
```rscript
  p.true.tox <- matrix(c(0.05,0.15,0.30,0.45,0.55,  0.15,0.30,0.45,0.55,0.65,  0.30,0.45,0.55,0.65,0.75),nrow = 3, byrow = TRUE)
  p.true.eff <- matrix(c(0.05,0.25,0.50,0.55,0.60,  0.25,0.50,0.55,0.60,0.65,  0.50,0.55,0.60,0.65,0.70),nrow = 3, byrow = TRUE)
  locrm12(p.true.tox,p.true.eff,target.tox=0.35,cutoff.tox=0.85,target.eff=0.2,cutoff.eff=0.9,nmax=51,cohortsize=3,ntrial=100)
  
  -----------------------output------------------------

$sel
     [,1] [,2] [,3] [,4] [,5]
[1,]    0    7   16    3    0
[2,]    5   30    9    0    0
[3,]   21    9    0    0    0

$tox
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.17 0.50 2.37 1.45 0.01
[2,] 0.68 3.03 2.06 0.29 0.02
[3,] 2.06 2.61 0.60 0.02 0.00

$eff
     [,1] [,2] [,3] [,4] [,5]
[1,] 0.14 0.97 3.64 1.89 0.04
[2,] 1.21 5.13 2.37 0.29 0.01
[3,] 3.86 3.05 0.55 0.01 0.00

$pts
     [,1] [,2] [,3] [,4] [,5]
[1,] 3.60 3.63 7.35 3.24 0.09
[2,] 4.44 9.93 4.23 0.48 0.03
[3,] 7.35 5.67 0.93 0.03 0.00

$ntox
[1] 15

$neff
[1] 25

$npts
[1] 51

$p.stop
[1] 0
```

  # Authors and Reference
  * Jingyi Zhang, Fangrong Yan, Nolan A. Wages, and Ruitao Lin
  * Zhang, J., Yan, F., Wages, N. A., and Lin, R. (2023) “Local Continual Reassessment Methods for Dose Finding and Optimization in Drug-Combination Trials”, Statistical Methods in Medical Research, revision submitted to journal.

