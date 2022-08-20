import array as ar
from enum import EnumMeta

samenvatting = [""]*8

samenvatting[1] = (
'''
Variable = general propery of an object, allows to distinguish objects
Value = specific propery, interpretation for that variable

Measurement Levels
------------------
Qualtitative = Limited number of values, not necessarily numeric (bv: small, medium, large)
    Nominal = categories, e.g. gender, race, country, shape,...
    Ordinal = Order, rank e.g. military rank, level of education
Quantitative = Many values, often unique. Number + unit of measurement (often contain result of measurement)
    Interval = No fixed zero point => no proportions, e.g. °C, °F
    Ratio = Absolute zero point => proportions, e.g. distance (m), energy (J), weight (kg), temperature in Kelvin,...

Causal Relationships
--------------------
Cause = Independent variable
Consequence = Dependent variable
***WARNING!!! A relationship between variables does not necessarily indicate a CAUSAL relation!***

Sample and Population
---------------------
Population = the collection of all objects/people/... that you want to investigate
Sample = a SUBSET of the POPULATION from which measurements will be taken

Under certain circumstances, the results for a sample are representative for the population

Sampling method
---------------
1. Definition of population
2. Define sampling frame
3. Choice of sampling method (budget and time)

Random sample = every element from the population has an equal chance of being included in the sample
Non-random sample = convencience sampling = the elements for the sample are NOT randomly selected. Objects that can be collected EASILY are more likely to be included

Possible errors
---------------
Errors = measurements in a sample that deviate from the value in the entire population

Accidental      <=>     Systematic
Sampling error  <=>     Non-Sampling error

    Sampling errors
    ---------------
    1. Accidental sampling errors = pure coincidence
    2. Systematic sampling errors:
        a. Online Survey: people without internet are excluded
        b. Street Survey: only who is currently walking there
        c. Voluntary Survey: only interested parties participate
    
    Non-Sampling errors
    -------------------
    1. Accidental non-sampling errors = incorrectly ticked answers
    2. Systematic non-sampling errors:
        a. Poor or non-calibrated measuring instruments
        b. Value can be influenced by the fact that you measure
        c. Respondents lie
''')

samenvatting[2] = (
'''
Measure of Central Tendency: (what value is representative of the entire group?)
****************************
Mean or average (pdf p7)
------------------------
Arithmitic mean = the arithmetic mean is the sum of all values divided by the number of values

Mediaan (pdf p9)
----------------
= sort all values and pick the middle number
IF even number of values: average of the middle two

Mode (pdf p11)
--------------
= the value that appears most often in a dataset

Measure of Dispersion: (How large are the differences within the group?)
**********************

Range (pdf p13)
----------------
= absolute value of the difference between the highest and the lowest value

Quartiles (pdf p14)
-------------------
The quartiles of a sorted set of numbers are the three values that divide the set into 4 equally large subsets. Q1, Q2 and Q3

calculating:
    1. Sort the values:
        *When n is odd:
         Median (Q2) is the middle value
         Leave out the median, Q1 is the median of the first half, Q3 is the median of the second half
        
        *When n is even
         Median (Q2) is the average of the two middle values
         Q1 is the median of the first half, Q3 is the median of the second half
         two middle values that form Q2 are included in the dataset for calculating Q1 and Q3

Variance and standard deviation (pdf p16)
*****************************************
Variance
--------
= the mean squared difference between the values of a data set and the arithmetic mean

Standard deviation
------------------
= the square root of the variance

***NOTE***
Providing center value (e.g. mean, median,...) is never sufficient. (What is the dispersion?)

*********************************************************************************
*************************************SUMMARY*************************************
*********************************************************************************
|-----------------------|------------------|------------------------------------|
|   Measurement Level   |   Center         |   Spread Distribution              |
|-----------------------|------------------|------------------------------------|
|   Qualtitative        |   Mode           |   /                                |
|-----------------------|------------------|------------------------------------|
|   Quantitative        |   Average/Mean   |    Variance, Standard Deviation    |
|                       |   Median         |    Range, Interquartile Range      |
|-----------------------|------------------|------------------------------------|

Summary of Symbols (population vs sample) => pdf p26

Data Visualisation
*******************

Chart type overview:
--------------------
|-----------------------|------------------|
|   Measurement Level   |   Chart Type     |
|-----------------------|------------------|
|   Qualtitative        |   Bar chart      |
|-----------------------|------------------|
|   Quantitative        |   Boxplot        |
|                       |   Histogram      |
|                       |   Density plot   |
|-----------------------|------------------|

Avoid using a Pie chart; comparing angles is harder than comparing length

Interpretation of charts, tips:
1. Label the axes
2. Add a clear title
3. Name the unit
4. Add a label that clarifies the chart

Anscombe's Quartet (pdf p35):
= four completely different datasets with the same measurements of central tendency and dispersion

'''
)

samenvatting[3] = (

'''
Probability Distribution of a Sample
************************************

Probabilities represent beliefs of how likely it is that a certain event wil happen when performing a certain experiment

Observations:
    *Probabilties are numbers assigned to sets
    *These sets are part of some ALL-ENCOMPASSING set = the UNIVERSE (denoted greek capital omega)

Axioms of Probability
---------------------

1.  Probabilities are non-negative: P(A)>=0 for each A 
2.  The universe has probability 1: P(OMEGA) = 1
3.  When A and B are disjoint events (i.e. A intersect B = 0)
    => P(A union B) = P(A) + P(B)
    => SUM RULE

Properties of Probabilities
---------------------------
= derived from the three axioms of probability

1. COMPLEMENT RULE: for each A it holds that
        P(/A) = 1 - P(A)
2. The impossible event has probability zero: P(0) = 0
3. GENERAL SUM RULE: P(A union B) = P(A) + P(B) - P(A intersect B)

Independent events = the occurence of one event (or knowing that the event occurred) does not change the probability of the other event occurring
=> A and B are independent when: P(A intersect B) = P(A)P(B)

Random Variable => Probability Distribution Function

Expectation of a random variable (pdf p37)

Variance of a random variable (pdf p39)

Continuous Random Variable (pdf p41)
************************************

=> takes on an uncountably infinite number of possible values
so the probability that X equals to some A exactly is not useful, because this is always zero
=> probability that X takes on a value in some interval [a,b]
=> this can be found by integrating the probability density function of the random variable

CONTINUOUS PROBABILITY DISTRIBUTION (pdf p43)
*   Normal Distribution
*   Exponential Distribution
*   continuous Uniform Distribution

Standard Normal distribution (pdf p47)

Python functions
----------------

import scipy.stats

For normal distribution with mean m and standard deviation s: (TIP: draw Gauss curve to visual this)
|----------------------------------|--------------------------------------|
|   Function                       |   Purpose                            |
|----------------------------------|--------------------------------------|
|   norm.pdf(x, loc=m, scale=s)    |   Probability density at x           |
|----------------------------------|--------------------------------------|
|   norm.cdf(x, loc=m, scale=s)    |   Left-tail probability P(X < x)     |
|----------------------------------|--------------------------------------|
|   norm.sf(x, loc=m, scale=s)     |   Right-tail probability P(X < x)    |
|----------------------------------|--------------------------------------|
|   norm.isf(1-p, loc=m, scale=s)  |   p% of observations are expected    |
|                                  |   to be lower than result            |
|----------------------------------|--------------------------------------|

Exponentional Distribution (pdf p58): values for an exponential random variable occur when there are fewer large values and more small values

Continuous Uniform Distribution (pdf p60): the density function is constant where every value has an equal chance of occurring

The Central Limit Theorem (pdf p63)
***********************************
=> if the size of the sample is sufficiently large, the probability distribution of the sample mean will approximate a normal distribution, regardless of the probability distribution of the underlying population

Consider:
    *Population with expected value µ and standard deviation sigma
    *random sample of n observations

=> probability distribution of the sample mean (_x) will approximate a normal distribution with:
    *mean µ_x = µ
    *standard deviation sigma_x = sigma/sqrt(n)

NOTE: the larger the sample, the better the probability distribution of _x will approximate the expected value of the population, µ

Point Estimate
--------------
A point estimate for a population parameter is a formula or equation that allows us to calculate a value to estimate that parameter

Confidence Interval
-------------------
A confidence interval is an equation or formula that allows us to construct an interval that will contain the parameter to be estimated with a certain level of confidence

    ***Large Sample*** (pdf p67)
    ------------------
    Give a sample with mean _x
    => looking for interval = [_x - b, _x + b] for which we can say with a level of confidence (1 - alpha) of e.g. 95% that µ is inside this interval

    P(_x - b < µ < _x + b) = 1 - alpha = 0.95

    calculate the z-score for _x: z= (_x - µ) /(sigma / sqrt(n))

    calculate z_alpha/2 for which: P(-z_alpha/2 < z < z_alpha/2) = 1 - alpha = 0.95 (=> alpha = 1-0.95 = 0.05 => alpha/2 = 0.05/2 = 0.025)
    => P(z < z_aplha/2) = 1 - alpha/2 = 0.975

    z_alpha/2 = stats.norm.isf(1-0.975) = +-1.96

    ***Small Sample*** (pdf p71)
    ------------------
    != central limit theorem no longer valid

    If population X has a normal distribution (X ~ Nor(µ, sigma)) and you take a small sample with mean _x and standard deviation s, then:

    t = (_x - µ) / (sigma / sqrt(n))

    t_alpha/2 = t variant of the z_alpha/2, so same calculation method, but df must be added (n-1) = stats.t.isf(1-0.975, df=4)

    will behave as a t-distribution with n - 1 degrees of freedom (df)

    import scipy.stats

    For a t-distribution with df degrees: 
    |------------------------|--------------------------------------|
    |   Function             |   Purpose                            |
    |------------------------|--------------------------------------|
    |   t.pdf(x, df=d)       |   Probability density at x           |
    |------------------------|--------------------------------------|
    |   t.cdf(x, df=d)       |   Left-tail probability P(X < x)     |
    |------------------------|--------------------------------------|
    |   t.sf(x, df=d)        |   Right-tail probability P(X < x)    |
    |------------------------|--------------------------------------|
    |   t.isf(1-p, df=d)     |   p% of observations are expected    |
    |                        |   to be lower than result            |
    |------------------------|--------------------------------------|

Hupothesis testing (2e pdf van H3):
***********************************

Statistical Hypothesis Testing (pdf2 p5)
------------------------------

Hypothesis = idea that has yet to be proven: statement regarding numeric value of a population parameter
Hypothesis Test = verification of a statement about the values of one or multiple population parameters
Null Hypothesis (H0) = base hypothesis, we start with assuming it is true
Alternative Hypothesis (H1, Ha) = conclusion if the null hypothesis is unlikely to be true

Elements of a testing procedure (pdf2 p6)
-------------------------------

Test Statistic = the value that is calculated from the sample
Region of Acceptance = the region of values SUPPORTING the null hypothesis
Critical Region / Region of Rejection = the region of values REJECTING the null hypothesis
Significance Level = the probability of rejecting a true null hypothesis H0

Testing procedure (pdf2 p7)
-----------------
1.  Formulate both hypotheses (H0 and H1)
2.  Determine the significance level (alpha)
3.  Calculate the test statistic
4.  Determine the critical region or the probability value
5.  Draw conclusions

Probability Value (pdf2 p13)
-----------------
p-value is the probability, if the null hypothesis is true, to obtain a value for the test statistic that is at least as extreme as the observed value

p-value < alpha => reject H0; the discovered value of _x is too extreme
p-value >= alpha => do not reject H0; the discovered value of _x can still be explained by coincidence

Cricital Region (pdf2 p16)
---------------
= the collection of all values of the test statistic for which we can reject the null hypothesis

Critical value g

|------------------------|----------------------------------------------------------|
|   Goal                 |   Test regarding the value of the population mean µ      |
|                        |   using a sample of n independent values                 |
|------------------------|----------------------------------------------------------|
|   Prerequisite         |   The population has a random distribution, n is         |
|                        |   sufficiently large                                     |
|------------------------|----------------------------------------------------------|
|   Test Type            |   Two-tailed     |   Left-tailed     |   Right-tailed    |
|------------------------|----------------------------------------------------------|
|   H0                   |   µ = µ0             µ = µ0              µ = µ0          |
|   H1                   |   µ != µ0            µ < µ0              µ > µ0          |
|   Critical Region      |   |_x| > g           _x < -g             _x > g          |
|   Test statistic       |   z = (_x - µ0) / (sigma / sqrt(n))                      |
|------------------------|----------------------------------------------------------|
Tabel: Summary of Testing Procedures (pdf2 p22)

Requirements for z-test (pdf2 p23)
-----------------------
    *   The sample needs to be random
    *   The sample size needs to be sufficiently large (n>=30)
    *   The test statistic needs to have a normal distribution
    *   The standard deviation of the population, sigma, is known

Student's t-test (pdf2 p35)
----------------
If requirements for a z-test are not met, e.g.:
* Sample size too small
* Population stdev (sigma) unknown
IF the variable is normally distributed => t-test

T-test (pdf2 p36)
------
Determine critical value:
    g = µ +- t x s / sqrt(n)

    * t-value is derived from the Student's t-distribution, based on the number of degrees of freedom, n-1
    * Use function t.isf() in Python
    * Otherwhise same procedure as the z-test

Errors in Hypothesis Tests (pdf2 p38)
--------------------------

|------------------------|---------------------------------------------------------------------|
|   Conclusion           |                          Reality                                    |
|                        |            H0 True                |           H1 True               |
|------------------------|---------------------------------------------------------------------|
|   H0  not rejected     |       Correct inference           |  Type II error (false negative) |
|   H0  rejected         |  Type I error (false positive)    |          Correct inference      |
|------------------------|---------------------------------------------------------------------|

P(type I error) = alpha (=significance level)
P(type II error) = beta
Calculating beta is NOT trivial, but if alpha declining then beta rising
'''

)

def help():
    commands = {'searchString(str)','overzichtHs()','printTheorie()','printH1()'}

    for i in commands:
        print(i)

def overzichtHs():
    print(
'''
H1 = Basisbegrippen, steekproefonderzoek
H2 = Analyse van 1 variabele
H3 = De centrale limietstelling, betrouwbaarheidsintervallen, hypothesetoetsen
H4 = Analyse van 2 kwalitatieve variabelen
H5 = Analyse van kwalitatieve vs kwantitatieve variabele
H6 = Analyse van 2 kwantitatieve variabelen
H7 = Tijdserie-analyse
'''
    )

def searchString(str):
    for index, s in enumerate(samenvatting):
        if str in s:
            print(f"H{index}")

def printTheorie():
    for index, s in enumerate(samenvatting):
        print(f'Hoofdstuk {index}:\n{s}')

#-------H1------

def printH1():
    
    print(samenvatting[1])

def printH2():
    
    print(samenvatting[2])
