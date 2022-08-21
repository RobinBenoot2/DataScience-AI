import array as ar
from enum import EnumMeta

samenvatting = [""]*9

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

samenvatting[4] = (
'''
Bivariate analysis:
*******************
= determning whether there is an association between two stochastic variables (X and Y)

Association = you can predict (to some extent) the value of Y from the value of X
    * X = independent variable
    * Y = Dependent variable
***IMPORTANT! Finding an association does NOT imply a causal relation!***

Overview (pdf p6):
|-----------------------|------------------|------------------------------------|
|   Independent         |   Dependent      |   Test/Metric                      |
|-----------------------|------------------|------------------------------------|
|   Qualtitative        |   Qualtitative   |   X²-test                          |
|                       |                  |   Cramér's V                       |
|-----------------------|------------------|------------------------------------|
|   Qualtitative        |   Quantitative   |   two-sample t-test                |
|                       |                  |   Cohen's d                        |
|-----------------------|------------------|------------------------------------|
|   Quantitative        |   Quantitative   |    -                               |
|                       |                  |    Regression, correlation         |
|-----------------------|------------------|------------------------------------|

QUALITATIVE - QUALITATIVE variate:
----------------------------------

Data Visualization (pdf p10):
------------------
* Clustered bar chart
* A mosaic plot

Contingency tables = table in matrix format that displays the (multivariate) frequency distribution of the variables
Marginal totals = the totals of each value of the variables

Expected values = if there is no difference (so this imply an association), we expect the same ration in each column of the contingency table

measuring dispersion = how far is the observed value o from the expected e:
=>  (o - e)²
    --------
       e

The chi-squared statistic
-------------------------

The sum of all thes values is notated:

    X² = sum((o_i - e_i)² / e_i)

    * X = Greek letter chi
    * o_i = number of observations in the i'th cell of the contingency table
    * e_i = expected frequency
    * Small value => NO association
    * Large value => association

Chi-squared value is dependent of the table size
=> Cramér's V is independent of the table size

    V = sqrt(X² / (n(k - 1)))

    with:
        * n = the number of observations
        * k = min(numRow, numCols)

    |------------------------|--------------------------------------|
    |   Cramér's V           |   Interpretation                     |
    |------------------------|--------------------------------------|
    |   +-0                  |   no association                     |
    |   +-0.1                |   weak association                   |
    |   +-0.25               |   moderate association               |
    |   +-0.5                |   strong association                 |
    |   +-0.75               |   very strong association            |
    |   +-1                  |   complete association               |
    |------------------------|--------------------------------------|

Chi-squared test for independence
----------------------------------
 = Alternative to Cramér's V to investigate association between qualtitative variables

 Chi-squared distribution in Python
-----------------------------------
import scipy.stats

For a X²-distribution with df degrees of freedom:

    |------------------------|--------------------------------------|
    |   Function             |   Purpose                            |
    |------------------------|--------------------------------------|
    |   chi2.pdf(x, df = d)  |   Probability density for x          |
    |   chi2.cdf(x, df = d)  |   Left-tail probability P(X < x)     |
    |   chi2.sf(x, df = d)   |   Right-tail probability P(X > x)    |
    |   chi2.isf(x, df = d)  |   p% of observations is expected     |
    |                        |   to be lower than this value        |
    |------------------------|--------------------------------------|

    Test procedure
    --------------
    1. Formulate hypotheses:
        * H0: there is NO association (X² is "small")
        * H1: there is an association (X² is "large")
    2. Choose significance level, e.g. alpha = 0.05
    3. Calculate the test statistic, X²
    4. Use df = (numRow - 1) x (numCol - 1) and either:
        * Determine critical value g so P(X² > g) = alpha
        * Calculate the p-value
    5. Draw conclusion:
        * X² < g: do not reject H0  <=>     X² > g: reject H0
        * p > alpha: do not reject  <=>     p < alpha: reject H0

Test for independence in Python
-------------------------------
Function to calculate X² and p-value from a contingency table:

observed = pd.crosstab(index, columns)
chi2, p, df, expected = stats.chi2_contingency(observed)

Goodness-of-fit test
--------------------
The X² test can also be used to determine whether a sample is representative for the population

    X² = sum((o_i - e_i)² / e_i)

Value of X²:
    * Small     =>  distribution is representative
    * Large     =>  distribution is NOT representative
X² measures the degree of conflict with the null hypothesis

* The test statistic X² follows the X² distribution
* Critical value g from the X² distribution: this is dependent on the number of degrees of freedom (df). In general df = k -1; with k the number of categrories
* The critical value g for a given significance level alpha and number of degrees of freedom df can be calculated in Python using the function isf()
    P(X² < g) = 1 - alpha

    Test procedure
    --------------
    1. Formulate hypotheses
        * H0: sample is representative for the population
        * H1: sample is not representative for the population
    2. Choose significance level, e.g. alpha = 0.05
    3. Calculate the test statistic, X²
        * Critical area: calculate g so that P(X² < g) = 1 - alpha
        * Probabiltiy value: calculate p = 1 - P(X < X²)
    4. Conclusion (the test is always right-tailed):
        * X² < g: do not reject H0  <=>     X² > g: reject H0
        * p > alpha: do not reject  <=>     p < alpha: reject H0

Goodness-of-fit test in Python (pdf p36)
------------------------------

observed = np.array(*[OBSERVED VALUES]*)
expected_p = np.array(*[EXPECTED PERCENTAGES]*)
expected = expected_p * sum(observed)
chi2, p = stats.chisquare(f_obs = observed, f_exp = expected)

Standardized Residuals (pdf p44)
----------------------
indicate which classes make the greatest contribution to the value of X²

              o_i - n * π_i
    r_i = -----------------------
          sqrt(n * π_i (1 - π_i))

* r_i ∈ [-2,2]  => "acceptable" values
* r_i < -2      => underrepresented
* r_i > 2       => overrepresented

Cochran's rules (pdf p45)
---------------
In order to apply the X²-test, the following conditions must be met (Rule of Cochran)
    1. For all categories, the expected frequency e must be greater than 1
    2. In a maximum of 20% of the categories, the expected frequency e may be less than 5
'''
)

samenvatting[5] = (
'''
Bivariate analysis:
*******************
QUALITATIVE - QUALITATIVE variate:
----------------------------------

Data Visualization (pdf p6):
------------------
* Chart types for quantitative data
* Grouped by qualitative variable

Suitable chart types:
* Grouped boxplot (pdf p9)
* Grouped violin plot (pdf p10)
* Grouped density plot (pdf p11)
* Bar chart with error bars (only makes sense for normally distributed data) (pdf p14)

Two-sample t-test (pdf p15)
-----------------
(is the sample mean of two samples significantly differen?)

    1. Independent samples (pdf p17)
    ----------------------
        Testing procedure:
        a.  Hypotheses:
            * H0: µ1 - µ2 = 0
            * H1: µ1 - µ2 < 0
        b.  Significance level: e.g. alpha = 0.05
        c.  Test statistic:
            * _x - _y
            * _x = estimation for µ1 (control group)
            * _y = estimation for µ2 (intervention group)
        d.  Calculate p
        e.  Conclusion by compare p with the significance level alpha
                p > alpha: do not reject  <=>     p < alpha: reject H0

        Calculation in Python:
        control = np.array(*[CONTROL VALUES]*)
        intervention = np.array(*[INTERVENTION VALUES]*)
        stats.ttest_ind(a=control, b=intervention, alternative='less', equal_var=False)

            Result:
            Ttest_indResult(statistic= , pvalue= )

    2. Paired samples (pdf p23)
    -----------------
        Testing procedure:
        a.  Hypotheses:
            * H0: _(x - y) = 0
            * H1: _(x - y) < 0
        b.  Significance level: e.g. alpha = 0.05
        c.  Test statistic:
            * _(x - y)
            * x = estimation for µ1 (control group)
            * y = estimation for µ2 (intervention group)
        d.  Calculate p
        e.  Conclusion by compare p with the significance level alpha
                p > alpha: do not reject  <=>     p < alpha: reject H0

        Calculation in Python:
        regular = np.array(*[REGULAR VALUES]*)
        additives = np.array(*[ADDITIVES VALUES]*)
        stats.ttest_rel(regular, additives, alternative='less')

            Result:
            Ttest_relResult(statistic= , pvalue= )

Effect size (pdf p28)
-----------
= is a metric which expresses how great the difference between two groups is

    * Control group   <=>     intervention group
    * can be used in addition to hypothesis test
    * Several definitoins, here => Cohen's d

    Cohen's D: (pdf p29)
    ----------

        _x1 - _x2
    d = ---------
            s

    with:
        * _x1, _x2, the sample means
        * s, standard deviation of both groups combined:

                ((n1 - 1) * s1² + (n2 - 1) * s2²)
        s = sqrt(-------------------------------)
                (           n1 + n2 - 2         )

        with:
            * n1, n2, the sample sizes
            * s1, s2, the standard deviation of both groups
    
    Interpretation of Cohen's d: (pdf p30)

    |--------|-----------------|
    |   |d|  |   Effect Size   |
    |--------|-----------------|
    |   0.01 |   Very small    |
    |   0.2  |   Small         |
    |   0.5  |   Average       |
    |   0.8  |   Large         |
    |   1.2  |   Very Large    |
    |   2.0  |   Huge          |
    |--------|-----------------|

    Interpretation for educational sciences, see pdf p30-31
'''
)

samenvatting[6] = (
'''
Bivariate analysis:
*******************
QUANTITATIVE - QUANTITATIVE variate:
------------------------------------

Data Visualization (pdf p6):
------------------
Scatter plot:
* X-axis: independent variable
* Y-axis: dependent variable
* Each point corresponds to an observation

Linear Regression (pdf p8)
-----------------
With regression we will try to find a CONSISTENT and  SYSTEMATIC relationship between two QUANTITATIVE variables
1.  Monotonic: consistent direction to the relationship between the two variables: 
        INCREASING   <=>     DECREASING
2.  Non-monotonic: value of dependent variable changes systematically with value of independent variable, but the direction is not consistent

=> Linear regression = linear relationship between an independent and dependent variable

Characteristics:
    * Presence: is there a relationship?
    * Direction: increasing or decreasing?
    * Strength of the relationship: strong, moderate, weak, nonexistent,...

Method of least squares (pdf p16)

Covariance: (pdf p19)
-----------
= a measure that indicates whether a relationship between two variables is increasing or decreasing

Cov(X, Y) = 1 /(n -1) * sum((x_i - _x)(y_i - _y))

* Cov > 0: increasing
* Cov +- 0: no relationship
* Cov < 0: decreasing
***NOTE: Covariance of population (denominator n) <=> sample (denominator n - 1)***

Pearson's correlation coefficient (pdf p25)
---------------------------------
Pearson's product-moment correlation coefficient R is a measure for the strength of a linear correlation between x and y

    Cov(X, Y)
R = ---------
    σ_x * σ_y

R ∈ [-1, 1]

Coefficient of determination (pdf p28)
----------------------------
= R²: explains the percentage of the variance of the observed values relative to the regression line

R²: percentage variance observations explained by the regression line
1 - R²: percentage variance observatoins NOT explained by regression

Interpretation of R and R² values (pdf p29)
---------------------------------

    |---------------|--------------|------------------------|-------------------|
    |   |R|         |   R²         |    Explained variance  |   Interpretation  |
    |---------------|--------------|------------------------|-------------------|
    |   < 0.3       |   < 0.1      |    < 10%               |   very weak       |
    |   0.3 - 0.5   |   0.1 - 0.25 |    10 - 25%            |   weak            |
    |   0.5 - 0.7   |   0.25 - 0.5 |    25 - 50%            |   moderate        |
    |   0.7 - 0.85  |   0.5 - 0.75 |    50 - 75%            |   strong          |
    |   0.85 - 0.95 |   0.75 - 0.9 |    75 - 90%            |   very strong     |
    |   > 0.95      |   > 0.9      |    > 90%               |   exceptional(!)  |
    |---------------|--------------|------------------------|-------------------|

Considerations:
    * The correlation coefficient only looks at the relationship between two variables. Interactions with other variables are not considered
    * The correlation coefficient explicitly does not assume a causal relationship
    * Pearson's correlation coefficient only expresses linear relationships
'''
)

samenvatting[7] = (
'''
Time series and predictions (pdf p5)
***************************
Time series is a sequence of observations of some variable over time
Many decisions in business operations depend on a forecast of some quantity

Time series are a statistical problem: observations vary with time

Time series components (pdf p8)
----------------------
* Level
* Trend
* Seasonal fluctuations
* Cyclic patterns
* Random noise (residuals)

Time series models (pdf p10)
------------------
The simplest mathematical model:
    X_t = b + epsilon_t                   (1)

    * X_t = estimate for time series, at time t
    * b = the level (a constant), based on observations x_t
    * epsilon_t = random nois. We assume that epsilon_t +- Nor(µ = 0; sigma)

We could also assume that there is a linear relationship:
    X_t = b_0 + b_1 * t + epsilon_t       (2)

    * b_0 = level
    * b_1 = trend

Equation 1 and 2 are special cases of the polynomial case:
    X_t = b_0 + b_1 * t + b_2 * t² + ... + b_n * t ^ n + epsilon _t     (3)

General expression time series: (pdf p13)
    X_t = f(b_0, b_1, b_2, b_n, t) + epsilon_t                          (4)

    Assumptions:
    * Two components of variability:
        1. the mean of the predictions changes with time
        2. the variations to this mean vary randomly
    * The residuals of the model (X_t - x_t) have a constant variance in time (homoscedastic)

Estimating the parameters:  (pdf p14)
    Make predictions based on the time series model:
        1. select the most suitable model
        2. estimation for parameters b_i(i: 1,..., n) based on observations
    The estimations b_i are selected so that they approximate the observed values as close as possible

Moving average (pdf p19)
--------------
= is a series of averages (means) of the last m observations

    * Notation: SMA
    * Hide short-term fluctuations and show long-term trends
    * Parameter m is the time window

    SMA(t) = sum(x_i / m) from i=k to t

    with k = t - m + 1

Weighted moving average (pdf p22)
-----------------------
    * For SMA, the weights of the observations are equal
    * For a weighted moving average (WMA), more recent observations gain relatively more weight
    * A specific form of this is single exponential smoothing or the
      EXPONENTIAL MOVING AVERAGE (EMA)

      X_t = alpha * x_t-1 + (1 - alpha) * X_t-1                   (6)

      With alpha the smoothing constant (0 < alpha < 1), and t >=3

      This is only valid from t = 3. So for value x_2 must choose suitable value:
        * x_2 = x_1
        * x_2 = 1/m sum(x_i) (so the mean of the first m observations)
        * x_2 equal to a specific objective

Exponential smoothing (pdf p26)
---------------------
Older observations have an exponentially smaller weight

The speed at which the old observations are "forgotten" depends on the value of alpha.
Alpha close to 1    =>  old observations are quickly forgotten
Alpha close to 0    =>  old observations are take more time to ooze away

Forecasting: (pdf p28)
Forecast for time t + m (m time units in the "future") take always the last estimate of the level:
    F(t + m) = X_t

Double Exponential Smoothing (aka Holt-method) (pdf p29)
----------------------------
Basic exponential smoothing does not work well if there is a trend in the data, the errors keep getting bigger
=> add additional term to model the trend: b_t for the estimation of the trend at time t > 1:

        X_t = alpha * x_t-1 + (1 - alpha) * (X_t-1 + b_t-1)
        b_t = beta * (X_t - X_t-1) + (1 - beta) * b_t-1

        with 0 < aplha < 1 and 0 < beta < 1
        * b_t is an estimate for the slope of the trend line
        * added to the first equation to ensure that the trend is followd
        * X_t - X_t-1 is positive or negative, this corresponds to an increasing/decreasing trend

        Options for selecting intial values:
        * X_1 = x1
        * b_1 = X_2 - x_1
        * b_1 = 1/3 [(x_2 - x_1) + (x_1 - x_2) + (x_4 - x_3)]
        * b_1 = (x_n - x_1) / (n - 1)

Predicting (forecasting) (pdf p32)
Forecast for time t + m:
    F(t + m) = X_t + m * b_t

Triple Exponential Smoothing (aka Holt-Winter's Method) (pdf p34)
----------------------------
If there is a recurring (seasonal) pattern

Notation:
    * L: length of the seasonal cycle (number of time units)
    * c_t: term that models the seasonal variations
    * gamma: smoothing factor for the seasonal variation

    X_t = alpha * (X_t / C_t-L) + (1 - alpha)(X_t-1 + b_t-1)        Smoothing
    b_t = beta * (X_t - X_t-1) + (1 - beta) * b_t-1                 Trend Smoothing
    c_t = gamma * (x_t/X_t) + (1 - gamma) * c_t-L                   Seasonal Smoothing

Prediction: (pdf p36)
Prediction at time t + m:
    F(t + m) = (X_t + m * b_t) * c_t-L+m

Quality of a time series model (pdf p38)
------------------------------
Compare forecast results with actual observations, when they become available:
    * Mean absolute error: MAE = 1 / m * sum(abs(x_i - F_i)) for i = t + 1 to t + m
    * Mean squared error: MSE = 1 / m * sum((x_i - F_i)²) for i = t + 1 to t + m

If square root of MSE is well below standard deviation over all observations, you have a good model!
'''
)

samenvatting[8] = (
'''
H3:
* Right or left tailed z-test: (oefening: confituurpotten; controleren van de instelling van de machine correct is)
  t-test if sample is too small or standard deviation is unknown

H4:
* Chi²:
    1. Chi² for independence: (oefening: antwoorden op enquête verschillende voor M als voor V)
    2. Chi² goodness of fit: (oefening: komt de sample overeen met de theoretisch te verwachten waarde)

H5:
* two sample t-test
    1. independent: reactiesnelheid na innemen van medicijn
    2. paired: benzine met additieven

H6:
* regressie (flipperlenght ifv body mass) + CoV + R + R²

H7:
* time series
'''
)

def help():
    commands = {'searchString(str)','overzichtHs()','printTheorie()','printH(hoofdstuk)','overzichtTs()'}

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

def printH(hoofdstuk):
    
    print(samenvatting[hoofdstuk])

def overzichtTs():
    print(samenvatting[8])

