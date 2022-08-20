# Importing the necessary packages
import numpy as np  # "Scientific computing"
import scipy.stats as stats  # Statistical tests
import pandas as pd  # Data Frame
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt  # Basic visualisation
from statsmodels.graphics.mosaicplot import mosaic  # Mosaic diagram
import seaborn as sns  # Advanced data visualisation
import altair as alt  # Alternative visualisation system
import math
from scipy.stats.contingency import association
from matplotlib import cm
from sklearn.linear_model import LinearRegression
import sys
from statsmodels.tsa.holtwinters import ExponentialSmoothing


import matplotlib
import statsmodels
print(matplotlib.__version__)
print(statsmodels.__version__)
if(matplotlib.__version__ !='3.4.3'):
  !pip install matplotlib==3.4.3
  !pip install --upgrade statsmodels
  
!if cd DataScience-AI; then git pull; else git clone 'https://github.com/RobinBenoot2/DataScience-AI.git'; fi

path_m = 'DataScience-AI'

sys.path.insert(0,path_m)
import theorie as th
importlib.reload(th)

# --------------------------H 1 - 2 --------------------------
def printCommands():
    commands = ('set index: df.set_index([id])',
                'Drop column: df.drop("column name", axis="columns")',
                'Select adjacent columns: df.iloc[:, 2:4]',
                'Select multiple columns: df[["col1", "col2]]',
                'Observation with row number 5 (counting from zero): df.iloc[5])',
                'The first 4 observations df.iloc[0:4]',
                'Select observations where the value of Age is less than 18 df[col < 18]',
                'Select all boys younger than 10: df.query("(Sex=="male") and (Age < 18)")',
                'cleaned = df.dropna() # Drop any row that has at least one missing value',
                'cleaned = df.dropna(how="all") # Drop any row that has all missing values',
                'df = df.fillna(value={Age : avg_age, Cabin : AAA})',
                'air_quality["london_mg_per_cubic"] = air_quality["station_london"] * 1.882',
                'Mode for all the variables in the DataFrame: tips.mode()'
                'Mode for a specific variable: tips["day"].mode()'
                )
    for i in commands:
        print(i)


def connectDrive():
    from google.colab import drive
    drive.mount('/content/drive')


def openCsv(path, seperator):
    # /content/drive/MyDrive/ds/dsai-en-labs-main/data/ais.csv
    source = pd.read_csv(path, sep=seperator)
    source.head()
    return source


def printGeneralInfo(source):
    print(source.head())
    print(f"Number of rows: {len(source)}")
    print(f"Number of columns: {len(source.columns)}")
    print(f"The shape of the Data Frame is: {source.shape}")
    print("*" * 50)
    print(source.info())
    print("*" * 50)
    print(source.dtypes)
    print("*" * 50)
    print(source.dtypes.value_counts())
    print("*" * 50)
    print(source.describe())
    print("*" * 50)
    print(source.count())
    print("*" * 50)
    print(source.value_counts())


def printSpecificInfo(source, col):
    # Centrality and dispersion deasures
    # Mean, standard deviation & friends
    print(f"Mean:                {source[col].mean()}")
    print(f"Standard deviation:  {source[col].std()}")  # Pay attention: n-1 in the denominator
    print(f"Variance:            {source[col].var()}")  # Pay attention: n-1 in the denominator
    print(f"Skewness:            {source[col].skew()}")
    print(f"Kurtosis:            {source[col].kurtosis()}")
    # Median & co
    print(f"Minimum:   {source[col].min()}")
    print(f"Median:    {source[col].median()}")
    print(f"Maximum:   {source[col].max()}")
    percentiles = [0.0, 0.25, 0.5, 0.75, 1.0]
    print("Percentiles", percentiles, "\n", source[col].quantile(percentiles))
    print("Inter Quartile Range:", source[col].quantile(.75) - source[col].quantile(.25))
    print(f"Range :    {source[col].max() - source[col].min()}")


def printColumnInfo(source):
    print(source.unique())
    print("*" * 50)
    print(source.describe())


def convertColToQualitative(source, objectList):
    for o in objectList:
        source[o] = source[o].astype('category')


def convertToOrdinal(source, orderArray):
    return source.astype(CategoricalDtype(categories=orderArray, ordered=True))


def seabornPlots():
    # Bar chart in Seaborn: catplot() with 'kind = "count"''
    print('sns.catplot(data = ap, kind="count",  x = "DataSize", hue="PersistenceType");')
    # Visualisation using a box plot (Seaborn)
    print('sns.boxplot(data=ap, x="Time", y ="DataSize", hue="PersistenceType")'),
    # Violin plot (Seaborn)
    print('sns.violinplot(data = tips, x = "tip");'),
    # probability density (kde or kernel density estimation).
    print('sns.kdeplot(x = tips["tip"]);'),
    # histogram
    print('sns.displot(x = tips[col], kde=True);'),
    print('g = sns.FacetGrid(data=ap, col="DataSize", hue="PersistenceType")  g.map(sns.kdeplot,"Time")  g.add_legend()')


def groupByExample():
    print('''cols=['mean','std']
    df = ap['Time']
    df.describe()[cols]
    df = ap.groupby('DataSize')['Time'].describe()[cols]
    df2 = ap.groupby(['DataSize', 'PersistenceType'])['Time'].describe()[cols]''')


# -------------------------- H2 --------------------------


def standardNormalDist():
    # Take 200 values for the X-axis, between -4 and 4, evenly spaced
    x = np.linspace(-4, +4, num=201)
    y = stats.norm.pdf(x, 0, 1)
    # Plot the probability density function (pdf) for these X-values
    plt.plot(x, y)


def normalDist(m, s):
    x = np.linspace(m - 4 * s, m + 4 * s, num=201)
    plt.plot(x, stats.norm.pdf(x, loc=m, scale=s))

# -------------------------- Module 3 - The Central Limit Theorem --------------------------
# Denk aan een sample nemen in een fabrieksomgeving om het bvb het gewicht in te schatten per confituurpot
# 1. Formulate both hypotheses (𝐻0 and 𝐻1)
# 2. Determine the significance level (𝛼)
# 3. Calculate the test statistic
# 4. Determine the critical region or the probability value
# 5. Draw conclusions

def leftAndRightTailedProbability(val, m, s):
    cdf = stats.norm.cdf(val, loc=m, scale=s)
    sf = stats.norm.sf(val, loc=m, scale=s)
    isf = stats.norm.isf((1 - .8), loc=m, scale=s)
    # In Python, you can use the SciPy-function norm.cdf()
    # to calculate the left tail probability  P(X<x)  or  P(Z<z)  (also called the cumulative distribution).
    print(f"left tail cdf: {cdf:.3f}")
    # In order to calculate the right tail probability, we use the norm.sf() function
    # (defined as 1 - cdf, also called the survival function, hence the function name)
    print(f"right tail sf: {sf:.3f}")
    # Another type of question: under what value will 80% of observations be?
    # To calculate this we would need the inverse function of cdf().
    # However, it does not exist in SciPy. We do have the inverse function of sf(),
    # though, which is called isf().
    # We can find the result by calculating the reaction time above which 20% of the values lie:
    print(f"probability for 80% of values: {isf:.3f}")
    twoTailedProbability(cdf, sf, m, s)
    return cdf, sf, isf


def twoTailedProbability(high, low, m, s):
    if high <= 1 and low <= 1:
        high = stats.norm.isf(high, loc=m, scale=s)
        low = stats.norm.isf(low, loc=m, scale=s)
    if high < low :
        temp = low
        low = high
        high = temp
    # In Python, you can use the SciPy-function norm.cdf()
    # to calculate the left tail probability  P(X<x)  or  P(Z<z)  (also called the cumulative distribution).
    print(f"two tailed cdf: {stats.norm.cdf(high, loc=m, scale=s) - stats.norm.cdf(low, loc=m, scale=s):.3f}")
    # X-values
    dist_x = np.linspace(m - 4 * s, m + 4 * s, num=201)
    # Y-values for drawing the Gauss curve
    dist_y = stats.norm.pdf(dist_x, m, s)
    # Plot the Gauss-curve
    plt.plot(dist_x, dist_y)
    #plt.plot(dist_x, stats.norm.cdf(dist_x, m, s))
    # Fill the area left of x
    plt.fill_between(dist_x, 0, dist_y, where=(dist_x <= high) & (dist_x > low), color='lightblue')
    # Show x with a green line
    plt.axvline(stats.norm.isf(low, loc=m, scale=s), color="green")
    plt.axvline(stats.norm.isf(high, loc=m, scale=s), color="green")


def getParams(source, percentage):
    n = source.count()
    m = source.mean()
    s = source.std()
    alpha = 1 - percentage / 100  # 1 - alpha is the confidence level
    return n, m, s, alpha


# So we state with a confidence level of 95% that the reaction speed of the superheroes is between 4.37 and 6.03 ms.
def confidenceIntervals(n, m, s, percentage, populationStdKnown):
    # n = sample count, m = sample mean, s = population std of sample std
    alpha = 1 - percentage / 100  # 1 - alpha is the confidence level
    # We then find the  z  score between which 95% of all values lie with a standard normal distribution.
    if n > 29 and populationStdKnown:
        res = stats.norm.isf(alpha / 2)
        print("z-score: %.3f" % res)
    else: # t-test when std (or var) of population is unknown or sample size < 30
        res = stats.t.isf(alpha / 2, df=n - 1)
        print("t-score: %.3f" % res)
    # We use this to determine the values to the left and right of the sample mean
    # between which we expect 95% of the values to fall for the probability distribution
    # that we get from the central limit theorem.
    lo = m - res * s / np.sqrt(n)
    hi = m + res * s / np.sqrt(n)
    print("The region of values supporting the null hypothesis = Confidence interval: [%.3f, %.3f]" % (lo, hi))
    # X-values
    twoTailedProbability(hi, lo, m, s)


# -------------------------- Module 4: Bivariate analysis - 2 qualitative variables --------------------------
# Chi-squared ( χ2 ) and Cramér's V:
# are statistics that can help us to determine whether there is an association between two qualitative (categorical) variables.
# The chi-squared test for independence :
# To answer the question of when the value of chi-square is sufficient to assume an association between two variables, we can use the chi-square independence test.


def convertPopulationStdToSampleStd(populationStd, sampleSize):
    return populationStd / math.sqrt(sampleSize)


def calculateZscoreFromSampleAverage(sampleMean, m, s):
    print("z: %.3f" % (sampleMean - m) / s)


# The  z -test is used to confirm or refute an assumption about the (unknown) population mean, based on a sufficiently large sample.
def zTest(h0, h1, n, sm, s, alpha, m0):
    # Properties of the sample:
    # n             # sample size
    # sm            # sample mean
    # s             # population standard deviation (assumed to be known)
    # a             # significance level (chosen by the researcher)
    # m0            # hypothetical population mean (H0)
    cdf, sf, isf = leftAndRightTailedProbability(sm, m0, s / np.sqrt(n))
    if h0 < h1 and n >= 30:
        print('Using right tailed z-test')
        p = sf
        g = stats.norm.isf(alpha, loc=m0, scale=s / np.sqrt(n))
    elif h0 > h1 and n >= 30:
        print('Using left tailed z-test')
        p = cdf
        g = stats.norm.isf(1-alpha, loc=m0, scale=s / np.sqrt(n))
    elif h0 < h1 and n < 30:
        print('Using right tailed t-test')
        p = stats.t.sf(sm, loc=m0, scale=s/np.sqrt(n), df=n-1)
        g = stats.t.isf(alpha, loc=m0, scale=s/np.sqrt(n), df=n-1)
    elif h0 > h1 and n < 30:
        print('Using left tailed t-test')
        p = stats.t.cdf(sm, loc=m0, scale=s/np.sqrt(n), df=n-1)
        g = stats.t.isf(1 - alpha, loc=m0, scale=s/np.sqrt(n), df=n-1)
    print("p-value: %.3f" % p)
    print("Critical value g ≃ %.3f" % g)
    if p < alpha:
        print("p < a: reject H0")
        print('The probability value is smaller than the significance level, so we can reject the null hypothesis.')
    else:
        print("p >= a: do not reject H0")
        print('The probability value is equal or higher than the significance level, so we cannot reject the null hypothesis.')
    if sm < g and h0 < h1:
        print("sample mean = %.3f < g = %.3f: do not reject H0" % (sm, g))
        print('x is inside the normal region, so we cannot reject the null hypothesis.')
    elif sm > g and h0 < h1:
        print("sample mean = %.3f > g = %.3f: reject H0" % (sm, g))
        print('x  is inside the critical region, so we can reject the null hypothesis. Therefore, we can assume that binding recommendation on continuation of studies does increase the success rate.')
    elif sm > g and h0 > h1:
        print("sample mean = %.3f > g = %.3f: do not reject H0" % (sm, g))
        print('x is inside the normal region, so we cannot reject the null hypothesis.')
    elif sm < g and h0 > h1:
        print("sample mean = %.3f < g = %.3f: reject H0" % (sm, g))
        print('x  is inside the critical region, so we can reject the null hypothesis. Therefore, we can assume that binding recommendation on continuation of studies does increase the success rate.')
    g_under = stats.norm.isf(1-alpha/2, loc=m0, scale=s/np.sqrt(n))
    g_upper = stats.norm.isf(alpha/2, loc=m0, scale=s/np.sqrt(n))
    print("Confidence interval: [%.3f, %.3f]" % (g_under, g_upper))
    plot(m0, s/np.sqrt(n), sm, g)


def plot(m, s, sm, g):
    # Gauss-curve plot:
    # X-values
    dist_x = np.linspace(m - 4 * s, m + 4 * s, num=201)
    # Y-values for the Gauss curve
    dist_y = stats.norm.pdf(dist_x, m, s)
    fig, dplot = plt.subplots(1, 1)
    # Plot the Gauss-curve
    dplot.plot(dist_x, dist_y)
    # Show the hypothetical population mean with an orange line
    dplot.axvline(m, color="orange", lw=2)
    # Show the sample mean with a red line
    dplot.axvline(sm, color="red")
    # Fill the acceptance area in light blue
    dplot.fill_between(dist_x, 0, dist_y, where=dist_x <= g, color='lightblue')
    plt.show()
    plt.clf()


def cramersV(table):
    val = association(table, method="cramer")
    print('Cramers v: %.3f' % val)
    if val < 0.1:
        return 'no association'
    elif val < 0.3:
        return 'weak association'
    elif val < 0.5:
        return 'moderate association'
    elif val < 0.8:
        return 'strong association'
    elif val < 0.99:
        return 'very strong association'
    else:
        return 'complete association'


# Chi-squared ( χ2 ) and Cramér's V are statistics that can help us to determine whether there is an association between two qualitative (categorical) variables.
# The reasoning goes as follows:
# if there is no association between Gender and Survey,
# then we expect the proportions of the values of Survey to be the same for all values of Gender.
# In other words, for both women and men, the same percentage of respondents will give the same answer to the question.
def qualitativeIndependence(observed):
    chi2, p, df, expected = stats.chi2_contingency(observed)
    print('Chi-squared : %.3f' % chi2)
    print('Degrees of freedom: %d' % df)
    print('the number q for which the right tail probability is exactly 5 percent: %.3f' % stats.chi2.isf(0.05, df))
    print('P-value : %.3f' % p)
    print('Cramers v: %s' % cramersV(observed))


# is this sample representative of the population? Does each type occur in the sample in proportion to the expected percentage in the population as a whole?
# komt de sample overeen met de theoretisch te verwachten waarde
def goodnessOfFitTest(arr_observed, arr_expected):
    alpha = 0.05               # Significance level
    n = sum(arr_observed)          # Sample size
    k = len(arr_observed)          # Number of categories
    dof = k - 1                # Degrees of freedom
    expected = arr_expected * n  # Expected absolute frequencies in the sample
    g = stats.chi2.isf(alpha, df=dof)  # Critical value
    # Goodness-of-fit-test in Python:
    chi2, p = stats.chisquare(f_obs=arr_observed, f_exp=expected)
    print("Significance level  ⍺ = %.3f" % alpha)
    print("Sample size         n = %d" % n)
    print("k = %d; df = %d" % (k, dof))
    print("Chi-squared        χ² = %.3f" % chi2)
    print("Critical value      g = %.3f" % g)
    print("p-value             p = %.3f" % p)


# -------------------------- Module 5 - Bivariate analysis: qualitative vs quantitative --------------------------
# In this module, we discuss the case where the independent variable is qualitative and the dependent is quantitative. Some typical examples of research question in this case:
# Within a certain species of animals, are male individuals significantly larger than females?
# Does a new vaccine protect against the disease like it's supposed to?
# Does a certain study method like "retrieval practice" actually improve learning outcomes (i.e. student's grades)?


#Effect size
# If we want to know whether two groups are significantly different,
# we can use a statistical test like the two sample  t -test.
# The result of a statistical test is generally either "true" or "false",
# depending on the  p -value and the chosen significance level
def cohen_d(a, b):
    na = len(a)
    nb = len(b)
    pooled_sd = np.sqrt(((na-1) * a.std(ddof=1)**2 +
                          (nb-1) * b.std(ddof=1)**2) / (na + nb - 2) )
    res = (b.mean() - a.mean()) / pooled_sd
    print('Cohen_d: %.3f' % res)
    var = math.fabs(res)
    if var < 0.05:
        return 'Very small'
    elif var < 0.25:
        return 'Small'
    elif var < 0.6:
        return 'Average'
    elif var < 0.9:
        return 'Large'
    elif var < 1.4:
        return 'Very large'
    else:
        return 'Huge'


# de onafhankelijke t-test bij een kwantitatieve en kwalitatieve vergelijking
# returns statistic, pvalue
def independent_Ttest(a, b):
    # Parameters a and b are the two groups to be compared.
    # alternative='less' indicates that we want to test for the alternative hypothesis that the mean of the control group is less than the mean of the treatment group.
    # Finally, by setting equal_var=False, it is not assumed that both groups have the same standard deviation.
    print(stats.ttest_ind(a=a, b=b, alternative='less', equal_var=False))
    return stats.ttest_ind(a=a, b=b, alternative='less', equal_var=False)


#The t-test for paired samples
# In this variant of the  t -test, a measurement is taken on each element of the sample, one time before and one time after an intervention.
# The aim is to determine whether the intervention had a significant effect on the measurement.
def dependent_Ttest(a, b):
    print(stats.ttest_rel(a, b, alternative='less'))    # kan ook alt: 'more' of two-sided' zijn!
    return stats.ttest_rel(a, b, alternative='less')


# -------------------------- Module 6. Bivariate analysis of 2 quantitative variables --------------------------


#params: two columns
def linearRegression(colA, colB):
    x = colA.values.reshape(-1, 1)
    y = colB
    weight_model = LinearRegression().fit(x, y)
    print(f"Regression line: ŷ = {weight_model.intercept_:.3f} + {weight_model.coef_[0]:.3f} x")
    sns.regplot(x=colA, y=colB)
    return weight_model


# Depending on the value of  |R|  (or  R2 ), you can draw a conclusion about the strength of the linear relation between the two variables:
def correlationCoefficient(colA, colB):
    cor = np.corrcoef(colA, colB)[0][1]
    Coefficient_of_determination = cor**2
    print(f"R = {cor}")
    print(f"R² = {Coefficient_of_determination}")
    if Coefficient_of_determination < 0.1:
        print('Very weak linear relation between the two variables')
    elif Coefficient_of_determination < 0.25:
        print('Weak linear relation between the two variables')
    elif Coefficient_of_determination < 0.5:
        print('Moderate linear relation between the two variables')
    elif Coefficient_of_determination < 0.75:
        print('Strong linear relation between the two variables')
    elif Coefficient_of_determination < 0.9:
        print('Very strong linear relation between the two variables')
    else:
        print('Exceptionally strong linear relation between the two variables')

def movingAverage(df, col, val):
    df[col + '_' + str(val)] = df[col].rolling(val).mean()
    df.plot( y=[col, col + '_' + str(val)], figsize=[10,5])


def holtWinters(source, col, range):
    train = source['col'][:range]
    test = source['col'][range:]
    model = ExponentialSmoothing(train, trend='add', seasonal='mul', seasonal_periods=12, freq='MS').fit()
    train.plot(legend=True, label='train')
    test.plot(legend=True, label='test')
    model.fittedvalues.plot(legend=True, label='fitted')