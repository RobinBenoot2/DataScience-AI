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
