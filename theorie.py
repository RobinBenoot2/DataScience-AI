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
