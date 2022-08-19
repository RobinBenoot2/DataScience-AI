import pandas as pd

def help():
    commands = {'printTheorie()','printH1'}

    for i in commands:
        print(i)

def printTheorie():
    printH1()

#-------H1------

def printH1():
    print('Variable = general propery of an object, allows to distinguish objects')
    print('Value = specific propery, interpretation for that variable\n\n')

    print("\u0332".join("Measurement Levels"))

    print('Qualtitative = Limited number of values, not necessarily numeric (bv: small, medium, large)')
    print('Quantitative = Many values, often unique. Number + unit of measurement (often contain result of measurement')


