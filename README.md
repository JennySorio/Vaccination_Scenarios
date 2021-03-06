# Vaccination Scenarios

Program in Matlab associated with the generation of vaccination scenarios against COVID-19.

This program was carried out as part of the research project on epidemiological models of COVID-19 and the impact of vaccination in this pandemic. 

The project members are Eduardo Massad, Jorge P. Zubelli, Vinicius Albani and Jennifer Loria.

## Repository contents

There are two folders: Chicago and New York. In the New York folder there are two types of program;  one that considers the age range of the population and another in which it is not considered.  

Each folder also has the data used from 01-Mar-2020 to 28-Nov-2020 (https://www.chicago.gov/city/en/sites/covid-19/home.html and https://www1.nyc.gov/site/doh/covid/covid-19-data.page). The programs present in each one are very similar. The main differences are related to the data of each one.

To create part of the code of the new model proposed in this project, we took as an initial guide some known codes of the SEIR model (https://cs.uwaterloo.ca/~paforsyt/SEIR.html).

## How to use this repository

The objectives of this program are:

1. Estimate predictions about cases, hospitalizations and deaths associated with Covid-19.
2. Estimate different scenarios where the vaccination campaigns started on different dates.

To meet the first objective, we must use the files that start with the name **Bootstraping** and **mySEIR**, which, estimate the predictions in addition they estimate all the parameters associated with the model. We suggest to run first the files **mySEIR** because they contain the portion of codes where the datas are imported.

To fulfill the second objective, we use the files that start with the name **Vaccination** and **EvaluatingPaths**.

The functions that are in the program whose name begins with **seir**, are the functions that generate as output some parameters defined in the SEIR-like model proposed in the article "The Impact of Covid-19 Vaccination Delay: A Modelling Study for Chicago and NYC Data". 

The functions that start with **ObjFun** are those related to the generation and calibration of the vector beta. For this, are estimated the mildly infective population  and b.

The function **basic_reproduction_rate_beta** generates the time-dependent effective reproduction number.


## Contacts

* Vinicius Albani: v.albani@ufsc.br

* Jennifer Loría: jennyls@impa.br

