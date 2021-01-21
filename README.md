# Vaccination Scenarios

Program in Matlab associated with the generation of vaccination scenarios against COVID-19.

This program was carried out as part of the research project on epidemiological models of COVID-19 and the impact of vaccination in this pandemic. Its members are: Eduardo Massad, Jorge P. Zubelli, Vinicius Albani and Jennifer Loria.

## Repository contents

There are two folders: Chicago and NYC (New York City). The programs present in each one are very similar. The main differences are related to the data format of each one.

In each folder we have two types of programs; one that considers the age range of the population and another in which it is not considered. Alem, each folder also has the data used (fuente bibliografica). 
FECHAS DE LOS DATOS.

* Bootstrapping_20201027:
* Bootstrapping_20201104:
* Bootstrapping_20201204:
*
*
*

## How to use this repository

The objectives of this program are:

1. Estimate predictions about cases, hospitalizations and deaths associated with Covid-19.
2. Estimate different scenarios where the vaccination campaigns started on different dates.

To meet the first objective, we must use the files that start with the name **Bootstraping** and **mySEIR**, which, in addition to estimating the predictions, estimates all the parameters associated with the model.

To fulfill the second objective, we use the files that start with the name **Vaccination** and **EvaluatingPaths**.

The functions that are in the program whose name begins with **seir**, are the functions that generate as output the derivatives defined in the SEIR-like model proposed in the article "The Impact of Covid-19 Vaccination Delay: A Case Study with Chicago and NYC Data".

Taking between its inputs one vector with the population susceptible (S), vaccinated (V), exposed (E), asymptomatic and infective (IA), mildly infective (IM), severely infective (IS), critically infective (IC), recovered (R), and deceased (D).

The functions that start with **ObjFun** are those related to the generation of the beta rate (Explicar quien es la tasa beta!!!)


## Contacts

* Vinicius Albani: v.albani@ufsc.br

* Jennifer Loría: jennyls@impa.br

