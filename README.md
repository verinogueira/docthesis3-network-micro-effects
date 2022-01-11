# docthesis3-network-micro-effects: README
Replicate analysis of "Production network structures and micro effects of macro shocks,'' the third chapter of my doctoral thesis available at the University of Essex's [Repository](http://repository.essex.ac.uk/31206/). 

****
WORK IN PROGRESS: ignore below here
****



Data sources are Lee and Lee (2016), Barro and Lee (2001), The World Bank (2017), Feenstra, Inklaar, and Timmer (2015). The model is calibrated using parameters from Bils and Klenow (2000), Montenegro and Patrinos (2014), Acemoglu and Autor (2011) and Ciccone and Peri (2005). All sources are presented in the paper too. 

The construction of the dataset, analysis and exhibits use Stata. Code for data cleaning and analysis is provided as part of the replication package. The main file runs all of the code to generate the data for the 10 figures and 3 tables in the paper. The replicator can expect the code to run for a couple of minutes.


--------
## Data Availability and Provenance Statements

### Rights, license and availability

- [X] I certify that the author of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [X] The data are licensed under a GNU General Public License v3.0. See LICENSE.txt for details.
- [X] All data **are** publicly available.


### Data Sources

#### UK Data Service

The paper uses data from the Annual Population Survey Household Dataset, which is distributed by the UK Data Service (Office for National Statistics, 2015). Data is subject to a redistribution restriction, but can be freely downloaded from https://ukdataservice.ac.uk/. Search for `SN 7750` then click at `Access data`. Fill out the registration form, including a brief description of the project, and agree to the conditions of use. Note: "the data files themselves are not redistributed" and other conditions. Save the file in the directory `RawData`. 

Datafile: `RawData/apsp_jd14_eul_pwta17.dta` (not provided)


Office for National Statistics, Social Survey Division. (2015). Integrated Household Survey, January - December, 2014. [data collection]. UK Data Service. SN: 7839, DOI: 10.5255/UKDA-SN-7839-1




#### Lee and Lee (2016)

Data on school enrolment comes from Lee and Lee (2016), hereafter Lee-Lee, which
covers 111 countries in five-year intervals, from 1820 to 2010. The version 1.0 used
in this chapter comes from the most recent update available. The sample selected
is restricted to the male population to comply with the Mincerian returns used. 

Data can be downloaded from https://github.com/barrolee/BarroLeeDataSet/blob/master/LeeLee/LeeLee_v1.dta. 
Save the file in the directory `Input/`. 

Datafile: `Input/LeeLee_v1.dta` (not provided)



#### Barro and Lee (2001)

The bulk of the data on schooling comes from Barro and Lee (2001), hereafter Barro-
Lee, presenting educational data for 146 countries, ranging from 1950 to 2010. The
version used is the 2.1. Once more, only the male population subsample was selected.

Data can be downloaded from https://barrolee.github.io/BarroLeeDataSet/BLData/BL2013_M_v2.2.dta. 
Save the file in the directory `Input/`. 

Datafile: `Input/BL2013_M_v2.2.dta` (not provided)

#### World Bank Indicators

Additionally, information on the theoretical duration of educational levels was retrieved
from the World Bank Indicators platform (The World Bank, 2017). The data
are provided by the United Nations Educational, Scientific, and Cultural Organization
(UNESCO) Institute for Statistics is available yearly for 247 countries since 1970. I use
the variable for the number of grades (years) in secondary education (Duration).


Secondary education, duration (years): https://data.worldbank.org/indicator/SE.SEC.DURS
Choose "Download -> Excel" at the end of the menu on the right-hand side. Save file in `Input/`
directory. Rename it `WB_duration_raw.xlsx`.

Datafile: `Input/WB_duration_raw.xlsx` (provided)

#### Penn World Table

Finally, the national accounting variables are collected from version (9.0) of the
Penn World Table (PWT) (Feenstra, Inklaar, and Timmer, 2015). 


Data can be downloaded from https://www.rug.nl/ggdc/docs/pwt90.dta. 
Save the file in the directory `Input/`. 

Datafile: `Input/pwt90.dta` (not provided)

#### Bils and Klenow (2000)

The coefficients of the Mincerian equations are the cross-country averages
as reported and used by Bils and Klenow's (2000) Table B1. The values equal 9.9 per cent for
the return to schooling, and 0.0512 and -0.00071 per cent for the experience-earnings
profile level and squared, respectively. These are typed directly into the equations of the code.

##### Montenegro and Patrinos (2014)

The Mincerian coefficients used to calculate skill premia come from Montenegro and Patrinos (2014). The paper's Annex Table 1 is converted to Excel and available in the directory `Input/`. Please cite the original data providers when using this data. 

Datafile: `Input/MP_mincer.xls` (provided)

##### Elasticity of substitution

The values for the parameter representing the elasticity of substitution between unskilled and
skilled labour come from Acemoglu and Autor (2011) and Ciccone and Peri (2005). The range between 1 and 2 used in the sensitivity analysis derives from the former while the standard value of 1.4 is obtained in the latter.


### Dataset list

| Data file | Source | Notes    |Provided |
|-----------|--------|----------|---------|
| `Input/LeeLee_v1.dta` | Lee and Lee (2016) | https://github.com/barrolee/BarroLeeDataSet/blob/master/LeeLee/LeeLee_v1.dta | No |
| `Input/BL2013_M_v2.2.dta` | Barro and Lee (2001) | https://barrolee.github.io/BarroLeeDataSet/BLData/BL2013_M_v2.2.dta  | No |
| `Input/pwt90.dta` | Penn World Table v9.0 | https://www.rug.nl/ggdc/docs/pwt90.dta | No |
| `Input/WB_duration_raw.xlsx` | Unesco/WB | Duration of secondary (SE.SEC.DURS) | Yes |
| `Input/MP_mincer.xls` | Montenegro and Patrinos (2014) | Annex Table 1 | Yes |
| `Output/WB_duration_raw.dta` | Unesco/WB | Cleaned WB data | Yes |
| `Output/MP_mincer.dta` | Montenegro and Patrinos (2014) | Cleaned MP data | Yes |


---------------------------
## Computational requirements


### Software Requirements

- Stata (code was last run with version 14)
  - `sxpose` (as of 2021-09-23)   
  - `sutex` (as of 2021-09-23)   
  - `estout` (as of 2021-09-23)  
  - `parmest` (version to Stata 11 *not any newer*)   
  - the program `DoFiles/000_BKSK_all.do` has the code to install all dependencies locally, if necessary.


### Memory and Runtime Requirements

Approximate time needed to reproduce the analyses on a standard 2021 desktop machine: less than 10 minutes.

The code was last run on a **dual-core (i5-3210M) Intel-based laptop with 6 GB of RAM using Windows 10 Pro version 20H2**. 


----------------------------
## Description of code

The program `DoFiles/000_BKSK_all.do` will run all the programs in `DoFiles/000_BKSK_all`. It also presents the installation of the required packages.

### License for Code

The code is licensed under a GNU General Public License. See [LICENSE.txt](LICENSE.txt) for details.

---------------------------
## Instructions to Replicators

- Download the data files referenced above. Each should be stored in the subdirectory `Input/`, in the format that you download them in. No further action is needed on the replicator's part.
- Edit `DoFiles/000_BKSK_all.do` to adjust the default path
- Check required packages and install whenever necessary using instructions in the first block of `DoFiles/000_BKSK_all.do`
- Run `DoFiles/000_BKSK_all.do` to produce all the dataset, analysis and exhibits of the paper in sequence. Data will be saved in the subdirectory `Output/` while the figures and tables will be saved in the subdirectory `Exhibits/`


---------------------------

## List of tables and programs


The provided code reproduces:

- [X] All numbers provided in text in the paper
- [X] All tables and figures in the paper


| Table #    | Program                  | Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|
| Table 1           | DoFiles/8_BKSK_exhibits.do    |  7           | Tab1_descriptive.tex                 || As indicated in the code
| Table 2           | DoFiles/8_BKSK_exhibits.do |     87      | Tab12_BKBJ.tex                       || As indicated in the code
| Table 3           | DoFiles/8_BKSK_exhibits.do |     99     | Tab13_BKACES.tex                      || As indicated in the code


---
## References

ACEMOGLU, DARON, AND DAVID AUTOR (2011): "Skills, Tasks and Technologies: Implications for Employment and Earnings," vol. 4 of *Handbook of Labor Economics*,
chap. 12, pp. 1043–1171. Elsevier.

BARRO, ROBERT J., AND JONG WHA LEE (2001): "International Data on Educational Attainment: Updates and Implications," *Oxford Economic Papers*, 53(3), 541–63.

BILS, MARK, AND PETER J. KLENOW (2000): "Does Schooling Cause Growth?," *American Economic Review*, 90(5), 1160–1183.

CICCONE, ANTONIO, AND GIOVANNI PERI (2005): "Long-Run Substitutability Between More and Less Educated Workers: Evidence from U.S. States, 1950-1990," *The Review
of Economics and Statistics*, 87(4), 652–663.

FEENSTRA, ROBERT, ROBERT INKLAAR, AND MARCEL TIMMER (2015): "The Next Generation of the PennWorld Table," *American Economic Review*, 105(10), 3150–82.

LEE, JONG-WHA, AND HANOL LEE (2016): "Human capital in the long run," *Journal of Development Economics*, 122(C), 147–169.

MONTENEGRO, CLAUDIO E., AND HARRY ANTHONY PATRINOS (2014): "Comparable estimates of returns to schooling around the world," *Policy Research Working Paper
Series* 7020, The World Bank.

THE WORLD BANK (2017): "World Development Indicators - Education: Outcomes," https://data.worldbank.org/indicator/SE.SEC.DURS, Accessed 19 Jul 2018.



---

## Acknowledgements

This README complies to the structure and content that have been suggested by [Endorsers](https://github.com/social-science-data-editors/template_README/blob/development/Endorsers.md) of several journals of Social Sciences. Their template is available [here](https://github.com/social-science-data-editors/template_README/blob/master/template-README.md).
