# Norovirus-Reporting-Patterns
Contains R code supporting publication titled "Understanding Norovirus Reporting Patterns in England: A Mixed Model Approach".

### Prerequisites

* Basic knowledge of R programming language.
* R or RStudio is installed on your machine.

### How to use

* Clone/download the repository.<br>
		git clone https://github.com/NikolaOndrikova/Norovirus-Reporting-Patterns.git
* Set your working directory from R/RStudio.<br>
		setwd('Location_on your_machine/Norovirus-Reporting-Patterns/')
* Run master.R from R.<br>
		source('./master.R') 

**Note:** The data used here are synthetic, i.e. we fitted a model to the original data and generated a new data set. For this reason, the results differ from those presented in the paper, and you are likely to see some warning messages (e.g. convergence). 


### References
Meyer S, Held L, Höhle M. Spatio-Temporal Analysis of Epidemic Phenomena Using the R Package surveillance. J Stat Softw. 2017;77:1–55. doi:10.18637/jss.v077.i11. <br>
Meyer S, Held L. Incorporating social contact data in spatio-temporal models for infectious disease spread. Biostatistics. 2017;18:338–51.
