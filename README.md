### Blanquart.2016
This is code used for a modeling sub-analysis in Francois Blanquart's analysis of set point viral load in Rakai, Uganda.

### Brief explanation
This is a stochastic, agent-based simulation of HIV epidemiology and evolution. To run, source the "Run.new.Oct2016.R" from within an R session. This R wrapper calls the .C file. To visualize key output, source the "Plots.R" script. The ".. .alpha ... .c" files incorporate different values of an evolutionary parameter that can be edited on line 810 of each file. Make sure the .R wrapper script calls the correct .c source (line 163 of the .R script).

### References
This model has been used for the following studies, and is described in more detail within each manuscript:

An HIV Epidemic Model Based on Viral Load Dynamics: Value in Assessing Empirical Trends in HIV Virulence and Community Viral Load
Joshua T. Herbeck, John E. Mittler, Geoffrey S. Gottlieb, James I. Mullins

Published: 19 June 2014
http://dx.doi.org/10.1371/journal.pcbi.1003673
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003673

Evolution of HIV virulence in response to widespread scale up of antiretroviral therapy: a modeling study
Joshua T. Herbeck, John E. Mittler, Geoffrey S. Gottlieb, Steven M. Goodreau, James T. Murphy, Anne Cori, Michael Pickles, Christophe Fraser

Published: 3 October 2016
http://dx.doi.org/10.1093/ve/vew028 
http://ve.oxfordjournals.org/content/2/2/vew028

### Contact
Questions can be sent to Josh Herbeck, herbeck@uw.edu.
