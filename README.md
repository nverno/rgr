Fitting residual growth rates
=============================

Residual growth rate (rGR) = observed growth / predicted growth

rGR is calculated differently for use in the growth and hazard models.  In the growth
models, predicted growth is a upper quantile model (i.e. the predictive model is
fit to the top performers in size classes).  For hazards, the primary objective is
to have rGR values that are not significantly correlated with local relative size (LRS).  
Thus, various model (linear, power, segemented linear, polynomial) were tried for each 
plot (separately by time period as well).  

Issues
------

* rGR for hazards is still correlated in two plots

Dependencies
------------

### Scripts

* [functions.R](http://github.com/ghandi9000/functions)

### Repositories

* [data](http://github.com/ghandi9000/data)
* [hazards](http://github.com/ghandi9000/rgr) 
* [sdp](http://github.com/ghandi9000/sdp)

Data
----

* 

Directory Layout
----------------

* data-prep, data-trans contain scripts to clean and perform manipulations such as transformation between wide/long or constructing size classes from continuous variables.
* functions.R stores store user-defined functions
* etc.

References
----------

Papers, codes