# Github Repositorium to Reproduce the Assessment of the Pauling rules

This is the github repositorium for our publication that assesses the Pauling rules on around 5000 experimental oxides from the [Materials Project](http://materialsproject.org/) and around 5000 experimental oxides from the [Crystal Open Database](http://www.crystallography.net/cod/). This code has been written by J. George with contributions and ideas from D. Waroquiers, D. Di Stefano, G.-M. Rignansese, and G. Hautier.

One can run the overall analysis by running the scripts "Final_Run.py", "Final_Run_exp.py", and "Final_Run_binaries.py" in the folder "Assessment". Also, make sure that the classes used in this script are accessible (export the Python path). They are in the files "PaulingRules.py" and "Classes_for_statistics.py". The requirements.txt includes the Python libraries needed to run the scripts: one will need [pymatgen](http://pymatgen.org/) and some other standard Python packages. This will only work with Python 3.x. To make sure that everything is working correctly, please run the tests in "tests". This might take about 1 h.

Based on this analysis, we will subsequently publish a package that will allow you to analyse the Pauling rules for your own structures. Please have a look at this repositorium for future information.
