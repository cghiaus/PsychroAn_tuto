# PsychroAn_tuto - Computational psychrometric analysis of HVAC systems: tutorials

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/cghiaus/dm4bem/blob/main/LICENSE)

**Contents**

*Interactive web pages (using Voilà)*
1. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT01_moist_air_prop.ipynb) [T01_moist_air_prop.ipynb](T01_moist_air_prop.ipynb) Moist air properties.
2. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT02_load_winter.ipynb) [T02_loads_winter.ipynb](T02_load_winter.ipynb) Thermal loads: winter.
3. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT03_mix.ipynb) [T03_mix.ipynb](T03_mix.ipynb) Adiabatic mixing and isentalpic condensation.
4. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT04_va_hum.ipynb) [T04_va_hum.ipynb](T04_heat_va_hum.ipynb) Heating and vapor humidification.
5. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT05_ad_hum.ipynb) [T05_ad_hum.ipynb](T05_heat_ad_hum.ipynb) Heating and adiabatic humidification.
6. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/PsychroAn_tuto/HEAD?urlpath=%2Fvoila%2Frender%2FT06_cool.ipynb) [T06_cool](T06_cool.ipynb) Cooling with dehumidification.

**Annex**
- [Design procedure for mass flow rate of supply air](Annex01.ipynb)

## Elementary processes

The models presented in the tutorials are formed by a system of equations constructed with blocks of elementary processes (Table 1).
> Table 1. Models of the elementary processes ([Ghiaus 2022](https://hal.archives-ouvertes.fr/hal-03484064/document)).
> ![](Figures/elementary_processes.svg)

## Bibliography

Ghiaus, C. (2014). Linear algebra solution to psychometric analysis of air-conditioning systems. *Energy, 74*, 555-566. DOI: [10.1016/j.energy.2014.07.021](https://doi.org/10.1016/j.energy.2014.07.021)

Ghiaus, C. (2016). Analyse psychrométrique des systèmes de climatisation. *Revue générale du Froid & du Conditionnement d’air*, pp.38-42. [hal-03379788](https://hal.archives-ouvertes.fr/hal-03379788/document)

Ghiaus, C. (2021). PsychroAn_cool: Psychrometric analysis of cooling systems as a control problem. In Journal of Building Performance Simulation (0.0.0, Vol. 15, Number 1, pp. 21–38). *Zenodo*. DOI: [10.5281/zenodo.5236450](https://doi.org/10.5281/zenodo.5236450)

Ghiaus, C. (2022) Computational psychrometric analysis as a control problem: case of cooling and dehumidification systems, *International Journal of Building Performance Simulation, 15*(1), pp. 21-38, DOI: [10.1080/19401493.2021.1995498](https://doi.org/10.1080/19401493.2021.1995498), (open access preprint [hal-03484064](https://hal.archives-ouvertes.fr/hal-03484064/document))

## Support

- [LaTeX/Mathematics](https://en.wikibooks.org/wiki/LaTeX/Mathematics)

- [LaTex Equations](https://latex.codecogs.com/eqneditor/editor.php)

- [Latex Table generator](https://www.tablesgenerator.com/markdown_tables#)

- [Anaconda cheetsheet](https://docs.continuum.io/anaconda/user-guide/cheatsheet/)

- [Jupyter Notebook cheatsheet](https://medium.com/ibm-data-science-experience/markdown-for-jupyter-notebooks-cheatsheet-386c05aeebed)

- [NumPy for MATLAB users](http://mathesaurus.sourceforge.net/matlab-numpy.html)

## Questions
1. What is a mathematical model of a physical process?
2. What is the difference between *correlation* and *causation*?
3. Define the *inputs* and the *outputs* of a physical system.
4. Define the *inputs* and the *outputs* of a computational model.
5. Based on the definition of the inputs and outputs of a physical system and the inputs and outputs of a computational model, define the direct problem.
6. Describe the inverse problem of model identification.
7. Describe the inverse problem of control.
8. Explain the citation: "Divide each difficulty into as many parts as is feasable and necessary to solve it" by [René Descartes](https://www.goodreads.com/quotes/565817-divide-each-difficulty-into-as-many-parts-as-is-feasible), Discourse on Method, 1637.
9. Explain the citation: "Whenever I run into a problem I can't solve, I always make it bigger. I can never solve it by trying to make it smaller, but if I make it big enough, I can begin to see the outlines of a solution." by [Dwight David Eisenhower](https://engine-for-change.com/quote-of-the-week-make-it-big-enough/) in the context of numerical simulation of physical systems.
10. Explain the type of inputs of numerical models used for simulation: know knowns, unknown knowns, known unknowns, unknown unknowns ([Donald Rumsfeld](https://en.wikipedia.org/wiki/There_are_unknown_unknowns) or [Rumsfeld matrix](https://medium.com/@andreamantovani/known-knowns-known-unknowns-unknown-unknowns-leadership-367f346b0953)).
11. Define a *white box*, a *black box* and a *gray box* model.
12. What are the aims of the indoor climat control?
13. The functions of an Air Handling Unit (AHU).
14. Why the air vector becomes widly used in buildings?
15. Definition of enthalpy. Sensible enthalpy and latent enthalpy
16. Definition of moist air.
17. Saturation pressure of water vapor.
18. Humidity ratio (or absolute humidity); why is it important?
19. Relative humidity; why is it important?
20. Moist air enthalpy.
21. Dew temperature and dew point.
22. Dry bulb and wet bulb temperature.
23. Write the sensible heat balance and the latent heat balance of a thermal zone in steady-state.
24. Write the equations for the sensible and the latent load of a building in steady-state (considering the transfer through the walls, thermal bridges, by infltration and the auxiliary loads). 
25. Write the equations for the sensible heat balance and the latent heat balance and represent the process on a psychrometric chart for:
  - a *thermal zone*;
  - a *mixing process*;
  - a *heating process*;
  - a *dry cooling process*;
  - a *cooling with dehumidification process*;
  - a *vapor humidification process*.
  - an *adiabatic humidification process*.

