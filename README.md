# Air Quality Analytic Collaborative Framework (AQACF) Simulations of Air Pollution Health Effects


aq_health.py is the primary Python script. It was developed by Sudip Chakraborty based on the original Matlab script by Omar Nawaz.


The jupyter notebook serves as a tutorial on how to create an interactive app using adjoint sensitivites of three pollutants (PM2.5, O3, and NO2) to their precursor species. We combine these sensitivities with emissions data from the US EPA's NEI for base year 2011. 



To run this notebook you'll need to install jupyter notebooks: https://jupyter.org/install and then enter the "jupyter notebook la_aq_health_v1.ipynb" command into terminal.



To run this notebook as a standalone app you'll need to install voila: https://voila.readthedocs.io/en/stable/install.html and then enter the "voila la_aq_health_v1.ipynb" command into terminal. You can remove the markdown cells in the notebook to only display the app.



You'll also need the following python libraries to run this notebook:

netCDF4: https://unidata.github.io/netcdf4-python/

ipywidgets: https://ipywidgets.readthedocs.io/en/stable/user_install.html

xarray: https://docs.xarray.dev/en/stable/getting-started-guide/installing.html

numpy: https://numpy.org/install/

matplotlib: https://matplotlib.org/stable/users/installing/index.html




Files used in this notebook are accessible at: http://adjoint.colorado.edu/~jplaqacf/v1/

For further questions about the notebook, email: muna9068@colorado.edu

Original notebook is from: https://github.com/omarnawaz/la-aq-health/
