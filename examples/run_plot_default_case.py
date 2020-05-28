# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:py3d]
#     language: python
#     name: conda-env-py3d-py
# ---

# %%
#pylint: disable=unnecessary-semicolon

# %% [markdown]
# # Test functionality for the default case

# %%
import matplotlib.pyplot as plt
#import numpy as np

# %%
# %matplotlib notebook

# %% [markdown]
# To allow this nb to work with binder, we can pull pyspctral2 from GitHub manually, as a workaround until Python packaging has been done.
#
# Locally, we want to to be sure to run this nb in the `examples` subdirectory of the repo, i.e., from its location. 

# %% language="bash"
# #rm -rf pyspctral2
# git clone https://github.com/zmoon92/pyspctral2.git
# #rm -rf pyspctral2/.git  # don't need the Git history

# %% [markdown]
# Because there is an `__init__.py` in the root level of the repo, we can now import pyspctral2 with a normal import command.

# %%
import pyspctral2 as psp2

# %% [markdown]
# Examining `model`'s docstring:

# %%
print(psp2.model.__doc__.strip())

# %% [markdown]
# ## Examine case params and then run

# %%
m = psp2.model()  # no args => default settings (based on those found in the original C run program)
#m = psp2.model(tau500=0.4, watvap=2.0, ozone=-1.0)
print('lat, lon:', m.lat, m.lon)
print('datetime:', m.dt)
print('temp, press:', m.temp, m.press)
print('tau500, watvap, ozone:', m.tau500, m.watvap, m.ozone)

# %%
m.run()

# %% [markdown]
# ## Test provided plotting routine

# %%
psp2.plot_spectrum();  # plots and returns handle so `;` to silence returned one in the Notebook

# %% [markdown]
# ## Test spectrum correction by ground-based measurement

# %%
measured_val = 350  # W/m^2
measured_units = 'E'
measured_bnds = (0.4, 0.7)  # PAR/Vis band (um)

m.correct(measured_units=measured_units, 
          measured_val=measured_val, 
          measured_bnds=measured_bnds)

psp2.plot_spectrum(output_type='corrected');
