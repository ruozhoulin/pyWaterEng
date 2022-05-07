# %%
from hydrology import FrequencyAnalysis as fqa
from hydraulics.Pipeflow import Circle
from hydraulics.Pipeflow import BernoulliEquation
# %%
Re = Circle(diameter=0.75).reynolds_number(3.0, kinematicViscosity=0.856e-6)
eq = BernoulliEquation(l=300)
hl = eq.darcy_weisbach(hl=None, f=0.015, d=0.075, u=3)
hl
