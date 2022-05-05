# %%
from hydrology import FrequencyAnalysis as fqa


sr = [0, 15, 191, 54, 37, 28, 21, 15, 10, 6, 3, 0]
excessRainfall = [20, 40, 0]
hg = hdg.Hydrograph(time, sr)
# uh = UH.UnitHydrograph(time, sr, duration=2, timestep=2)
hg.plot()
# %%
uh = hdg.UnitHydrograph(time, sr, duration=6, timestep=6)
uh.plot()

# %%
uh.area()
# %%
