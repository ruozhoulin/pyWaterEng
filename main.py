# %%
from hydrology import Hydrograph as hdg


time = list(range(0, 14+1, 2))
sr = [0, 12, 42, 70, 51, 22, 10, 0]
excessRainfall = [20, 40, 0]

uh = hdg.UnitHydrograph(time, sr, duration=2, timestep=2)
# uh.baseflow = 1
uh.surface_runoff(excessRainfall)
# %%
uh4 = uh.duration_increase(newDuration=4)
uh4.plot()
uh4
# %%
sCurve = uh4.s_curve(maxTime=25)
sCurve.plot()
sCurve
# %%
shiftDuration = 2
sCurve2 = uh4.shift_curve(sCurve, shiftDuration)
# sCurve.plot(), sCurve2.plot()
hdg.Hydrograph.plot_multi([sCurve, sCurve2], name='S-Curve')

# %%
uh4.duration_decrease(newDuration=2, maxTime=25)

# %%
uh.area()
# %%
excessRainfall = [0, 0, 5.1]  # delete all appending zeros
flow = [0, 0, 25, 40, 60, 35, 0, 0, 0]
values = hdg.estimate_unit_hydrograph(excessRainfall, flow)
values
# %%
times = [0, 1, 2, 3, 4]
data = values + [0, 0, 0]
uh = hdg.UnitHydrograph(times, data, duration=1, timestep=1)
excessRainfall = [5.1]
uh.surface_runoff(excessRainfall)
# %%
