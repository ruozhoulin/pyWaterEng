from typing import Callable
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt


class Hydrograph(OrderedDict):

    def __init__(self, times: list, values: list, timestep,
                 name='Hydrograph',
                 unitTime='Unknown', unitValue='Unknown'):
        super().__init__(zip(times, values))
        self.name = name
        self.unitTime = unitTime
        self.unitValue = unitValue
        self.timestep = timestep

    @property
    def times(self):
        return list(self.keys())

    @property
    def data(self):
        return list(self.values())

    def extend_time(self, maxTime):
        # extend remainder by zero.
        assert maxTime > self.times[-1]
        newTimes = range(self.times[-1]+self.timestep,
                         maxTime+self.timestep, self.timestep)
        for t in newTimes:
            self[t] = 0

    def new_hydrograph(self, values: list, name):
        # it will automatically extend time based on value
        hydrograph = self.copy()
        hydrograph.name = name
        if len(values) > len(self.times):
            maxTime = self.times[-1] + len(values) - len(self.times)
            hydrograph.extend_time(maxTime)
        # assign values
        for i, t in enumerate(hydrograph.times):
            hydrograph[t] = values[i]
        return hydrograph

    def __repr__(self):
        splitLine = '========================'
        header = '[%s]\nTime\tValue\n' % (self.name)
        data = ['{}\t{}\n'.format(i, v) for i, v in self.items()]
        data = ''.join(data)
        return splitLine+'\n'+header+data+splitLine

    def copy(self):
        # copy argument
        newOne = Hydrograph(self.times, self.data,
                            self.timestep)
        # copy attribute
        newOne.__dict__.update(self.__dict__)
        return newOne

    def plot(self):
        x = self.times
        y = self.data
        plt.plot(x, y, 'k')
        plt.plot(x, y, 'bo')
        plt.title(self.name)
        plt.xlabel('Time (%s)' % (self.unitTime))
        plt.ylabel('Discharge (%s)' % (self.unitValue))
        plt.tight_layout()


def plot_Hydrographs(hydrographs: list[Hydrograph], name: str, alpha=0.7):
    colors = ['b', 'g', 'r', 'k']
    linestyles = ['-', '--', '-.', ':']
    for i, hydrograph in enumerate(hydrographs):
        x = hydrograph.times
        y = hydrograph.data
        color = colors[i % len(colors)]
        linestyle = linestyles[i % len(linestyles)]
        param = {'color': color,
                 'alpha': alpha,
                 'linestyle': linestyle,
                 'marker': 'o'}
        plt.plot(x, y, label=hydrograph.name, **param)
    plt.title(name)
    plt.xlabel('Time (%s)' % (hydrographs[0].unitTime))
    plt.ylabel('Discharge (%s)' % (hydrographs[0].unitValue))
    plt.tight_layout()
    plt.legend()


class UnitHydrograph(Hydrograph):
    def __init__(self, times: list, values: list, duration: int,
                 timestep: int,
                 baseflow=0) -> None:
        super().__init__(times, values, timestep, unitTime='hour', unitValue='m3/s')
        self.name = 'Unit Hydrograph'
        self.duration = duration
        self.timestep = timestep
        self.baseflow = baseflow
        # test time step
        for i in range(len(times)-1):
            assert times[i] + timestep == times[i+1]

    def surface_runoff(self, excessRainfall) -> Hydrograph:
        NofResult = self.minimum_result_number(len(excessRainfall))
        # column vector
        U = np.matrix(list(self.data)).T
        # create matrix of q
        nrows = NofResult
        ncols = len(self.data)
        P = UnitHydrograph.rainfall_matrix(excessRainfall, nrows, ncols)
        # compute surface runoff
        surfaceRunoff = P * U + self.baseflow
        surfaceRunoff_list = surfaceRunoff.flatten().tolist()[0]
        return self.new_hydrograph(surfaceRunoff_list,
                                   name='Surface Runoff')

    def duration_increase(self, newDuration):
        rate = newDuration / self.duration  # > 1.0
        rain = [1/rate] * int(rate)
        hdg = self.surface_runoff(rain)  # hydrograph
        unitHydrograph = UnitHydrograph(
            hdg.times, hdg.data, duration=newDuration,
            timestep=self.timestep,
            baseflow=self.baseflow)
        return unitHydrograph

    def duration_decrease(self, newDuration, maxTime=0):
        """_summary_

        Args:
            duration (int): new duration that is smaller than original duration
            maxTime (int, optional): _description_. Defaults to 0.

        Returns:
            _type_: _description_
        """
        sh1 = self.s_curve(maxTime=maxTime)
        sh2 = self.shift_curve(sh1, newDuration)
        times = sh1.times
        values = [(sh1[i] - sh2[i]) / newDuration for i in times]
        unitHydrograph = UnitHydrograph(
            times, values, duration=newDuration,
            timestep=self.timestep,
            baseflow=self.baseflow)
        return unitHydrograph

    def s_curve(self, maxTime=0) -> Hydrograph:
        values = []
        if maxTime == 0:
            times = self.times
        else:
            time = list(self.times)
            times = list(range(time[0], maxTime, self.timestep))

        for t in times:
            # first UH
            value = self.get(t, 0)
            # other UH
            while t >= self.duration:
                t -= self.duration  # next UH
                value += self.get(t, 0)
            values.append(value)
        # multiply duration (h) to ...
        values = [value * self.duration for value in values]
        sHydrograph = self.new_hydrograph(values, 'S-Curve')
        return sHydrograph

    def shift_curve(self, sCurve: Hydrograph, shiftDuration: int) -> Hydrograph:
        # hour
        hdg = sCurve.copy()  # the new hydrograph after shifting
        tmin, tmax = hdg.times[0], hdg.times[-1]
        # shift: start at the last record
        for i in range(tmax, tmin+shiftDuration-self.timestep, -self.timestep):
            hdg[i+self.timestep] = hdg[i]
        # shift: fill the starting point with zero
        for i in range(tmin, tmin+shiftDuration+self.timestep, self.timestep):
            hdg[i] = 0
        return hdg

    def area(self):
        """_summary_
        m^3/s -> km2
        Returns:
            _type_: km2
        """
        # km^2 * mm = m^3/s * h
        area = 0.36 * sum(self.data) * self.timestep
        return area

    def check_result_number(self, NexcessRainfall, NofResult):
        assert NofResult >= self.minimum_result_number(
            NexcessRainfall), 'Need larger number of results'

    def minimum_result_number(self, NexcessRainfall):
        # Number of items in output surface runoff (NofResult) should
        #  be large enought, this won't impact results.
        Nuh = len(self.times)
        return Nuh + NexcessRainfall - 1

    @staticmethod
    def rainfall_matrix(excessRainfall: list, nrows: int, ncols: int):

        p = np.zeros((nrows, ncols))
        for i in range(ncols):
            for j, rain in enumerate(excessRainfall):
                p[i+j, i] = rain
        return p


def estimate_unit_hydrograph(excessRainfall: list, flow: list) -> list:
    """Estimate the unit hydrograph from **observed flow data** and corresponding **excess rainfall data**.

    Args:
        excessRainfall (list): _description_
        flow (list): _description_

    Returns:
        list: _description_
    """
    # column vector
    Q = np.matrix(flow).T
    # create matrix of q
    nrows = len(flow)
    ncols = nrows - len(excessRainfall) + 1
    p = UnitHydrograph.rainfall_matrix(excessRainfall, nrows, ncols)
    # compute unit hydrograph
    u = np.linalg.inv(p.T.dot(p)).dot(p.T).dot(Q)
    return u.flatten().tolist()[0]


class StorageDischarge:
    def __init__(self, area, excessRainfall: list, baseflow=0, timestep=1) -> None:
        self.area = area  # catchment area (km2)
        self.timestep = timestep  # time step of computation (hour)
        self.unitTime = 'hour'
        times = list(range(0, len(excessRainfall), self.timestep))
        self.excessRainfall = Hydrograph(times, excessRainfall, timestep=1,
                                         name='Excess Rainfall',
                                         unitTime=self.unitTime, unitValue='mm')
        self.baseflow = baseflow

    def compute_storage(self):
        pass

    def excess_rainfall_in_flow(self):
        # from mm to m^3/s
        unit = 1e-3 * self.area * 1e6 / (self.timestep * 60 * 60)
        newRain = [rain * unit for rain in self.excessRainfall.data]
        excessRainfall = self.excessRainfall.new_hydrograph(
            newRain, name='Excess Rainfall (m^3/s)')
        return excessRainfall


class LinearResevoirModel(StorageDischarge):
    def __init__(self, k, area, excessRainfall: list, baseflow=0, timestep=1) -> None:
        super().__init__(area, excessRainfall, baseflow, timestep)
        # Storge = k * discharge
        self.k = k

    def compute_storage(self, flowrate: float):
        """_summary_

        Args:
            flowrate (float): flow rate in m3/s

        Returns:
            _type_: storage in m3
        """
        return self.k * flowrate

    def direct_runoff(self, maxTime):
        excessRainfall = self.excess_rainfall_in_flow()
        excessRainfall.extend_time(maxTime)
        # at the beginning, everything is zero.
        assert excessRainfall[0] == 0
        storages = [0]
        flowrates = [0]
        tDelta = self.timestep*60*60  # hour to second
        for i in range(1, maxTime):  # skip the first one
            i1 = excessRainfall[i-1]
            i2 = excessRainfall[i]
            s = storages[i-1]
            q = flowrates[i-1]

            qNew = (s+((i1+i2)-q)*tDelta/2) / (self.k+tDelta/2)
            sNew = self.compute_storage(qNew)
            flowrates.append(qNew)
            storages.append(sNew)
        times = list(range(len(excessRainfall)))
        storages = Hydrograph(times, storages, timestep=self.timestep,
                              name='Storage',
                              unitTime=self.unitTime, unitValue='m^3')
        flowrates = Hydrograph(times, flowrates, timestep=self.timestep,
                               name='Runoff',
                               unitTime=self.unitTime, unitValue='m3/s')
        return flowrates, storages

    def runoff(self, maxTime, digits=2):
        runoff, storages = self.direct_runoff(maxTime=maxTime)
        for time, value in runoff.items():
            runoff[time] = round(value + self.baseflow, digits)
        return runoff


class NonLinearResevoirModel(StorageDischarge):
    pass


class RationalMethod:
    def __init__(self, c, area) -> None:
        self.c = c  # runoff coefficient
        self.area = area  # (km2)
        self.tc = 0.0  # time of concertation (min)
        self.L = 0.0
        self.slope = 0.0  # unit in (%)
        self.n = 0.0

    def time_concertation(self, intensity):
        assert self.L * self.n != 0
        assert self.slope != 0
        tc = 6.94 * (self.L * self.n)**0.6 / (intensity**0.4 * self.slope**0.3)
        return tc

    def critical_time(self, rainfallFunc: Callable[[float], float]):
        # compute when t==tc
        t = 10.0  # initial guess
        while True:
            intensity = rainfallFunc(t)
            tc = self.time_concertation(intensity=intensity)
            if abs(tc-t) < 1e-3:
                break
            t = tc
        return tc

    def flow(self, intensity: float):
        """_summary_

        Args:
            intensity (float): Rainfall intensity

        Returns:
            float: discharge (m3/s)
        """
        # rainfall intensity
        # mm/hour
        assert self.tc != 0
        q = 1/3.6 * self.c * intensity * self.area
        return q

    def flow_partial(self, intensity, t):
        if t < self.tc:
            rate = (t / self.tc)
        else:
            rate = 1
        assert self.tc != 0
        q = self.flow(intensity) * rate
        return q


def plot_flows(times, flows: list[list], alpha=0.7):
    colors = ['b', 'g', 'r', 'k']
    linestyles = ['-', '--', '-.', ':']
    for i, flow in enumerate(flows):
        color = colors[i % len(colors)]
        linestyle = linestyles[i % len(linestyles)]
        param = {'color': color,
                 'alpha': alpha,
                 'linestyle': linestyle,
                 'marker': 'o'}
        plt.plot(times, flow, label='flow %d' % (i))
    flowSum = [sum(items) for items in zip(*flows)]
    plt.plot(times, flowSum, label='flow sum')

    plt.title('Hahah')
    plt.xlabel('Time (min)')
    plt.ylabel('Discharge (m3/s)')
    plt.tight_layout()
    plt.legend()

    peaks = [max(flow) for flow in flows] + [max(flowSum)]
    return peaks
