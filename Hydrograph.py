from cProfile import label
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt


class Hydrograph(OrderedDict):

    def __init__(self, times: list, values: list,
                 name='Hydrograph',
                 unitTime='Unknown', unitValue='Unknown'):
        super().__init__(zip(times, values))
        self.name = name
        self.unitTime = unitTime
        self.unitValue = unitValue

    @property
    def times(self):
        return list(self.keys())

    @property
    def data(self):
        return list(self.values())

    def __repr__(self):
        splitLine = '========================\n'
        header = '[%s]\nTime\tValue\n' % (self.name)
        data = ['{}\t{}\n'.format(i, v) for i, v in self.items()]
        data = ''.join(data)
        return splitLine+header+data+splitLine

    def copy(self):
        # copy argument
        newOne = self.__class__(self.times, self.data)
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

    @staticmethod
    def plot_multi(hydrographs: list, name: str):
        for hydrograph in hydrographs:
            x = hydrograph.times
            y = hydrograph.data
            plt.plot(x, y, label=hydrograph.name)
            plt.plot(x, y, 'o')
        plt.title(name)
        plt.xlabel('Time (%s)' % (hydrograph.unitTime))
        plt.ylabel('Discharge (%s)' % (hydrograph.unitValue))
        plt.tight_layout()
        plt.legend()


class UnitHydrograph(Hydrograph):
    def __init__(self, times: list, values: list, duration: int,
                 timestep: int,
                 baseflow=0) -> None:
        super().__init__(times, values, unitTime='hour', unitValue='m3/s')
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

    def extend_times(self, maxTime):
        times = list(self.times)
        assert maxTime > times[-1]
        newTimes = range(times[0], maxTime, self.timestep)
        return newTimes

    def extend_times_by_index(self, NofIndex):
        times = list(self.times)
        assert NofIndex > len(times)
        newTimes = [times[0] + i*self.timestep for i in range(NofIndex)]
        return newTimes

    def new_hydrograph(self, values: list, name) -> Hydrograph:
        # it will automatically extend time based on value
        times = self.times
        if len(values) > len(times):
            times = self.extend_times_by_index(len(values))
        hydrograph = Hydrograph(times, values,
                                name=name)
        return hydrograph

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
