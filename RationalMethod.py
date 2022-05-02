# %%
class RationalMethod:
    def __init__(self, c, area) -> None:
        self.c = c  # runoff coefficient
        self.area = area  # (km2)
        self.tc = 0  # time of concertation (min)
        self.L = 0
        self.slope = 0  # unit in (%)
        self.n = 0

    def time_concertation(self, i):
        assert self.L * self.n != 0
        assert self.slope != 0
        tc = 6.94 * (self.L * self.n)**0.6 / i**0.4 * self.slope**0.3
        return tc

    def flow(self, i):
        # mm/hour
        assert self.tc != 0
        q = 1/3.6 * self.c * i * self.area
        return q

    def flow_partial(self, i, t):
        assert t < self.tc
        assert self.tc != 0
        q = self.flow(i) * (t / self.tc)
        return q


# %%
rm1 = RationalMethod(0.6, 0.12)
rm1.tc = 24
rm2 = RationalMethod(0.8, 0.1)
rm2.tc = 12
rm1.flow(80), rm2.flow(80), rm1.flow_partial(90, 20)

# %%
