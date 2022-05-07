class CrossSection:
    def area(self):
        pass

    def wetted_perimeter(self, depthRatio: float = 1) -> float:
        pass

    def hydraulic_radius(self, depthRatio: float = 1) -> float:
        pass

    def hydraulic_diameter(self, depthRatio: float = 1) -> float:
        # Despite what the name may suggest, the hydraulic diameter
        #  is not twice the hydraulic radius, but four times larger.
        return 4 * self.hydraulic_radius(depthRatio=depthRatio)

    def reynolds_number(self, velocity: float,
                        kinematicViscosity: float = 1e-6):
        """_summary_

        Args:
            velocity (float): m/s
            kinematicViscosity (float, optional): m2/s. Defaults to 1e-6.

        Returns:
            _type_: _description_
        """
        pass


class Circle(CrossSection):
    def __init__(self, pi=3.14159, theta=90,  **kwargs) -> None:
        super().__init__()
        self.pi = pi
        self.theta = theta  # default is vertical
        self.r = kwargs.pop('radius', 0)  # radius
        self.d = kwargs.pop('diameter', 0)  # diameter
        if self.r != 0 and self.d != 0:
            assert self.r*2 == self.d
        elif self.r != 0:
            self.d = self.r*2
        elif self.d != 0:
            self.r = self.d/2
        else:
            raise ValueError('Please input a positve radius or diameter.')

    def area(self, depthRatio=1):
        # ! not finished
        return self.r**2 * self.pi

    def wetted_perimeter(self, depthRatio=1):
        # D*pi when pipe is full.
        # ! not finished
        return self.r*2*self.pi

    def hydraulic_radius(self, depthRatio=1):
        # ! not finished
        # D/4 when pipe is full.
        return self.area()/self.wetted_perimeter()

    def reynolds_number(self, velocity: float,
                        kinematicViscosity: float = 1e-6):
        """_summary_

        Args:
            velocity (float): m/s
            kinematicViscosity (float, optional): m2/s. Defaults to 1e-6.

        Returns:
            _type_: _description_
        """
        Re = velocity * self.d / kinematicViscosity
        return Re

    def Ixx(self):
        return self.pi * self.r**4 / 4

    def yD(self, yc):
        # yc != hc
        # location of center of pressure
        return yc + self.Ixx() / (yc * self.area())


class Conduit:
    def __init__(self) -> None:
        self.length = 0
        self.slope = 0
        self.area = 0

    def compute_discharge(self):
        pass

    def compute_velocity(self):
        pass

    def compute_primieter(self):
        pass

    def compute_area(self):
        pass


class CircularPipe(Conduit, Circle):
    def __init__(self) -> None:
        super().__init__()


class TrapezoidalOpenChannel(Conduit):
    def __init__(self) -> None:
        super().__init__()


class BernoulliEquation:
    def __init__(self, g=9.8, rou=1e3, gamma=0) -> None:
        self.g = g  # gravity, N/kg
        self.rou = rou  # density, kg/m3
        if gamma == 0:
            self.gamma = rou * g
        else:
            self.gamma = gamma

    def compute(self, z1, u1, p1, z2, u2, p2,
                hf=0.0, hl=0.0, hp=0.0, ht=0.0):
        """_summary_

        Args:
            z1 (_type_): _description_
            u1 (_type_): _description_
            p1 (_type_): _description_
            z2 (_type_): _description_
            u2 (_type_): _description_
            p2 (_type_): _description_
            hf (float, optional): friction head loss. Defaults to 0.0.
            hl (float, optional): local head loss. Defaults to 0.0.
            hp (float, optional): pump head. Defaults to 0.0.
            ht (float, optional): turbine head. Defaults to 0.0.

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        # z1, u1, p1 is known variables
        assert len([i for i in [z1, u1, p1] if i == None]) == 0
        # There is one and only one unknown variable
        assert len([i for i in [z2, u2, p2, hf, hl, hp, ht] if i == None]) == 1
        rhs = self.total_head(z1, u1, p1)

        if z2 == None:
            htotal = hf+hl+ht-hp
            unknown = rhs - \
                self.velocity_head(u2) - self.pressure_head(p2) - htotal
        elif p2 == None:
            htotal = hf+hl+ht-hp
            unknown = (rhs - z2 - self.velocity_head(u2)-htotal)*self.gamma
        elif u2 == None:
            htotal = hf+hl+ht-hp
            unknown = ((rhs-z2-self.pressure_head(p2)-htotal)*2*self.g)**0.5
        elif hl == None:
            unknown = rhs - self.total_head(z2, u2, p2) - hf - ht + hp
        else:
            raise ValueError()
        return unknown

    def total_head(self, z, u, p):
        rhs = z + self.velocity_head(u) + self.pressure_head(p)
        return rhs

    def pressure_head(self, p):
        return p / self.gamma

    def velocity_head(self, u):
        return u**2 / (2*self.g)

    def darcy_weisbach(self, l, hl, f, d, u):
        """compute head loss due to roughness for circular pipe only

        Args:
            f (_type_): f is checked from the Moody diagram
            d (_type_): diameter (m)
            u (_type_): Velocity (m/s)

        Returns:
            _type_: _description_
        """
        # assert self.l > 0, 'Please assign length of the pipe.'
        # There is one and only one unknown variable
        assert len([i for i in [hl, f, d, u] if i == None]) == 1
        # Darcy-Weisbach equation
        if hl == None:
            unknown = f * l/d * u**2 / (2*self.g)
        elif f == None:
            unknown = hl / (l/d * u**2 / (2*self.g))
        else:
            raise ValueError()
        return unknown

    def local_head_loss(self, k, u):
        return k * u**2 / (2 * self.g)
