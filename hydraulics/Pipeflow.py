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


class CircularPipe(Conduit):
    def __init__(self) -> None:
        super().__init__()


class TrapezoidalOpenChannel(Conduit):
    def __init__(self) -> None:
        super().__init__()
