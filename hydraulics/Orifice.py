def orifice(Cd: float, area: float, h: float, g: float = 9.8):
    """_summary_

    Args:
        Cd (float): coefficent of discharge
        area (float): area of orifice (m2)
        h (float): head (m)
        g (float, optional): gravity acceleration (m/s2)

    Returns:
        _type_: flow rate (CMS)
    """
    Q = Cd * area * (2*g*h)**0.5
    return Q
