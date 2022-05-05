from math import exp


def Horton_method(t, f0, fc, k):
    """_summary_

    Args:
        t (_type_): time
        f0 (_type_): Starting infiltration rate
        fc (_type_): Final constant infiltration rate
        k (_type_): Decay constant with dimension T^-1

    Returns:
        _type_: Infiltration rate at time t
    """
    ft = fc + (f0-fc)*exp(-k*t)
    return ft


def Green_Ampt_method(t):
    pass
