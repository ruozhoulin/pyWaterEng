def ppm(values):
    n = len(values)
    rankings = [sorted(values).index(x) for x in values]
    p = [(rank-0.4)/(n+0.2)
         for rank in rankings]
    return p

