import pandas as pd
from scipy.stats import pearson3
from scipy.optimize import minimize, basinhopping


class FrequencyAnalysis:
    def __init__(self, years, discharge, type='China') -> None:
        assert len(years) == len(discharge)
        self.table = pd.DataFrame({'year': years,
                                   'discharge': discharge})
        self.compute_probability(type)

    def compute_probability(self, type='China'):
        values = self.table['discharge']
        n = len(values)
        # start from 0, while the rank in formula start from 1.
        valuesSorted = sorted(values, reverse=True)
        rankings = [valuesSorted.index(x)+1 for x in values]
        if type == 'China':
            # ! not finished.
            probabilities = [(rank)/(n+1)
                             for rank in rankings]
        elif type == 'Australia':
            probabilities = [(rank-0.4)/(n+0.2)
                             for rank in rankings]
        else:
            raise ValueError('Wrong type.')
        self.table['rank'] = rankings
        self.table['P'] = probabilities

    def fit(self):
        x0 = [0.1, 1, 1]  # initail value
        res = minimize(self.fitted_function, x0)
        skew, loc, scale = res.x
        print(self.fitted_function([skew, loc, scale]))
        self.skew = skew
        self.loc = loc
        self.scale = scale
        return res
#     def fitted_function(self, params):
#         skew, loc, scale = params
#         probabilityFitted = [pearson3.sf(q, skew, loc=loc, scale=scale)
#                              for q in self.table['discharge']]
#         error = sum([abs(p1-p2)**2 for p1, p2 in
#                      zip(self.table['P'], probabilityFitted)])
#         return error

    def fitted_function(self, params):
        skew, loc, scale = params
        discharges = [pearson3.ppf(probbility, skew, loc=loc, scale=scale)
                      for probbility in self.table['P']]
        error = sum([abs(p1-p2)**2 for p1, p2 in
                     zip(self.table['discharge'], discharges)])/len(discharges)
        error += sum([1e9 for q in discharges if q <= 0])
        return error

    def estimate_discharge(self, probability):
        #assert hasattr(self, 'param'), 'Call `fit()` before this function.'
        discharge = pearson3.ppf(
            probability, self.skew, loc=self.loc, scale=self.scale)
        return discharge
