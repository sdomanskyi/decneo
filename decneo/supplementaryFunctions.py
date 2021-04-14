''' Supplementary functions
'''

from .general import *
from .genes import *

import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects import r as R, IntVector, FloatVector, DataFrame

try:
    utils = importr('gravity')
except Exception as exception:
    print(exception)
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1)
    utils.install_packages('gravity')
finally:
    utils = importr('gravity')

def getNonParametricPValue(labels, values, random_seed = 0, printResults = True):

    '''Markers localization p-value calculation:
    Poisson pseudo-maximum likelihood estimation (PPML) by J.M.C. Santos Silva & Silvana Tenreyro, 2006.
    Implemented in R in "gravity: Estimation Methods for Gravity Models" at:
    https://rdrr.io/cran/gravity/man/ppml.html
    '''

    np.random.seed(random_seed)
    #np.random.shuffle(labels)

    dataf = DataFrame({'label': IntVector(tuple(labels)), 'distance': FloatVector(tuple(values))})
    fit = R('function(x) ppml(dependent_variable="label", distance="distance", additional_regressors=NULL, robust=TRUE, data=x)')(dataf)

    # Deviance is -2.*log_likelihood
    altDeviance = list(fit[9].items())[0][1]
    nullDeviance = list(fit[11].items())[0][1]
    p_value = scipy.stats.chi2.sf(nullDeviance - altDeviance, 1)

    if printResults:
        print('Non-parametric method:', '\n\t',
              #'  Robust PPML based (a.k.a. QMLE) deviances and their difference test (chi2 p-value):\n\t', 
              #'Null deviance:\t', np.round(nullDeviance, 1), '\n\t',
              #'Alternative deviance:\t', np.round(altDeviance, 1), '\n\t',
              'p-value:\t', '%.1e' % p_value, '\n')

    return p_value

def getLogisticRegressionPValue(labels, values, random_seed = 0, printResults = True):

    '''Logistic regression based Deviances and their difference test (chi2 p-value)'''

    data = values[:, None]

    clf = LogisticRegression(random_state=random_seed).fit(data, labels)

    clf.fit(data, labels)

    prprob = clf.predict_proba(data)

    alt_log_likelihood = -log_loss(labels, prprob, normalize=False)
    null_prob = sum(labels) / float(labels.shape[0]) * np.ones(labels.shape)
    null_log_likelihood = -log_loss(labels, null_prob, normalize=False)
    altDevianceLogit = -2.*alt_log_likelihood
    nullDevianceLogit = -2.*null_log_likelihood

    p_value_logit = scipy.stats.chi2.sf(nullDevianceLogit - altDevianceLogit, 1)

    AUC_ROC = roc_auc_score(~labels, data)

    if printResults:
        print('Logistic regression based method:', '\n\t',
                #'  Logistic regression based deviances and their difference test (chi2 p-value):\n\t', 
                #'Null deviance:\t', np.round(nullDevianceLogit, 1), '\n\t',
                #'Alternative deviance:\t', np.round(altDevianceLogit, 1), '\n\t',
                'p-value:\t', '%.1e' % p_value_logit, '\n\t',
                'AUC_ROC:', np.round(AUC_ROC, 2))

    return p_value_logit, AUC_ROC

def getHypergeometricTestPValue(unionSize, firstSubsetSize, secondSubsetSize, overlapSize):

    return scipy.stats.hypergeom(unionSize, firstSubsetSize, secondSubsetSize).sf(overlapSize-1)
