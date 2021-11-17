import numpy as np
import matplotlib.pyplot as plt
from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker
from matplotlib import rc
from scipy import optimize
from data_get import get_testdata
import scipy.stats as stats

def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
    # calculate mean for each row and phenotype mean
    matrix_means = np.nanmean(meth_matrix, axis=1)
    #print('meth_matrix',meth_matrix)
    #print('matrix means',matrix_means)
    phenotype_mean = np.mean(phenotype)
    #print('phenotype mean',phenotype_mean)
    
    # subtract means from observed values
    # check: https://stackoverflow.com/questions/18691084/what-does-1-mean-in-numpy-reshape
    transformed_matrix = meth_matrix - matrix_means.reshape([-1,1])
    #print('transformed matrix',transformed_matrix)
    transformed_phenotype = phenotype - phenotype_mean
    #print('trans phe',transformed_phenotype)
    
    # calculate covariance
    covariance = np.nansum(transformed_matrix * transformed_phenotype, axis=1)
    #print('cov',covariance)
    variance_meth = np.sqrt(np.nansum(transformed_matrix ** 2, axis=1))
    #print('var meth',variance_meth)
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))
    #print('var pheno',variance_phenotype)
    return covariance / (variance_meth * variance_phenotype)

def r2(x,y):
    # return r squared
    return stats.pearsonr(x,y)[0] **2

# Plot Prediction: known_ages
def plot_prediction(sites, known_ages, predicted_ages, label = None):
    # define optimization function
    def func(x, a, b):
        return a * np.asarray(x)**0.5 + b
    def func_lin(x, a, b):
        return a * np.asarray(x) + b
    def func_log(x, a, b):
        return a * np.log(x) + b
    
    # fit trend line
    popt, pcov = optimize.curve_fit(func, [1 + x for x in known_ages], predicted_ages)
    #print('popt',popt,pcov)
    popt_lin, pcov_lin = optimize.curve_fit(func_lin, [1 + x for x in known_ages], predicted_ages)
    popt_log, pcov_log = optimize.curve_fit(func_log, [1 + x for x in known_ages], predicted_ages)
    
    # get r squared
    rsquared = r2(predicted_ages, func([1 + x for x in known_ages], *popt))
    #print('rsquared',rsquared)
    rsquared_lin = r2(predicted_ages, func_lin([1 + x for x in known_ages], *popt_lin))
    rsquared_log = r2(predicted_ages, func_log([1 + x for x in known_ages], *popt_log))
    
    # initialize plt plot
    plt.subplots(figsize=(12,12))
    plt.plot(sorted(known_ages), func(sorted([1 + x for x in known_ages]), *popt), 'r--')
    plt.plot(sorted(known_ages), func_lin(sorted([1 + x for x in known_ages]), *popt_lin), 'b--')
    plt.plot(sorted(known_ages), func_log(sorted([1 + x for x in known_ages]), *popt_log), 'g--')
    
    # scatter plot
    plt.scatter(known_ages, predicted_ages, marker='o', alpha=0.8, color='g')
    
    # display formula
    eqn = 'f(x) = '+str(round(popt[0],2)) + r'$\sqrt{x}$' ' + ' + str(round(popt[1],2))+r'$, R^{2}=$'+str(round(rsquared,2))
    plt.text(70,-40, eqn,fontsize=20,color='r')
    eqn_lin = 'f(x) = '+str(round(popt_lin[0],2))+'x' ' + '+str(round(popt_lin[1],2))+r'$, R^{2}=$'+str(round(rsquared_lin,2))
    plt.text(70,-55, eqn_lin,fontsize=20,color='b')
    eqn_log = 'f(x) = '+str(round(popt_log[0],2))+r'$\log{x}$' ' + '+str(round(popt_log[1],2))+r'$, R^{2}=$'+str(round(rsquared_log,2))
    plt.text(70,-70, eqn_log,fontsize=20,color='g')
    
    # plt.scatter(predicted_ages, sites, marker='o', alpha=0.8, color='r')
    plt.title(label, fontsize=18)
    plt.xlabel('Chronological Age', fontsize=16)
    plt.ylabel('Epigenetic State', fontsize=16)
    plt.savefig('myplot.png')

def plot_testdata(file,pcc):
    # retrieve the training and testing data
    test_data = get_testdata(file)
    #print(test_data)
    # unpack the training and testing data
    test_samples, test_cpg_sites, test_ages, test_methylation_values = test_data
    #train_samples, train_cpg_sites, train_ages, train_methylation_values = train_data
    #print(test_methylation_values)
    #print('test ages',test_ages)
    abs_pcc_coefficients = abs(pearson_correlation(test_methylation_values, test_ages)) 
    
    # site selection based on user defined pcc values
    testing_sites = np.where(abs_pcc_coefficients > pcc)[0]
    
    # initialize the EPM model 
    epm = EpigeneticPacemaker(iter_limit=100, error_tolerance=0.00001)
    
    # fit the model using the data
    epm.fit(test_methylation_values[testing_sites,:], test_ages)
    
    # generate predicted ages using the test data
    test_sites = test_methylation_values[testing_sites,:]
    #print('test sites',test_sites)
    test_predict = epm.predict(test_sites)
    #print('predict',test_predict)
    #print('sites',len(test_sites))
    #print('predict',len(test_predict))
    #print('age',len(test_ages))
    #test_ages = test_ages[testing_sites]
    # print('test ages',test_ages)
    
    # plot the model results 
    plot_prediction(test_sites, test_ages, test_predict, "EpigeneticPacemaker")
    
    num_sites = len(test_sites)
    num_predicts = len(test_predict)
    return num_sites,num_predicts,test_ages,test_predict
