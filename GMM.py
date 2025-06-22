# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 10:59:59 2017

@author: harryholt

GMM.py

Purpose:
    - Generate the GMM object in reduced space using X_train
    - Retaing information on the classes in reduced and reconstructed space
    - Assigned labels and posterior probabilities to the full training dataset X
    
"""
import pickle
from sklearn import mixture
import numpy as np
import time
import ClassProperties
import Print

start_time = time.perf_counter()

def create(address, runIndex, n_comp, cov_type):
    print("GMM.create")
    """ Takes the training dataset and creates the GMM object """
    # load col_reduced
    col_reduced = None
    col_reduced = Print.readColreduced(address, runIndex)
    col_reduced_array = np.arange(col_reduced)
    
    # load training data in reduced pca space
    lon_train, lat_train, dynHeight_train, X_train_array, varTime_train = \
      None, None, None, None, None
    lon_train, lat_train, dynHeight_train, X_train_array, varTime_train = \
      Print.readPCAFromFile_Train(address, runIndex, col_reduced)
    
    # calculate GMM Object
    gmm, gmm_weights, gmm_means, gmm_covariances = \
      None, None, None, None
    gmm, gmm_weights, gmm_means, gmm_covariances = \
      GaussianMixtureModel(address, runIndex, n_comp, X_train_array, cov_type)
    
    """ Print the information on the classes to a file """
    class_number_array = np.arange(0,n_comp).reshape(-1,1)
    Print.printGMMclasses(address, runIndex, class_number_array, gmm_weights, \
                          gmm_means, gmm_covariances, col_reduced_array, 'reduced')
    
###############################################################################
def apply(address, runIndex, n_comp):
    print("GMM.apply")
    # load col_reduced value
    col_reduced = None
    col_reduced = Print.readColreduced(address, runIndex)
    
    # load full data array - X
    lon, lat, dynHeight, X_array, varTime = None, None, None, None, None
    lon, lat, dynHeight, X_array, varTime = Print.readPCAFromFile(address, \
      runIndex, col_reduced)
    
    # load GMM object
    gmm = None
    with open(address+'Objects/GMM_Object.pkl', 'rb') as input:
        gmm = pickle.load(input)
    
    # calculate the labels and probabilities of the profiles
    sortedLabels, labels, post_prob = None, None, None

    # predict classes based on training dataset, output labels
    labels = gmm.predict(X_array)  # expected shape (n_profiles, pca index)
    Print.printLabelsUnsorted(address, runIndex, lon, lat, dynHeight, varTime, labels)

    # sort labels by mean SST of each class
    # (do this during the "plot" stage now
#   allDF, sortedLabels, old2new = ClassProperties.main(address,runIndex,n_comp) 
#   Print.printLabels(address, runIndex, lon, lat, dynHeight, varTime, sortedLabels)
   
    # calculate posterior probabilities
    post_prob = gmm.predict_proba(X_array) # expected shape (n_profiles, classes)

    # needed for input of printPosteriorProb
    class_number_array = np.arange(0,n_comp).reshape(-1,1)
    
    # Print Labels and probabilities to file
    Print.printPosteriorProb(address, runIndex, lon, lat, dynHeight, \
                             varTime, post_prob, class_number_array)
    
def GaussianMixtureModel(address, runIndex, n_comp, X_train, cov_type):
    print("GMM.GaussianMixtureModel")
    gmm = None
    gmm = mixture.GaussianMixture(n_components = n_comp, \
                                  covariance_type = cov_type)
#    gmm = mixture.BayesianGaussianMixture(n_components = n_comp, \
#                                          covariance_type = cov_type)

    # use training dataset to "fit" Gaussian mixture model
    gmm.fit(X_train)
    
    # store the GMM object
    gmm_store = address+"Objects/GMM_Object.pkl"
    with open(gmm_store, 'wb') as output:
        gmmObject = gmm
        pickle.dump(gmmObject, output, pickle.HIGHEST_PROTOCOL)
    del gmmObject
    
    # means and covariances
    weights, means, covariances = None, None, None
    weights = np.squeeze(gmm.weights_)  # shape (n_components)
    # note: "weights" is the same for each col_red
    means = np.squeeze(gmm.means_)  # shape (n_components, n_features)
    covariances = abs(np.squeeze(gmm.covariances_))  # shape (n_components, n_features) 
    
    return gmm, weights, means, covariances

print('GMM runtime = ', time.perf_counter() - start_time,' s')
