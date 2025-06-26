# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 13:08:09 2017

@author: harryholt

Load.py

Purpose:
    - Load the data
    - Clean (remove Nan values from profiles and depths)
    - Centre and Standardise
    - Print data to store the results

"""
# Importing modules
import h5py
import numpy as np
import pickle
from sklearn import preprocessing
import time
import os
from Load_ECCO import decode_llc, decode_llc_grid, transp_tiles

import Print

start_time = time.perf_counter()

def main(address, filename_raw_data, runIndex, subsample_uniform, subsample_random,\
         subsample_inTime, grid, conc, fraction_train, inTime_start, inTime_finish,\
         fraction_nan_samples, fraction_nan_depths, cov_type, run_bic=False, use_ecco=False):
    print("Starting Load.main")
    """ Main function for module"""
    if use_ecco:
        lon, lat, dynHeight, Tint, Sint, varTime = load_ecco(filename_raw_data)
    else:
        lon, lat, dynHeight, Tint, Sint, varTime = load(filename_raw_data)
    print("Removing depths with high NaN counts")
    Tint, Sint, depth = \
        removeDepthFractionNan(Tint, Sint, fraction_nan_depths)
    print("Removing profiles with high NaN counts")
    lon, lat, dynHeight, Tint, Sint, varTime = \
        removeSampleFractionNan(lon, lat, dynHeight, Tint, Sint, \
        varTime, fraction_nan_samples)
    print("Dealing with remaining NaN values")
    Tint, Sint = dealwithNan(Tint, Sint)
    
    # at this point the data has been successfully cleaned.
    """ now we need to subselect the training data """
    print("Selecting subset of data that will be used as training data")
    Tint_train, Sint_train, varTime_train = None, None, None
    if subsample_uniform: # currently working 
        lon_train, lat_train, dynHeight_train, Tint_train, Sint_train, \
            varTime_train = uniformTrain(lon, lat, dynHeight, Tint, Sint, \
            varTime, depth, grid, conc)
    if subsample_random: # also working
        lon_train, lat_train, dynHeight_train, Tint_train, Sint_train, \
            varTime_train = randomTrain(lon, lat, dynHeight, Tint, Sint, \
            varTime, depth, fraction_train)
    if subsample_inTime:
        lon_train, lat_train, dynHeight_train, Tint_train, Sint_train, \
            varTime_train = inTimeTrain(lon, lat, dynHeight, Tint, Sint, \
            varTime, depth, inTime_start, inTime_finish)
    
    ## At this point we should have a training data set to go with the full data set
    """ Now we can centre and standardise the training data, and the whole data will follow """
    # NOTE: currently unsure how to include Sint in the centring
    # I also change the nomenaculature at this point in the code
    """ It is important that the training data set initialises the standardised object !! """
    print("Centre and standardise the training dataset")
    stand, stand_store, varTrain_centre = centreAndStandardise(address, runIndex, Tint_train)
    var_centre = stand.transform(Tint)  # Centre the full dataset based on the training data set
    
    """ Create the test dataset, by subtracting the set(var_train) from set(var)"""
    """
    lon_test = [x for x in lon if x not in set(lon_train)]
    lat_test = [x for x in lat if x not in set(lat_train)]
    dynHeight_test = [x for x in dynHeight if x not in set(dynHeight_train)]
    Tint_test = [x for x in Tint if x not in set(Tint_train)]
    Sint_test = [x for x in Sint if x not in set(Sint_train)]
    varTime_test = [x for x in varTime if x not in set(varTime_train)]
    
    varTest_centre = stand.transform(Tint_test)
    """
    # INFORMATION
    # varTrain_centre stores the training, standardised
    # var_centre stores the full standardised data set
    # varTest_centre stores the test dataset
    
    """ Now we can print the results of this process to a file for later use """
#    print("Starting Print")
    print("varTrain_centre.shape = ", varTrain_centre.shape)
    if not run_bic:
        Print.printLoadToFile(address, runIndex, lon, lat, dynHeight,\
                              Tint, var_centre, Sint, varTime, depth)
        Print.printLoadToFile_Train(address, runIndex, lon_train, \
                          lat_train, dynHeight_train, Tint_train, \
                          varTrain_centre, Sint_train, varTime_train, depth)
        Print.printDepth(address, runIndex, depth)
    
    if run_bic:
        return lon_train, lat_train, dynHeight_train, Tint_train, \
               varTrain_centre, Sint_train, varTime_train
    
###############################################################################
# Functions which Main uses
def load(filename_raw_data):
    print("Load.load")
    """ This function loads the raw data from a .mat file """

    # declare variables as empty arrays
    lon, lat, dynHeight, Tint, Sint, varTime = [], [], [], [], [], []

    # import lon, lat, Temperature, Salinity, times, and dynamic height explicitly
    variables = ['lon', 'lat', 'dynht300_1500', 'Tint', 'Sint', 'dectime']
    mat = h5py.File(filename_raw_data, variable_names = variables)
    lon = mat["lon"]
    lon = np.array(lon)[:,0]
    lat = mat["lat"]
    lat = np.array(lat)[:,0]
    dynHeight = mat["dynht300_1500"]
    dynHeight = np.array(dynHeight)[:,0]
    Tint = mat["Tint"]
    Tint = np.array(Tint)
    Sint = mat["Sint"]
    Sint = np.array(Sint)
    varTime = mat["dectime"]
    varTime = np.array(varTime)
    
#   print("Shape of variable = ", Tint.shape)
#    print("axis 0 = ", np.ma.size(Tint, axis=0)) # axis 0 should be ~290520
#    print("axis 1 = ", np.ma.size(Tint, axis=1)) # axis 1 should be ~400
    
    return lon, lat, dynHeight, Tint, Sint, varTime

def load_ecco(base_path, dic_file="DIC.0000002232.data", num_variables=15):
    """Load ECCO model output in llc270 tile format."""
    filename = os.path.join(base_path, dic_file)
    data_raw = decode_llc(filename, num_variables)
    data_dic = data_raw[0]
    for ind in range(data_dic.shape[0]):
        data_dic[ind] = transp_tiles(data_dic[ind])
    data_dic[data_dic == 0] = np.nan

    xc = decode_llc_grid(os.path.join(base_path, 'XC.data'), 2)
    xc = transp_tiles(xc)
    yc = decode_llc_grid(os.path.join(base_path, 'YC.data'), 2)
    yc = transp_tiles(yc)

    lon = xc.flatten()
    lat = yc.flatten()
    Tint = data_dic.reshape(data_dic.shape[0], -1).T
    Sint = np.full_like(Tint, np.nan)
    dynHeight = np.full(lon.shape, np.nan)
    varTime = np.arange(lon.size, dtype=float)
    return lon, lat, dynHeight, Tint, Sint, varTime

def removeDepthFractionNan(VAR, VAR2, fraction_of):
    print("Load.removeDepthFractionNan")
    """ This function removes all depths will a given number of Nan values """

    delete_depth = []
    
    # create a pressure array to record the pressure levels dept
    depth_remain = 5 * np.arange(np.ma.size(VAR, axis=1)) 
    
    for i in range(np.ma.size(VAR, axis=1)):
        var_mean_nan = []
        var_nan_sum = 0
        var_mean_nan = np.isnan(VAR[:,i])
        var_nan_sum = var_mean_nan.sum()
        if var_nan_sum >= np.ma.size(VAR, axis=0)/fraction_of:
            delete_depth.append(i)
            #print(i," DELETED")
            
    delete_depth = np.squeeze(np.asarray(delete_depth))
    #print(delete_depth.shape)
    VAR = np.delete(VAR, (delete_depth), axis=1 )
    VAR2 = np.delete(VAR2, (delete_depth), axis=1 )
    depth_remain = np.delete(depth_remain, (delete_depth))
    print("Number of depths deleted above the 1/"+\
        str(fraction_of)+" criterion = ", delete_depth.size)
    #print("Depth pressure remain = ", depth_remain.shape)
    #print("VAR shape after half Sample removed = ", VAR.shape)
    return VAR, VAR2, depth_remain  

def removeSampleFractionNan(LON, LAT, dynHeight, VAR, VAR2, varTime, fraction_of):
    print("Load.removeSampleFractionNan")
    """ This function removes all profiles will a given number of Nan values """
    delete_sample = []
    sample_count = 0
    for i in range(np.ma.size(VAR, axis=0)):
        var_mean_nan = []
        var_mean_nan = np.isnan(VAR[i,:])
        var_nan_sum = var_mean_nan.sum()        
        if var_nan_sum >= np.ma.size(VAR, axis=1)/fraction_of:
            sample_count = sample_count + 1
            delete_sample.append(i)
            
    delete_sample = np.squeeze(np.asarray(delete_sample))
    #print(delete_sample.shape)
    VAR = np.delete(VAR, (delete_sample), axis=0 )
    VAR2 = np.delete(VAR2, (delete_sample), axis=0 )
    varTime = np.delete(varTime, (delete_sample))
    LON = np.delete(LON, (delete_sample))
    LAT = np.delete(LAT, (delete_sample)) 
    dynHeight = np.delete(dynHeight, (delete_sample))
    print("Number of samples deleted above the 1/"+\
        str(fraction_of)+" criterion = ", sample_count)
    #print("VAR shape after half Col removed = ", VAR.shape)
    #print("LON shape = ", LON.shape)
    #print("LAT shape = ", LAT.shape)
    return LON, LAT, dynHeight, VAR, VAR2, varTime

def dealwithNan(VAR, VAR2):
    print("Load.dealwithNan")
    """ Simple Linear interpolator to deal with th eremaining Nan valus """
    print("total number of values before nan interpolation = ", np.size(VAR))
    print("number of nans before nan interpolation = ",np.isnan(VAR).sum())

    # loop over array
    for i in range(np.ma.size(VAR, axis=0)):
        mask = []
        VAR_temp = []
        VAR_temp = VAR[i,:]
        mask = np.isnan(VAR_temp)
        VAR_temp[mask] = np.interp(np.flatnonzero(mask), \
                         np.flatnonzero(~mask), VAR_temp[~mask])
        VAR[i,:] = VAR_temp

    # loop over array    
    for i in range(np.ma.size(VAR2, axis=0)):
        mask2 = []
        VAR2_temp = []
        VAR2_temp = VAR2[i,:]
        mask2 = np.isnan(VAR2_temp)
        VAR2_temp[mask2] = np.interp(np.flatnonzero(mask2), \
                           np.flatnonzero(~mask2), VAR2_temp[~mask2])
        VAR2[i,:] = VAR2_temp
    print("number of nans after = ",np.isnan(VAR).sum())
    return VAR, VAR2

###############################################################################

def uniformTrain(lon, lat, dynHeight, VAR, VAR2, varTime, depth, grid, concentration):
    print("Load.uniformTrain")
    # array of integegers from -180 to +180
    i_array = np.arange(-180, 180, grid ,dtype=np.int32) 
    array_lon, array_lat, array_dynHeight, array_time = [], [], [], []
    count = 0
    # Loop over the Longitudinal values
    for i in i_array:
        indices_lon = []
        # indices of lon that are in the interval [i,i+grid]
        indices_lon = np.squeeze( np.nonzero((lon>=i)&(lon<i+grid) ) )  
        lon_temp_i = lon[indices_lon] # lon values that match the criteria
        lat_temp_i = lat[indices_lon] # lat values that match the criteria
        dynHeight_temp_i = dynHeight[indices_lon] # dynamic heights
        time_temp_i = varTime[indices_lon] # times that match
        
        if indices_lon.size == 0:
            #print("Indices_lon is null for range ", i, " + ", grid)
            continue
        j_array = []
        j_array = np.arange(min(np.floor(lat[indices_lon])), \
                  max(np.ceil(lat[indices_lon])), grid ,dtype=np.int32) 
        
        # Loop over the Latitudes for a given longitude
        for j in j_array:
            indices_lat = []
            indices_lat = np.squeeze( np.nonzero((lat_temp_i >= j) & \
                              (lat_temp_i < j + grid)) ) 
            
            lon_temp_j = lon_temp_i[indices_lat] # lon values
            lat_temp_j = lat_temp_i[indices_lat] # lat values 
            dynHeight_temp_j = dynHeight_temp_i[indices_lat] # dynamic heights
            time_temp_j = time_temp_i[indices_lat] # times
            #print("Indices: ", indices_lat.size , indices_lat)
            
            if indices_lat.size == 0:
                #print("Indices_lat is null for range ",i , j, " + ",grid)
                continue
            if indices_lat.size == 1:
                #print("Converting # into [#] for indices_lat")
                indices_lat = np.asarray([indices_lat])
                lon_temp_j = np.asarray([lon_temp_j])
                lat_temp_j = np.asarray([lat_temp_j])
                time_temp_j = np.asarray([time_temp_j])
                dynHeight_temp_j = np.asarray([dynHeight_temp_j])
                
            select = None   # Reset the select variable for every iteration
            select = np.random.randint(indices_lat.size, size = concentration)
            
            # Select the training sample
            lon_select = lon_temp_j[select] # single lon value for selected sample
            lat_select = lat_temp_j[select] # single lat value for selected sample
            dynHeight_select = dynHeight_temp_j[select] # single dynamic height value
            time_select = time_temp_j[select] # single time value for selected sample
            
            for col_depth in range(0,len(depth)):
                variable_col, variable2_col = [], []
                variable_col = VAR[:,col_depth]
                variable2_col = VAR2[:,col_depth]
                variable_train_lon = variable_col[indices_lon]
                variable2_train_lon = variable2_col[indices_lon]
                variable_train_lat = variable_train_lon[indices_lat]
                variable2_train_lat = variable2_train_lon[indices_lat]
                variable_train = variable_train_lat[select] # single depth sample
                variable2_train = variable2_train_lat[select]
                
                if col_depth == 0:
                    var_train_sample = variable_train.reshape(variable_train.size,1)
                    var2_train_sample = variable2_train.reshape(variable2_train.size,1)
                else:
                    var_train = variable_train.reshape(variable_train.size,1)
                    var_train_sample = np.hstack([var_train_sample, var_train])
                    
                    var2_train = variable2_train.reshape(variable2_train.size,1)
                    var2_train_sample = np.hstack([var2_train_sample, var2_train])
                    
                # var_train_sample has whole sample selected from criteria

            array_lon = np.append(array_lon,lon_select)
            array_lat = np.append(array_lat,lat_select)
            array_dynHeight = np.append(array_dynHeight, dynHeight_select)
            array_time = np.append(array_time, time_select)
            
            if i == i_array[0] and j == j_array[0]:
                var_train_array = var_train_sample
                var2_train_array = var2_train_sample
            else:
                var_train_array = np.vstack([var_train_array,var_train_sample])
                var2_train_array = np.vstack([var2_train_array,var2_train_sample])
            
            count = count + concentration
    
    array_lon = np.asarray(array_lon)
    array_lon = np.squeeze(array_lon.reshape(1,array_lon.size))
    array_lat = np.asarray(array_lat)
    array_lat = np.squeeze(array_lat.reshape(1,array_lat.size))
    array_dynHeight = np.asarray(array_dynHeight)
    array_dynHeight = np.squeeze(array_dynHeight.reshape(1,array_dynHeight.size))
    array_time = np.asarray(array_time)
    array_time = np.squeeze(array_time.reshape(1,array_time.size))    

    print("var_train_array.shape = ", var_train_array.shape)
    
    return array_lon, array_lat, array_dynHeight, var_train_array, \
           var2_train_array, array_time

###############################################################################
def randomTrain(lon, lat, dynHeight, Tint, Sint, varTime, depth, fraction_train):
    lon_rand, lat_rand, dynHeight_rand, Tint_rand, Sint_rand, varTime_rand = \
        None, None, None, None, None, None, None

    array_size = np.ma.size(Tint, axis = 0)   # expecting around 280,000
    number_rand = fraction_train * array_size
    
    indices_rand = None
    indices_rand = np.random.randint(0, high=array_size, size = number_rand)
    
    lon_rand = lon[indices_rand]
    lat_rand = lat[indices_rand]
    dynHeight_rand = dynHeight[indices_rand]
    Tint_rand = Tint[indices_rand, :]
    Sint_rand = Sint[indices_rand, :]
    varTime_rand = varTime[indices_rand]
    
    print("Tint_rand.shape = ", Tint_rand.shape)
    
    return lon_rand, lat_rand, dynHeight_rand, Tint_rand, Sint_rand, varTime_rand
###############################################################################
def inTimeTrain(lon, lat, dynHeight, Tint, Sint, varTime, depth, \
                inTime_start, inTime_finish):
    lon_train, lat_train, dynHeight_train, Tint_train, Sint_train, \
    varTime_train = None, None, None, None, None, None
    
    indices_time = None
    indices_time = np.logical_and(varTime > inTime_start, \
                                  varTime < inTime_finish)
    
    lon_train = lon[indices_time]
    lat_train = lat[indices_time]
    dynHeight_train = dynHeight[indices_time]
    Tint_train = Tint[indices_time, :]
    Sint_train = Sint[indices_time, :]
    varTime_train = varTime[indices_time]
    
    print("Tint_train.shape = ", Tint_train.shape)
    
    return lon_train, lat_train, dynHeight_train, Tint_train, \
           Sint_train, varTime_train

###############################################################################

def centreAndStandardise(address, runIndex, VAR):
    print("Load.centreAndStandardise")
    """ Function to creat a standardised object using the training data set """
    stand_store = address+"Objects/Scale_object.pkl"
    stand = preprocessing.StandardScaler()
    # VAR = Tint_train so that stand is initiaised with the traing dataset
    stand.fit(VAR)  
    var_stand = stand.transform(VAR)
    
    with open(stand_store, 'wb') as output:
        ScaleObject = stand
        pickle.dump(ScaleObject, output, pickle.HIGHEST_PROTOCOL)
    del ScaleObject
    
    return stand, stand_store, var_stand
    
print('Load runtime = ', time.perf_counter() - start_time,' s')
