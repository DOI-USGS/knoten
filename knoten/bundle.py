import numpy as np
import pandas as pd
import pvl
import os
import csmapi
import itertools
from math import floor

from pysis import isis
from ale.drivers import loads
from collections import OrderedDict
from knoten.csm import create_csm

def generate_sensors(cubes, directory=None, clean=False):
    """
    Generate a set of USGSCSM sensor models from a list of ISIS cube files

    Parameters
    ----------
    cubes     : str
              Directory/filename of a file containing ISIS cube file paths
    directory : str
              Output directory to save resulting json files. Defaults to the
              same directory as cube list file
    clean     : flag
              Option to delete json file outputs

    Returns
    -------
    sensors   : dictionary
              Dictionary mapping ISIS serial numbers to USGSCSM sensor models
    """
    if directory is None:
        directory = os.path.dirname(cubes)

    isd_files = []
    sensors = {}
    for line in open(cubes):
        basename = os.path.splitext(os.path.basename(line.strip()))[0]
        isd = os.path.join(directory, basename+'.json')
        isd_files.append(isd)
        with open(isd, 'w+') as f:
            f.write(loads(line.strip(), formatter='usgscsm'))

        sn = isis.getsn(from_=line.strip()).strip().decode('utf-8')
        sensors[sn] = create_csm(isd)

    if clean:
        for isd in isd_files:
            os.remove(isd)

    return sensors

def closest_approach(points, direction):
    """
    Compute the point of closest approach between lines.

    Parameters
    ----------
    points : ndarray
             An n x 3 array of points on each line
    direction : ndarray
                An n x 3 array of vectors that defines the direction of each line

    Returns
    -------
     : array
       The (x, y, z) point that is closest to all of the lines
     : ndarray
       The (x, y, z) covariance matrix that describes the uncertaintly of the
       point
    """
    num_lines = points.shape[0]
    design_mat = np.zeros((num_lines * 3, 3))
    rhs = np.zeros(num_lines * 3)
    for i in range(num_lines):
        point = points[i]
        line = direction[i] / np.linalg.norm(direction[i])
        design_mat[3*i:3*i+3] = np.identity(3) - np.outer(line, line)
        rhs[3*i:3*i+3] = np.dot(point,line) * line + point
    N_inv = np.linalg.inv(design_mat.T.dot(design_mat))
    closest_point = N_inv.dot(design_mat.T).dot(rhs)

    return closest_point, N_inv

def compute_apriori_ground_points(network, sensors):
    """
    Compute a priori ground points for all of the free points in a control network.

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio
    sensors : dict
              A dictionary that maps ISIS serial numbers to CSM sensors

    Returns
    -------
     : DataFrame
       The control network dataframe with updated ground points
    """
    for point_id, group in network.groupby('id'):
        # Free points are type 2 for V2 and V5 control networks
        if group.iloc[0]["pointType"] != 2:
            continue
        positions = []
        look_vecs = []
        for measure_id, row in group.iterrows():
            measure = csmapi.ImageCoord(row["line"], row["sample"])
            locus = sensors[row["serialnumber"]].imageToRemoteImagingLocus(measure)
            positions.append([locus.point.x, locus.point.y, locus.point.z])
            look_vecs.append([locus.direction.x, locus.direction.y, locus.direction.z])
        ground_pt, covar_mat = closest_approach(np.array(positions), np.array(look_vecs))
        covar_vec = [covar_mat[0,0], covar_mat[0,1], covar_mat[0,2],
                     covar_mat[1,1], covar_mat[1,2], covar_mat[2,2]]
        network.loc[network.id == point_id, ["aprioriX", "aprioriY", "aprioriZ"]] = ground_pt
        network.loc[network.id == point_id, ["adjustedX", "adjustedY", "adjustedZ"]] = ground_pt
        network.loc[network.id == point_id, ["adjustedX", "adjustedY", "adjustedZ"]] = ground_pt
        # We have to do a separate loop to assign a list to a single cell
        for measure_id, row in group.iterrows():
            network.at[measure_id, 'aprioriCovar'] = covar_vec
    return network

class CsmParameter:
    """
    Container class that describes a parameter for a CSM sensor model
    """

    def __init__(self, sensor, index):
        self.index = index
        self.name = sensor.getParameterName(index)
        self.type = sensor.getParameterType(index)
        self.units = sensor.getParameterUnits(index)
        self.value = sensor.getParameterValue(index)

    def __repr__(self):
        return f'{self.index} {self.name.strip()} ({self.type}): {self.value} {self.units}'

def get_sensor_parameters(sensor, set="adjustable"):
    """
    Get a set of the CSM parameters for a CSM sensor

    Parameters
    ----------
    sensor : CSM sensor
             The CSM sensor model
    set : str
          The CSM parameter set to get. Either valid, adjustable, or non_adjustable

    Returns
    -------
     : List
       A list of CsmParameters
    """
    if (set.upper() == "VALID"):
        param_set = csmapi.VALID
    elif (set.upper() == "ADJUSTABLE"):
        param_set = csmapi.ADJUSTABLE
    elif (set.upper() == "NON_ADJUSTABLE"):
        param_set = csmapi.NON_ADJUSTABLE
    else:
        raise ValueError(f"Invalid parameter set \"{set}\".")
    return [CsmParameter(sensor, index) for index in sensor.getParameterSetIndices(param_set)]

def compute_sensor_partials(sensor, parameters, ground_pt):
    """
    Compute the partial derivatives of the line and sample with respect to a set
    of parameters.

    Parameters
    ----------
    sensor : CSM sensor
             The CSM sensor model
    parameters : list
                 The list of  CsmParameter to compute the partials W.R.T.
    ground_pt : array
                The (x, y, z) ground point to compute the partial derivatives at

    Returns
    -------
     : ndarray
       The 2 x n array of partial derivatives. The first array is the line
       partials and the second array is the sample partials. The partial
       derivatives are in the same order as the parameter list passed in.
    """
    partials = np.zeros((2, len(parameters)))
    csm_ground = csmapi.EcefCoord(ground_pt[0], ground_pt[1], ground_pt[2])
    for i in range(len(parameters)):
        partials[:, i] = sensor.computeSensorPartials(parameters[i].index, csm_ground)
    return partials

def compute_ground_partials(sensor, ground_pt):
    """
    Compute the partial derivatives of the line and sample with respect to a
    ground point.

    Parameters
    ----------
    sensor : CSM sensor
             The CSM sensor model
    ground_pt : array
                The (x, y, z) ground point to compute the partial derivatives W.R.T.

    Returns
    -------
     : ndarray
       The 2 x 3 array of partial derivatives. The first array is the line
       partials and the second array is the sample partials. The partial
       derivatives are in (x, y, z) order.
    """
    csm_ground = csmapi.EcefCoord(ground_pt[0], ground_pt[1], ground_pt[2])
    partials = np.array(sensor.computeGroundPartials(csm_ground))
    return np.reshape(partials, (2, 3))

def compute_coefficient_columns(network, sensors, parameters):
    """
    Compute the columns for different coefficients

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio.
    sensors : dict
              Dictionary that maps ISIS serial numbers to CSM sensor models
    parameters : dict
                 Dictionary that maps serial numbers to lists of parameters to
                 solve for.

    Returns
    -------
     : OrderedDict
       Dictionary that maps serial numbers and point IDs to the column range
       their parameters are in the Jacobian matrix.
    """
    num_columns = 0
    coefficient_columns = OrderedDict()
    for serial in network["serialnumber"].unique():
        coefficient_columns[serial] = num_columns
        num_columns += len(parameters[serial])
        coefficient_columns[serial] = (coefficient_columns[serial], num_columns)
    for point_id in network["id"].unique():
        # Skip fixed points
        if network.loc[network.id == point_id].iloc[0]["pointType"] == 4:
            continue
        coefficient_columns[point_id] = num_columns
        num_columns += 3
        coefficient_columns[point_id] = (coefficient_columns[point_id], num_columns)
    return coefficient_columns

def compute_jacobian(network, sensors, parameters, coefficient_columns):
    """
    Compute the Jacobian matrix.

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio.
    sensors : dict
              Dictionary that maps ISIS serial numbers to CSM sensor models
    parameters : dict
                 Dictionary that maps serial numbers to lists of parameters to
                 solve for.
    coefficient_columns : OrderedDict
                          Dictionary that maps serial numbers and point IDs to
                          the column range their parameters are in the Jacobian
                          matrix.

    Returns
    -------
     : ndarray
       The Jacobian matrix
    """
    num_columns = max([col_range[1] for col_range in coefficient_columns.values()])
    num_rows = len(network) * 2
    jacobian = np.zeros((num_rows, num_columns))

    for i in range(len(network)):
        row = network.iloc[i]
        serial = row["serialnumber"]
        ground_pt = row[["adjustedX", "adjustedY", "adjustedZ"]]
        sensor = sensors[serial]
        params = parameters[serial]
        image_range = coefficient_columns[serial]
        point_range = coefficient_columns[row["id"]]
        jacobian[2*i : 2*i+2, image_range[0] : image_range[1]] = compute_sensor_partials(sensor, params, ground_pt)
        jacobian[2*i : 2*i+2, point_range[0] : point_range[1]] = compute_ground_partials(sensor, ground_pt)

    return jacobian

def compute_parameter_weights(network, sensors, parameters, coefficient_columns):
    """
    Compute the parameter weight matrix

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio.
    sensors : dict
              Dictionary that maps ISIS serial numbers to CSM sensor models
    parameters : dict
                 Dictionary that maps serial numbers to lists of parameters to
                 solve for.
    coefficient_columns : OrderedDict
                          Dictionary that maps serial numbers and point IDs to
                          the column range their parameters are in the Jacobian
                          matrix. Their parameters weights will have the same
                          ordering in the weight matrix.

    Returns
    -------
     : ndarray
       The parameter weight matrix
    """
    num_params = max([col_range[1] for col_range in coefficient_columns.values()])
    weight_mat = np.zeros((num_params, num_params))

    # Image parameters
    for sn, params in parameters.items():
        param_count = len(params)
        covar_mat = np.zeros((param_count, param_count))
        for a, b in itertools.product(range(param_count), range(param_count)):
            covar_mat[a, b] = sensors[sn].getParameterCovariance(params[a].index, params[b].index)
        col_range = coefficient_columns[sn]
        weight_mat[col_range[0]:col_range[1], col_range[0]:col_range[1]] = np.linalg.inv(covar_mat)

    # Point parameters
    for point_id, group in network.groupby('id'):
        ## If there is no covariance matrix, then just continue on
        point_covar = list(group.iloc[0]["aprioriCovar"])
        if len(point_covar) != 6:
            continue
        # The covariance matrix is stored as just one triangle, so we have
        # to unpack it.
        if len(point_covar) == 6:
            covar_mat = np.array(
                [[point_covar[0], point_covar[1], point_covar[2]],
                 [point_covar[1], point_covar[3], point_covar[4]],
                 [point_covar[2], point_covar[4], point_covar[5]]]
            )
            col_range = coefficient_columns[point_id]
            weight_mat[col_range[0]:col_range[1], col_range[0]:col_range[1]] = np.linalg.inv(covar_mat)

    return weight_mat

def compute_residuals(network, sensors):
    """
    Compute the error in the observations by taking the difference between the
    ground points groundToImage projections and measure values.

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio
    sensors : dict
             A dictionary that maps ISIS serial numbers to CSM sensors

    Returns
    -------
    V : np.array
       The control network dataframe with updated ground points
    """
    num_meas = len(network)
    V = np.zeros((num_meas, 2))

    for i in range(num_meas):
        row = network.iloc[i]
        serial = row["serialnumber"]
        ground_pt = row[["adjustedX", "adjustedY", "adjustedZ"]].values
        ground_pt = csmapi.EcefCoord(ground_pt[0], ground_pt[1], ground_pt[2])
        sensor = sensors[serial]
        img_coord = sensor.groundToImage(ground_pt)
        V[i,:] = [row['line'] - img_coord.line, row['sample'] - img_coord.samp]

    V = V.reshape(num_meas*2)
    return V

def compute_sigma(V, W_parameters, W_observations):
    """
    Computes the resulting standard deviation of the residuals for the current state of the bundle network.

    Parameters
    ----------
    V  :  np.array
          The control network dataframe with updated ground points
    W_parameters  :  ndarray
                     The parameter weight matrix (i.e.: sensor parameters and point weights)
    W_observations  :  ndarray
                     The observation weight matrix (i.e.: point weights)

    Returns
    -------
       : float64
         Standard deviation of the residuals

    """
    num_parameters = W_parameters.shape[0]
    num_observations = W_observations.shape[0]
    dof = num_observations - num_parameters
    VTPV = (V.dot(W_observations).dot(V))
    sigma0 = np.sqrt(VTPV/dof)
    return sigma0

def bundle_iteration(J, V, W_parameters, W_observations):
    """
    Parameters
    ----------
    J  :  ndarray
          The control network as a dataframe generated by plio.
    V  :  np.array
          The control network dataframe with updated ground points
    W_parameters  :  ndarray
                     The parameter weight matrix (i.e.: sensor parameters and point weights)
    W_observations  :  ndarray
                     The observation weight matrix (i.e.: measure weights)

    Returns
    -------
    N  :
    """

    N = J.T.dot(W_observations).dot(J) + W_parameters
    C = J.T.dot(W_observations).dot(V)
    dX = np.linalg.inv(N).dot(C)
    return N, dX

# For data snooping we need to calculate updated residuals
def compute_normalized_residual(J, V, N, W_parameters, W_observations):
    """
    Computes the normalized residual statistic for the data snooping method. Method derived from
    Forstner 1985 "The Reliability of Block Triangulation"

    Parameters
    ----------
    V  :  np.array
          The control network dataframe with updated ground points
    N  :

    W_parameters  :  ndarray
                     The parameter weight matrix (i.e.: sensor parameters and point weights)
    W_observations  :  ndarray
                     The observation weight matrix (i.e.: point weights)

    Returns
    -------
       : np.array
         Normalized residual statistic for the data snooping

    """
    sigma0 = compute_sigma(V, W_parameters, W_observations)
    Qxx = np.linalg.inv(N)
    Qvv = np.linalg.inv(W_observations) - J.dot(Qxx).dot(J.T)
    qvv = np.diagonal(Qvv)
    sigma_vi = sigma0*np.sqrt(qvv)
    wi = -V/sigma_vi

    return wi

def check_network(network):
    """
    Check that all control points in a network have at least 2 remaining measures.

    Parameters
    ----------
    network : DataFrame
              The control network as a dataframe generated by plio

    Returns
    -------
     : list
       List of measure indices that were masked out for being the only measure on a point.
    """
    bad_measures = []
    for point_id, group in network.groupby('id'):
        if len(group) < 2:
            for measure_index, _ in group.iterrows():
                bad_measures.append(measure_index)
    return bad_measures

def data_snooping(network, sensors, parameters, k=3.29, verbose=True):
    """
    Parameters
    ----------
    network  :  DataFrame
                The control network as a dataframe generated by plio
    sensors  :  dict
                A dictionary that maps ISIS serial numbers to CSM sensors
    parameters  : list
                 The list of  CsmParameter to compute the partials W.R.T.
    k  :  float64
          Critical value used for rejection criteria; defaults to Forstner's 3.29
          (or Baarda's 4.1??)
    verbose : bool
              If status prints should happen

    Returns
    -------
      :  list
      Indices of the network DataFrame that were rejected during data snooping
    """
    net = network
    net['mask'] = False

    rejected_indices = []
    awi = np.array([5, 5, 5, 5]) #initialize larger than k so you get into first iteration
    while (awi > k).any():

        # weight matrices
        coefficient_columns = compute_coefficient_columns(net[~net['mask']], sensors, parameters)
        num_parameters = max(col_range[1] for col_range in coefficient_columns.values())
        W_parameters = compute_parameter_weights(net[~net['mask']], sensors, parameters, coefficient_columns)
        num_observations = 2 * len(net[~net['mask']])
        W_observations = np.eye(num_observations)

        # bundle iteration (and set up)
        V = compute_residuals(net[~net['mask']], sensors)
        J = compute_jacobian(net[~net['mask']], sensors, parameters, coefficient_columns)
        sigma0 = compute_sigma(V, W_parameters, W_observations)
        N, dX = bundle_iteration(J, V, W_parameters, W_observations)

        # calculate test statistic
        wi = compute_normalized_residual(J, V, N, W_parameters, W_observations)
        awi = abs(wi)

        #find maximum
        imax = np.argmax(awi)
        if verbose:
            print(f'max wi = {awi[imax]}') # display

        if awi[imax] <= k:
            if verbose:
                print('Data Snooping Outlier Rejection Complete')
            break

        reject_index = floor(imax/2)
        reject = net.index[~net['mask']][reject_index]
        net.loc[reject, ['mask']] = True
        rejected_indices.append(reject)
        if verbose:
            print(f'max wi index = {imax}')
            print(f'max wi measure index = {reject_index}')
            print(f'rejecting measure {net.loc[reject, ["id", "serialnumber"]].values}')

        not_enough_measures = check_network(net[~net['mask']])
        if (not_enough_measures):
            for measure_index in not_enough_measures:
                if verbose:
                    print(f'single measure point {net.loc[measure_index, "id"]}')
                    print(f'rejecting measure {net.loc[measure_index, ["id", "serialnumber"]].values}')
                net.loc[measure_index, ['mask']] = True

        if verbose:
            print('')

    return rejected_indices
