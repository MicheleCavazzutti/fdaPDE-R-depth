## This file is part of fdaPDE, a C++ library for physics-informed
## spatial and functional data analysis.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

.DepthModelCtr <- setRefClass(
  Class = "DepthModel",
  fields = c(
    cpp_model  = "ANY", ## C++ model backend
    phi_function = "ANY", ## This phi function will be used to compute the weight function, when needed. If left to NULL,the identity function will be employed
    # Note that here I will eventually store objects used in future to save the C++ computations, even if it is not really required with this structure
    fun_representation_storage = "ANY", ## This matrix contains the Voronoi representation of the functions already computed. Will be used for the predict.
    # Note: probably the storage is not needed anymore, since we have already available the C++ matrix. In the future we will evaluate its elimination
  ),
  methods = c(
    fit = function() { # function(...) {
      # For the moment I assume that cpp_model has already been filled with all the data needed
      cpp_model$init() # This activates the model, see R_cpp to understand how other people behave with this method
      
      q_density_matrix <- cpp_model$get_density_vector ### Once in C++ we have computed the Voronoi representation of the functions, the density Q(p) is available. Therefore, we can evaluate the phi function here.
      cpp_model$set_phi_matrix(phi_function(q_density_matrix)) ### Apply the phi function to the empirical computed measures. the final weight will be computed inside Cpp class
      
      # Now solve the problem
      cpp_model$solve() # This action computes the integrated functional depth and all the quantities required, in Cpp
      # Store the already computed Voronoi representation of the functions to be used in predict.
      fun_representation_storage =  cpp_model$get_storage()
    },
    predict = function(f_pred, depth_types){
      # Prepare the matrices for prediction, in the same 
      f_pred_mask <- is.na(f_data)
      f_pred[f_pred_mask] <- rep(0,sum(f_pred_mask))
      
      # Set the depth types for prediction
      cpp_model$set_pred_depth_type(depth_types)
      
      output = model$predict(f_pred, f_pred_mask)
      
      print("Still to be concluded")
      }, # This function may be used to compute the depth of some new functions, w.r.t. the functions used in fit
    #domain = function(){ return(cpp_model$domain() }
    # ComputedDepths = function() { return(cpp_model$ComputedDepths()) }, # Solution, divided for univariate depth
    # ComputedHypo = function() { return(cpp_model$ComputedHypo()) }, # Modified Hypograph depth, used for outliers detection. Computed only if MHRD is required
    # ComputedEpi = function() { return(cpp_model$ComputedEpi()) } # Modified Epigraph depth, used for outliers detection. Computed only if MHRD is required
    # median = function() { return(cpp_model$median()) }, # Median, available after computation
    # FirstQuartile = function() { return(cpp_model$FirstQuartile()) }, # FirstQuartile, available after computation
    # ThirdQuartile = function() { return(cpp_model$ThirdQuartile()) }, # ThirdQuartile, available after computation
    # UpperFence = function() { return(cpp_model$UpperFence()) }, # UpperFence, available after computation
    # LowerFence = function() { return(cpp_model$LowerFence()) }, # LowerFence, available after computation
    output = function(){ print("To be implemented")} #Use the two output functions defined in the r_depth header } # This encapsulates all of the above (in a list), so that we expose a simpler routine
  )
)

# This is the generic function available to the public
setGeneric("Depth", function(f_data, locations, domain, depth_type, phi_function) { standardGeneric("Depth") })

.Depth <- function(f_data, locations, domain, depth_type, phi_function) {
### Here I need to treat the model I have (defined in a separate file) and to fill all the things that will be needed
### I have two cases: 1 data is a list of functions-locations couple, possibly missing
###                   2 data is a matrix of functions with locations separately, eventually missing
### For the moment I implement only the version with common locations, just to make it easy
  ### Evaluate the weights somehow on nodes (if possible, otherwise demand to C++, or default)
  
  ### Define the Cpp model
  ## extract local and embedding dimensions
  m <- ncol(domain$elements) - 1
  n <- ncol(domain$nodes)
  ## derive domain type
  if (m == 1 && n == 2) {
    model_ <- new(cpp_network_depth, domain)
  } else if (m == 2 && n == 2) {
    model_ <- new(cpp_2d_depth, domain)
  } else if (m == 2 && n == 3) {
    model_ <- new(cpp_surface_depth, domain)
  } else if (m == 3 && n == 3) {
    model_ <- new(cpp_3d_depth, domain)
  } else {
    stop("wrong input argument provided.")
  }
  
  ### DO SOMETHING TO REGULARIZE DATA and LOCATIONS
  # Still to be implemented
  ### Pass al the preprocessed data to the underlying C++ model
  
  # Transform the matrix with NA into a dense full matrix, with NA masl
  f_data_mask <- is.na(f_data)
  f_data[f_data_mask] <- rep(0,sum(f_data_mask))
  
  # Set the data
  model_$set_functional_data(f_data,f_data_mask) # f_data is a matrix of dimension n_stat_units x n_loc
  model_$set_locations(locations) # In the future will contain the union of the set of all the locations
  model_$set_depth_type(depth_type)
  
  domain$elements <- domain$elements - 1 ## perform index realignment for cpp handler
  model_$domain <- domain # Mesh object, used in C++
  
  if(is.null(phi_function)){
    phi_function <- function(values){return(values)} # Default identity \phi function
  }else{
    print("The phi function should be a positive function ")
  }
  
  # This function will be used afterwards after the model has been initialized
  # phi_function # functional object: needs to be a positive integrable function on Omega
  
  return(.DepthModelCtr(
    cpp_model = model_,
    phi_function = phi_function,
    fun_representation_storage = NULL
  ))
}

#' @export
setMethod(
  f = "Depth",
  signature = c(
    f_data="ANY", locations= "ANY", domain = "MeshObject", depth_type = "integer", phi_function = "missing" # function, to be set
  ),
  ## default to usual depth defined on # of missing data in node of the mesh
  definition = function(f_data, locations, domain, depth_type) {
    .Depth(f_data, locations, domain, depth_type, phi_function = NULL)
  }
)

#' @export
setMethod(
  f = "Depth",
  signature = c(
    f_data="ANY", locations= "ANY", domain = "MeshObject", depth_type = "integer", phi_function = "ANY" # function, to be set
  ),
  # In general, we may take a weight  function that is a positive function integrable on Omega
  definition = function(f_data, locations, domain, depth_type, phi_function) {
    .Depth(f_data, locations, domain, depth_type, phi_function)
  }
)