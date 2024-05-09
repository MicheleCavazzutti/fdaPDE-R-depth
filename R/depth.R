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

.DepthModel <- R6::R6Class(
  "DepthModel",
  private = list(
    model_  = NULL, ## C++ model backend
    phi_function_ = NULL ## This phi function will be used to compute the weight function, when needed. If left to NULL,the identity function will be employed
  ),
  public = list(
    initialize = function(domain, f_data, f_data_mask, locations, depth_type, phi_function) { 
      ### Define the C++ model
      ## extract local and embedding dimensions
      m <- ncol(domain$elements) - 1
      n <- ncol(domain$nodes)
      ## derive domain type
      if (m == 1 && n == 2) {
        # Deactivate due to Triangulation limitations # model_ <- new(cpp_network_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 2) {
        private$model_ <- new(cpp_2d_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 3) {
        # Deactivate due to Triangulation limitations # model_ <- new(cpp_surface_depth, get_private(domain)$mesh_)
      } else if (m == 3 && n == 3) {
        # Deactivate due to Triangulation limitations # model_ <- new(cpp_3d_depth, get_private(domain)$mesh_)
      } else {
        stop("wrong input argument provided.")
      }
      
      ### Set the data inside the C++ model
      private$model_$set_functional_data(f_data,f_data_mask) # f_data is a matrix of dimension n_stat_units x n_loc
      private$model_$set_locations(locations) # In the future will contain the union of the set of all the locations
      private$model_$set_depth_types(depth_type)
      
      # Set the C++ model and the phi_function to be evaluates
      private$phi_function_ = phi_function
    },
    init = function(){
      ### Initialization of the model: computing voronoi tessellation ad allowing for phi-function evaluation
      private$model_$init() # Needed to pass the phi-function evaluation to C++ class.
      
      # We extract from C++ the coverage density Q(p) (probability of a function to be observed in a Voronoi cell)
      q_density_vector <- private$model_$get_density_vector() 
      
      # Apply the phi function to the empirical computed measures. The final weight will be computed inside C++ class
      private$model_$set_phi_function_evaluation(private$phi_function_(q_density_vector)) 
    },
    fit = function() { 
      # Solve the problem : this action computes the integrated functional depth and all the quantities required, in C++, but does not return them
      private$model_$solve() 
    },
    predict = function(f_pred, depth_types){
      # Prepare the matrices for prediction, in the same fashion of the fit data
      f_pred_mask <- is.na(f_pred)
      f_pred[f_pred_mask] <- rep(0,sum(f_pred_mask))
      
      # Set the depth types for prediction
      private$model_$set_pred_depth_types(depth_types)
      
      private$model_$predict(f_pred, f_pred_mask)
    }, # This function may be used to compute the depth of some new functions, w.r.t. the functions used in fit
    phi_function = function(value){
      return(private$phi_function_(value))
    },
    solve = function(){
      private$model_$solve()
    },
    # domain = function(){ return(cpp_model$domain() }
    # ComputedHypo = function() { return(cpp_model$ComputedHypo()) }, # Modified Hypograph depth, used for outliers detection. Computed only if MHRD is required
    # ComputedEpi = function() { return(cpp_model$ComputedEpi()) } # Modified Epigraph depth, used for outliers detection. Computed only if MHRD is required
    # median = function() { return(cpp_model$median()) }, # Median, available after computation
    # FirstQuartile = function() { return(cpp_model$FirstQuartile()) }, # FirstQuartile, available after computation
    # ThirdQuartile = function() { return(cpp_model$ThirdQuartile()) }, # ThirdQuartile, available after computation
    # UpperFence = function() { return(cpp_model$UpperFence()) }, # UpperFence, available after computation
    # LowerFence = function() { return(cpp_model$LowerFence()) }, # LowerFence, available after computation
    IFD_fit = function(){ 
      # For the moment, the output just contains the evaluation of the IFD of fit functions
      return(private$model_$ifd_fit())
    }
    IFD_pred = function(){ 
      # For the moment, the output just contains the evaluation of the IFD of fit functions
      return(private$model_$ifd_pred())
    }
  )
)
 
# Public interface
#' @export
Depth <- function(f_data, locations, domain, depth_type, phi_function = NULL){
  ### Here I need to treat the model I have (defined in a separate file) and to fill all the things that will be needed
  ### I have two cases: 1 data is a list of functions-locations couple, possibly missing
  ###                   2 data is a matrix of functions with locations separately, eventually missing
  ### For the moment I implement only the version with common locations, just to make it easy
  ### Evaluate the weights somehow on nodes (if possible, otherwise demand to C++, or default)
  
  # Transform the matrix with NA into a dense full matrix, with NA mask
  f_data_mask <- is.na(f_data)
  f_data[f_data_mask] <- rep(0,sum(f_data_mask))
  
  # This function will be used after the model has been initialized
  # phi_function # functional object: needs to be a positive integrable function on Omega
  if(is.null(phi_function)){
    phi_function <- function(values){return(values)} # Default identity \phi function
  }else{
    print("The phi function should be a positive function \n")
  }
  
  print("At the moment, few checks on the parameters are provided. In future we will add further, detailed, checks.")
  
  # Build the R class, return it
  model = .DepthModel$new(domain, f_data, f_data_mask, locations, depth_type, phi_function)
  return(model)
}