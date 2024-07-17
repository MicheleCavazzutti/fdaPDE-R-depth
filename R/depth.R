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
    phi_function_ = NULL, ## This phi function will be used to compute the weight function, when needed. If left to NULL,the identity function will be employed
    initialized_ = FALSE, # Flag for model initialization
    solved_ = FALSE, # Flag for IFD computation for fit function
    predicted_ = FALSE, # Flag for IFD computation for pred function
    MHRD_fit_computed_ = FALSE, # This flag is need to understand whether the MHRD for fit is computed or not. In the latter case, mepi and mhypo are not available
    MHRD_pred_computed_ = FALSE # This flag is need to understand whether the MHRD for pred is computed or not. In the latter case, mepi and mhypo are not available
  ),
  public = list(
    initialize = function(domain, f_data, f_data_mask, locations, depth_types, phi_function) { 
      ### Define the C++ model
      ## extract local and embedding dimensions
      m <- ncol(domain$elements) - 1
      n <- ncol(domain$nodes)
      ## derive domain type
      if (m == 1 && n == 1) {
        # In future will be available, voronoi developed # private$model_ <- new(cpp_linear_depth, get_private(domain)$mesh_)
      } else if (m == 1 && n == 2) {
        # Deactivate due to Triangulation limitations # private$model_ <- new(cpp_network_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 2) {
        private$model_ <- new(cpp_2d_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 3) {
        # Deactivate due to Triangulation limitations # private$model_ <- new(cpp_surface_depth, get_private(domain)$mesh_)
      } else if (m == 3 && n == 3) {
        # Under test
        private$model_ <- new(cpp_3d_depth, get_private(domain)$mesh_)
      } else {
        stop("wrong input argument provided.")
      }
      
      ### Transform depth_types into a numeric type
      depth_types_num = NULL
      for(i in 1:length(depth_types)){
        num=0
        if(depth_types[i]=='MBD'){
          depth_types_num = c(depth_types_num,1) 
        }else{
          if(depth_types[i]=='FMD'){
            depth_types_num = c(depth_types_num,2)
          }else{
            if(depth_types[i]=='MHRD'){
              depth_types_num = c(depth_types_num,3)
              private$MHRD_fit_computed_ = TRUE
            }else{
              stop("Depth type should be a vector containing strings among 'MBD', 'FMD', 'MHRD'" )
            }
          }
        }
      }
      
      ### Set the data inside the C++ model
      private$model_$set_functional_data(f_data,f_data_mask) # f_data is a matrix of dimension n_stat_units x n_loc
      private$model_$set_locations(locations) # In the future will contain the union of the set of all the locations
      private$model_$set_depth_types(depth_types_num)
      
      # Set the C++ model and the phi_function to be evaluates
      private$phi_function_ = phi_function
    },
    init = function(){
      ### Initialization of the model: computing voronoi tessellation ad allowing for phi-function evaluation
      private$model_$init() # Needed to pass the phi-function evaluation to C++ class.
      
      # We extract from C++ the coverage density Q(p) (probability of a function to be observed in a Voronoi cell)
      q_density_vector <- private$model_$density_vector() 
      
      # Apply the phi function to the empirical computed measures. The final weight will be computed inside C++ class
      private$model_$set_phi_function_evaluation(private$phi_function_(q_density_vector)) 
      
      # Set the proper flag
      private$initialized_ = TRUE
    },
    solve = function(){
      if(!private$initialized_){
        stop("The model has not been initialized - run init()")
      }
      private$model_$solve()
      
      # Set the proper flag
      private$solved_ = TRUE
    },
    predict = function(f_pred, depth_types){
      
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      # Prepare the matrices for prediction, in the same fashion of the fit data
      f_pred_mask <- is.na(f_pred)
      f_pred[f_pred_mask] <- rep(0,sum(f_pred_mask))
      
      ### Transform depth_types into a numeric type
      depth_types_num = NULL
      for(i in 1:length(depth_types)){
        num=0
        if(depth_types[i]=='MBD'){
          depth_types_num = c(depth_types_num,1) 
        }else{
          if(depth_types[i]=='FMD'){
            depth_types_num = c(depth_types_num,2)
          }else{
            if(depth_types[i]=='MHRD'){
              depth_types_num = c(depth_types_num,3)
              private$MHRD_pred_computed_ = TRUE
            }else{
              stop("Depth type should be a vector containing strings among 'MBD', 'FMD', 'MHRD'" )
            }
          }
        }
      }
      
      # Set the depth types for prediction
      private$model_$set_pred_depth_types(depth_types)
      
      private$model_$predict(f_pred, f_pred_mask)
      
      # Set the proper flag
      private$predicted_ = TRUE
    }, # This function may be used to compute the depth of some new functions, w.r.t. the functions used in fit
    phi_function = function(value){
      return(private$phi_function_(value))
    },
    IFD_fit = function(){ 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      # For the moment, the output just contains the evaluation of the IFD of fit functions
      return(private$model_$ifd_fit())
    },
    IFD_pred = function(){ 
      if(!private$predicted_){
        stop("No predicted depeths available - run predict(...)")
      }
      
      # For the moment, the output just contains the evaluation of the IFD of fit functions
      return(private$model_$ifd_pred())
    },
    mhypo_fit = function() { 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      if(!private$MHRD_fit_computed_){
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for fit functions")
      }
      
      return(private$model_$mhypo_fit())
    },   # Modified Hypograph depth for fit function. Computed only if MHRD for fit is required
    mepi_fit = function() { 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      if(!private$MHRD_fit_computed_){
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for fit functions")
      }
      
      return(private$model_$mepi_fit())
    },     # Modified Epigraph depth for fit functions. Computed only if MHRD for fit is required
    mhypo_pred = function() { 
      if(!private$predicted_){
        stop("No predicted depeths available - run predict(...)")
      }
      
      if(!private$MHRD_pred_computed_){
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for fit functions")
      }
      
      return(private$model_$mhypo_pred())
    }, # Modified Hypograph depth for pred function. Computed only if MHRD for pred is required
    mepi_pred = function() {
      if(!private$predicted_){
        stop("No predicted depeths available - run predict(...)")
      }
      
      if(!private$MHRD_pred_computed_){
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for fit functions")
      }
      
      return(private$model_$mepi_pred()) 
    },    # Modified Epigraph depth for pred functions. Computed only if MHRD for  pred is required
    medians = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      medians = private$model_$medians()
      medians_mask = private$model_$medians_NA()
      
      medians[medians_mask]<-rep(NA,medians_mask) # Put to NA the missing values of the original functions
      
      return(medians)
    }, # Median, available after computation
    FirstQuartile = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      return(private$model_$first_quartile()) 
    }, # FirstQuartile, available after computation
    ThirdQuartile = function() { 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
    return(private$model_$third_quartile()) 
    }, # ThirdQuartile, available after computation
    UpperFence = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
    return(private$model_$up_whisker()) 
    }, # UpperFence, available after computation
    LowerFence = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
    return(private$model_$low_whisker()) 
    }, # LowerFence, available after computation
    outliers = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      return(private$model_$outliers())
    }
  )
)
 
# Public interface
#' @export
Depth <- function(f_data, locations, domain, depth_types, phi_function = NULL){
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
    warning("The phi function should be a positive function \n")
  }
  
  # Build the R class, return it
  model = .DepthModel$new(domain, f_data, f_data_mask, locations, depth_types, phi_function)
  return(model)
}