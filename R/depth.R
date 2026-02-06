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
    solved_ = FALSE, # Flag for depth computation for fit functions
    predicted_ = FALSE, # Flag for IFD computation for pred function
    MHRD_fit_computed_ = FALSE, # This flag is need to understand whether the MHRD for fit is computed or not. In the latter case, mepi and mhypo are not available
    MHRD_pred_computed_ = FALSE # This flag is need to understand whether the MHRD for pred is computed or not. In the latter case, mepi and mhypo are not available
  ),
  public = list(
    initialize = function(domain, f_data_list, f_data_mask_list, locations_list, depth_types, int_method, phi_function, region_of_interest, external_measures_vector) { 
      ### Define the C++ model
      ## extract local and embedding dimensions
      m <- ncol(domain$elements) - 1
      n <- ncol(domain$nodes)
      ## derive domain type
      if (m == 1 && n == 1) {
        # Not available yet # private$model_ <- new(cpp_linear_depth, get_private(domain)$mesh_)
      } else if (m == 1 && n == 2) {
        # Not available yet # private$model_ <- new(cpp_network_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 2) {
        private$model_ <- new(cpp_2d_depth, get_private(domain)$mesh_)
      } else if (m == 2 && n == 3) {
        private$model_ <- new(cpp_surface_depth, get_private(domain)$mesh_) # Remark: when Voronoi depth is required, external voronoi measures are ued in the integration
      } else if (m == 3 && n == 3) {
        private$model_ <- new(cpp_3d_depth, get_private(domain)$mesh_) # Remark: when Voronoi depth is required, external voronoi measures are ued in the integration
      } else {
        stop("wrong input argument provided.")
      }
      
      ### Transform depth_types into a numeric type
      depth_types_num <- numeric(length(depth_types))
      for (i in seq_along(depth_types)) {
        depth_types_num[i] <- switch(depth_types[i],
                                     "SD"   = 1,
                                     "FMD"  = 2,
                                     "MHRD" = {
                                       private$MHRD_fit_computed_ <- TRUE
                                       3
                                     },
                                     "DI-SD" = 4,
                                     "PDI-SD" = 5,
                                     # Default case
                                     stop("Depth type should be a vector containing strings among 'SD', 'FMD', 'MHRD', 'DI-SD', 'PDI-SD'.")
        )
      }
      
      ### Transform int_method in numeric value (0 --> "FEM-0", 1 --> "Voronoi")
      int_method_num <- ifelse(int_method=="Voronoi",-1,0)
      
      ### Set the data inside the C++ model
      private$model_$set_functional_data(f_data_list,f_data_mask_list) # f_data_list is a list of length n_train, where each element is a vector (possibly of different size). 
      private$model_$set_locations(locations_list) # Two cases are possible here: either a list of length 1, containing a single matrix of locations common to every functional datum, or a list of length n_train, with a locations matrix for each functional datum
      private$model_$set_depth_types(depth_types_num) # List of numbers indicating the depth types required
      private$model_$set_int_method(int_method_num) # integer indicating the type of integration required. Currently can be of value -1 (Voronoi depth) or 0 (FEM-0 depth)
      private$model_$set_roi(region_of_interest) # Two cases are possible here: a vector of dimension 1 with value -1 or a vector of arbitrary dimension with non-negative integre values
      
      # Just for 2.5D and 3D cases
      if (((m == 2 && n == 3) || (m==3 && n == 3)) && (int_method_num == -1) ){ if(length(external_measures_vector) == 1){stop("You need to provide external measures in 2.5D and 3D cases, when Voronoi integration is required")}}
      private$model_$set_external_voronoi_measures(external_measures_vector) # Has meaning only in the case of 2.5D 
      
      # Set the C++ model and the phi_function to be evaluates
      private$phi_function_ = phi_function
    },
    init = function(){
      ### Initialization of the model: set the data and compute the seed-based representation for the data (that is computing the Voronoi tessellation and computing the spatial averages in the Voronoi cell if Voronoi integration is required). 
      private$model_$init()
      
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
    predict = function(f_pred, locations_pred, depth_types){ ## Similarly to the fit case, accepts a list of vectors (representing the pred functions) and the associated list of locations.
      
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      ### Get from model the type of integration used to check the parameters
      int_method_num = private$model_$int_method()
      
      ### Check that functional data are either matrix or list of vectors
      if(!(is.matrix(f_pred) || (is.list(f_pred) && all(vapply(f_pred, is.numeric, TRUE))))){
        stop("f_pred must be either a matrix or a list of numeric vectors.")
      }
      
      if(!(is.null(locations_pred) ||  is.matrix(locations_pred) || (is.list(locations_pred) && all(vapply(locations_pred, is.matrix, TRUE))))){ 
        # Note: I am exploiting the fact that if locations_pred is null I'm not evaluating the other terms
        stop("locations_pred must be NULL, a matrix, or a list of matrices.")
      }
      
      if(is.list(f_pred) && !is.list(locations_pred)){
        stop("If f_pred is a list of vectors, locations must be a list of matrices. If all the functions in f_pred share the same locations, store them in a matrix.")
      }
      if(is.list(locations_pred) && !is.list(locations_pred)){
        stop("If locations is a list of matrices, f_pred must be a list of vectors.")
      }
      
      if(is.list(f_pred) && is.list(locations_pred)){
        if(length(f_pred) != length(locations_pred)){
          stop("f_pred and locations_pred must have the same number of elements when both are lists.")
        }
        for(i in seq_along(f_pred)){
          if(length(f_pred[[i]]) != nrow(locations_pred[[i]])){
            stop(sprintf(
              "Length mismatch: f_pred[[%d]] has length %d but locations_pred[[%d]] has %d rows.",
              i, length(f_pred[[i]]), i, nrow(locations_pred[[i]])
            ))
          }
        }
      }
      
      if(int_method_num == -1 && is.null(locations_pred)){
        stop("locations_pred must be provided (non-NULL) when int_method='Voronoi'.")
      }
      if(int_method_num == 0 && !is.null(locations_pred)){
        stop("locations_pred must be NULL when int_method='FEM-0'.")
      }
      
      if(is.matrix(f_pred) && is.matrix(locations_pred)){
        if(ncol(f_pred) != nrow(locations_pred)){
          stop(sprintf(
            "Column mismatch: f_pred has %d columns but locations_pred has %d rows.",
            ncol(f_pred), nrow(locations_pred)
          ))
        }
      }
      
      ### f_pred representation: represent f_pred with two lists, where in the first we store the values and in the second the NA_Masks
      if(is.matrix(f_pred)){
        f_pred_list      <- vector("list", nrow(f_pred))
        f_pred_mask_list <- vector("list", nrow(f_pred))
        
        for(i in 1:nrow(f_pred)){
          f_i <- f_pred[i,]
          mask_i <- is.na(f_i)
          f_i[mask_i] <- rep(0,sum(mask_i))  # replace NA with 0
          
          f_pred_list[[i]]      <- f_i
          f_pred_mask_list[[i]] <- mask_i
        }
      }else{
        f_pred_list      <- vector("list", length(f_pred))
        f_pred_mask_list <- vector("list", length(f_pred))
        
        for(i in 1:length(f_pred)){
          f_i <- f_pred[[i]]
          mask_i <- is.na(f_i)
          f_i[mask_i] <- rep(0,sum(mask_i))
          
          f_pred_list[[i]]      <- f_i
          f_pred_mask_list[[i]] <- mask_i
        }
      }
      
      ### Transform locations_pred into the standard list representation. If only one element is in the list, we need to compute Voronoi areas just once (if needed)
      if(is.null(locations_pred)){
        ### Locations is just a list with one element, the mesh nodes
        locations_list <- list(matrix(0,nrow=1,ncol=1)) # Fake locations list for FEM integration, ignored in C++
      } else if(is.matrix(locations_pred)){
        ### Locations is just a list with one element, the original locations_pred
        locations_list <- list(locations_pred)
      } else if(is.list(locations_pred)){
        locations_list <- locations_pred
      } else {
        stop("locations_pred must be NULL, a matrix, or a list of matrices.")
      }
      
      ### Transform depth_types into a numeric type
      depth_types_num <- numeric(length(depth_types))
      
      for (i in seq_along(depth_types)) {
        depth_types_num[i] <- switch(depth_types[i],
                                     "SD"   = 1,
                                     "FMD"  = 2,
                                     "MHRD" = {
                                       # Imposta il flag specifico per la predizione
                                       private$MHRD_pred_computed_ <- TRUE
                                       3
                                     },
                                     "DI-SD" = {
                                       stop("Double integral not implemented yet in predict")
                                     },
                                     "PDI-SD" = {
                                       stop("Partial Double integral not implemented yet in predict")
                                     },
                                     # Caso di errore se la stringa non è tra quelle permesse
                                     stop("Depth type should be a vector containing strings among 'SD', 'FMD', 'MHRD'")
        )
      }
      
      # Set the depth types for prediction
      private$model_$set_pred_depth_types(depth_types_num)
      
      private$model_$predict(f_pred_list, f_pred_mask_list, locations_list)
      
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
    f_fit_representaiton = function(){ 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      f_fit = private$model_$f_fit()
      f_fit_mask = private$model_$f_fit_NA()
      
      f_fit[f_fit_mask]<-rep(NA,sum(f_fit_mask)) # Put to NA the missing values
      
      return(f_fit)
    },
    f_pred_representaiton = function(){ 
      if(!private$predicted_){
        stop("No predicted depeths available - run predict(...)")
      }
      
      f_pred = private$model_$f_pred()
      f_pred_mask = private$model_$f_pred_NA()
      
      f_pred[f_pred_mask]<-rep(NA,sum(f_pred_mask)) # Put to NA the missing values
      
      return(f_pred)
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
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for pred functions")
      }
      
      return(private$model_$mhypo_pred())
    }, # Modified Hypograph depth for pred function. Computed only if MHRD for pred is required
    mepi_pred = function() {
      if(!private$predicted_){
        stop("No predicted depeths available - run predict(...)")
      }
      
      if(!private$MHRD_pred_computed_){
        stop("Epigraph and Hypograph indexes are available only if MHRD has been computed for pred functions")
      }
      
      return(private$model_$mepi_pred()) 
    },    # Modified Epigraph depth for pred functions. Computed only if MHRD for  pred is required
    medians = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      medians = private$model_$medians()
      medians_mask = private$model_$medians_NA()
      
      medians[medians_mask]<-rep(NA,sum(medians_mask)) # Put to NA the missing values of the original functions
      
      return(medians)
    }, # Median, available after computation
    FirstQuartile = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      first_quartile = private$model_$first_quartile()
      first_quartile_mask = private$model_$first_quartile_NA()
      
      first_quartile[first_quartile_mask]<-rep(NA,sum(first_quartile_mask)) # Put to NA the missing values

      return(first_quartile) 
    }, # FirstQuartile, available after computation
    ThirdQuartile = function() { 
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      third_quartile = private$model_$third_quartile()
      third_quartile_mask = private$model_$third_quartile_NA()
      
      third_quartile[third_quartile_mask]<-rep(NA,sum(third_quartile_mask)) # Put to NA the missing values
      
      return(third_quartile) 
    }, # ThirdQuartile, available after computation
    UpperFence = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      up_whisker = private$model_$up_whisker()
      up_whisker_mask = private$model_$up_whisker_NA()
      
      up_whisker[up_whisker_mask]<-rep(NA,sum(up_whisker_mask)) # Put to NA the missing values
      
      return(up_whisker) 
    }, # UpperFence, available after computation
    LowerFence = function() {
      if(!private$solved_){
        stop("The model has not been solved - run solve()")
      }
      
      low_whisker = private$model_$low_whisker()
      low_whisker_mask = private$model_$low_whisker_NA()
      
      low_whisker[low_whisker_mask]<-rep(NA,sum(low_whisker_mask)) # Put to NA the missing values
      
      return(low_whisker) 
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
Depth <- function(f_data, locations = NULL, domain, depth_types, int_method = 'Voronoi', phi_function = NULL, region_of_interest = NULL, external_measures_vector = NULL){
  ### Description of the inputs:
  ### - f_data: evaluations of the functional data in the locations. Can be two things: a matrix (in case only one set of locations is available
  ### or locations is NULL, that is FEM case) or a list() of vectors (in the case one wants to specify different locations for each functional datum
  ### available only in the Voronoi case)
  ### - locations: set of locations for the functional data. May be NULL (in this case the locations coincide with mesh nodes), 
  ### may be a single matrix (available only in the case of Voronoi, than all the functional data need to have the same length),
  ### may be a list of matrices (one for each functional datum, need to have the same length)
  ### - domain: mesh representing the problem
  ### - depth_types: a vector specifying the types of depths one wants to compute on the provided data. Computing different univariate depths does not bring any overhead.
  ### - int_method: can take value "Voronoi" or "FEM-0", indicated the type of integration one wants to perform
  ### - phi_function: type of weight function one would like to use in the depth integral weight
  ### - region_of_interest: set of elements that compose the area of interest w.r.t. which the Partial Double Integral (PDI-) depths are computed.
  ### If left to NULL, a surrounding area for each node is selected.
  ### - external_measures_vector: vector of length (number of nodes) that is reporting the Voronoi areas associated to the mesh nodes. Needed only in Voronoi integration and 2.5/3 dimensional problems.
  
  ### Check that functional data are either matrix or list of vectors
  if(!(is.matrix(f_data) || (is.list(f_data) && all(vapply(f_data, is.numeric, TRUE))))){
    stop("f_data must be either a matrix or a list of numeric vectors.")
  }
  
  if(!(is.null(locations) ||  is.matrix(locations) || (is.list(locations) && all(vapply(locations, is.matrix, TRUE))))){ 
    # Note: I am exploiting the fact that if locations is null I'm not evaluating the other terms
    stop("locations must be NULL, a matrix, or a list of matrices.")
  }

  if(!(int_method %in% c("Voronoi", "FEM-0"))){
    stop('int_method must be either "Voronoi" or "FEM-0".')
  }
  
  if(is.list(f_data) && !is.list(locations)){
    stop("If f_data is a list of vectors, locations must be a list of matrices.")
  }
  if(is.list(locations) && !is.list(locations)){
    stop("If locations is a list of matrices, f_data must be a list of vectors.")
  }
  
  if(is.list(f_data) && is.list(locations)){
    if(length(f_data) != length(locations)){
      stop("f_data and locations must have the same number of elements when both are lists.")
    }
    for(i in seq_along(f_data)){
      if(length(f_data[[i]]) != nrow(locations[[i]])){
        stop(sprintf(
          "Length mismatch: f_data[[%d]] has length %d but locations[[%d]] has %d rows.",
          i, length(f_data[[i]]), i, nrow(locations[[i]])
        ))
      }
    }
  }
  
  if(int_method == "Voronoi" && is.null(locations)){
    stop("locations must be provided (non-NULL) when int_method='Voronoi'.")
  }
  if(int_method == "FEM-0" && !is.null(locations)){
    stop("locations must be NULL when int_method='FEM-0'.")
  }
  
  if(is.matrix(f_data) && is.matrix(locations)){
    if(ncol(f_data) != nrow(locations)){
      stop(sprintf(
        "Column mismatch: f_data has %d columns but locations has %d rows.",
        ncol(f_data), nrow(locations)
      ))
    }
  }
  
  ### f_data representation: represent f_data with two lists, where in the first we store the values and in the second the NA_Masks
  if(is.matrix(f_data)){
    f_data_list      <- vector("list", nrow(f_data))
    f_data_mask_list <- vector("list", nrow(f_data))
    
    for(i in 1:nrow(f_data)){
      f_i <- f_data[i,]
      mask_i <- is.na(f_i)
      f_i[mask_i] <- rep(0,sum(mask_i))  # replace NA with 0
      
      f_data_list[[i]]      <- f_i
      f_data_mask_list[[i]] <- mask_i
    }
  }else{
    f_data_list      <- vector("list", length(f_data))
    f_data_mask_list <- vector("list", length(f_data))
    
    for(i in 1:length(f_data)){
      f_i <- f_data[[i]]
      mask_i <- is.na(f_i)
      f_i[mask_i] <- rep(0,sum(mask_i))
      
      f_data_list[[i]]      <- f_i
      f_data_mask_list[[i]] <- mask_i
    }
  }
  
  ### Transform locations into the standard list representation. If only one element is in the list, we need to compute Voronoi areas just once (if needed)
  if(is.null(locations)){
    ### Locations is just a list with one element, the mesh nodes
    if(is.null(domain$nodes) || !is.matrix(domain$nodes)){
      stop("domain$nodes must be a matrix when locations is NULL.")
    }
    locations_list <- list(domain$nodes)
  } else if(is.matrix(locations)){
    ### Locations is just a list with one element, the original locations
    locations_list <- list(locations)
  } else if(is.list(locations)){
    locations_list <- locations
  } else {
    stop("locations must be NULL, a matrix, or a list of matrices.")
  }
  
  # This function will be used after the model has been initialized
  # phi_function # functional object: needs to be a positive integrable function on Omega
  if(is.null(phi_function)){
    phi_function <- function(values){return(values)} # Default identity \phi function
  }else{
    warning("The phi function should be a positive function \n")
  }
  
  ### Check well posedness of the ROI, if any
  if (is.null(region_of_interest)) { # Set default value, that implies the contruction of a surrounding patch for each node
    region_of_interest <- 0
  } else {
    if (!is.numeric(region_of_interest) || any(region_of_interest %% 1 != 0)) {
      stop("region_of_interest must be a vector of positive integers if not left to NULL")
    }
    if (any(region_of_interest < 1)) {
      stop("region_of_interest must be a vector of positive integers if not left to NULL")
    }
    if (int_method == "FEM-0") {
      max_val <- nrow(domain$elements()) # In region of interest should be specified the values of the possible elements in the ROI
      if (any(region_of_interest > max_val)) {
        stop(
          "region_of_interest values should range from 1 to nrow(domain$elements()) if FEM-0 integration is required"
        ))
      }
    } else if (int_method == "Voronoi") {
      max_val <- nrow(domain$nodes()) # In region of interest should be specified the values of the nodes elements in the ROI
      if (any(region_of_interest > max_val)) {
        stop(
          "region_of_interest values should range from 1 to nrow(domain$nodes()) if Voronoi integration is required"
        ))
      }
    } else {
      stop("int_method non riconosciuto.")
    }
  }
  
  ### Cpp index alignment
  region_of_interest <- as.integer(region_of_interest) -1 # Shift the indexes to match C++ notation
  
  # Just for 2.5D objects and 3D objects
  if(is.null(external_measures_vector)){
    external_measures_vector <- as.numeric(rep(0,1)) # Default useless external measure vector
  }
  
  # Build the R class, return it
  model = .DepthModel$new(domain, f_data_list, f_data_mask_list, locations_list, depth_types, int_method, phi_function, region_of_interest, external_measures_vector)
  return(model)
}