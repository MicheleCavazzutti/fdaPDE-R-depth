// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __R_DEPTH_H__
#define __R_DEPTH_H__

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/utils/symbols.h>
#include <fdaPDE/models.h> 
using fdapde::models::DEPTH;

#include "mesh.h"
using fdapde::core::Triangulation;

namespace fdapde {
  namespace r {

    template <int M, int N> class R_DEPTH { // M and N need to be specified in instantiation before
    private:
      using DomainType = fdapde::core::Triangulation<M, N>;
  
      DomainType domain_ {};           // triangulation
      DEPTH<DomainType> model_; // statistical model to wrap
    public:
      R_DEPTH(Rcpp::Environment mesh) {
	// set domain
	using RDomainType = r::Mesh<DomainType::local_dim, DomainType::embed_dim>;
	SEXP meshptr = mesh[".pointer"];
	RDomainType* ptr = reinterpret_cast<RDomainType*>(R_ExternalPtrAddr(meshptr));
	domain_ = ptr->domain();
    
	// set model instance
	model_ = DEPTH<DomainType>(domain_); // to be harmonized with the construCtor of the class that will compute the depths, called DEPTH
      }
  
      // setters
      // Here I will have the setters for the data internal to the C++ class that will compute depth. First I need to construct it, and then to harmonize it
      void set_functional_data(const DMatrix<double>& f_data, const DMatrix<bool>& f_data_mask) 
      {model_.set_train_functions(f_data);
	model_.set_train_NA_matrix(f_data_mask);} // NBB substitute the "FUNCTIONAL_DATA" with an appropriate flag in the Cpp part.
      void set_phi_function_evaluation(const DVector<double>& phi_function_evaluation) { model_.set_phi_function_evaluation(phi_function_evaluation); } // Evaluated phi matrix in R
      void set_locations(const DMatrix<double>& locations) {model_.set_locations(locations);} // NBB check that the locations cannot be directly put inside mesh
      void set_depth_types(const DVector<int>& depth_type) {model_.set_depth_types(depth_type);} // NBB substitute the vector with the approriate structure in cpp part
      void set_pred_depth_types(const DVector<int>& depth_type) {model_.set_pred_depth_types(depth_type);} // NBB substitute the vector with the approriate structure in cpp part
  
      // getters: output management
      DVector<double> get_density_vector(){ return model_.get_density_vector(); }
      //DMatrix<double> get_output_matrices(){ return model_.get_output_matrices(); }
      DMatrix<double> get_ifd(){ return model_.IFD_fit(); }
      //DMatrix<double> get_storage(){ return model_.get_storage(); }
  
      // utilities
      void init() {
	// Setting data to the C++ class
	model_.init(); // init will ask to compute all the internals, such as the Voronoy representation
      }
      void solve() { model_.solve(); } // This part will solve the model, computing the reciprocal depths.
  
      DMatrix<double> predict(const DMatrix<double> & pred_data, DMatrix<bool> & pred_mask) {
	model_.set_pred_functions(pred_data);
	model_.set_pred_NA_matrix(pred_mask);
  
	model_.predict(); 
    
	DMatrix<double> output = model_.IFD_pred(); // Here will be both the ifd, the MEI if required, the others.
	return output;
      } // This part computes the predicted IFD, MEI, ...
  
      // destructor
      ~R_DEPTH() = default;
    };

  }
}

#endif // __R_DEPTH_H__
