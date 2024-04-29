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

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "headers/r_depth.h"

// Rcpp modules definition
using cpp_network_depth = R_DEPTH<1,2>;
RCPP_MODULE(cpp_network_depth) {
    Rcpp::class_<R_DEPTH<1,2>>("cpp_network_depth")
      .constructor<Rcpp::Environment, int>()
      .method("set_phi_matrix"     , &R_DEPTH<1,2>::set_phi_matrix     )
      .method("set_functional_data", &R_DEPTH<1,2>::set_functional_data)
      .method("set_locations"      , &R_DEPTH<1,2>::set_locations      )
      .method("set_depth_type"     , &R_DEPTH<1,2>::set_depth_type     )
      .method("get_density_vector" , &R_DEPTH<1,2>::get_density_vector )
      .method("get_output_matrices", &R_DEPTH<1,2>::get_output_matrices)
      .method("get_ifd"            , &R_DEPTH<1,2>::get_ifd            )
      .method("get_storage"        , &R_DEPTH<1,2>::get_storage        )
      .method("init"               , &R_DEPTH<1,2>::init               )
      .method("solve"              , &R_DEPTH<1,2>::solve              )
      .method("predict"            , &R_DEPTH<1,2>::predict            );
}

// Rcpp modules definition
using cpp_2d_depth = R_DEPTH<2,2>;
RCPP_MODULE(cpp_2d_depth) {
    Rcpp::class_<R_DEPTH<2,2>>("cpp_2d_depth")
      .constructor<Rcpp::Environment, int>()
      .method("set_phi_matrix"     , &R_DEPTH<2,2>::set_phi_matrix     )
      .method("set_functional_data", &R_DEPTH<2,2>::set_functional_data)
      .method("set_locations"      , &R_DEPTH<2,2>::set_locations      )
      .method("set_depth_type"     , &R_DEPTH<2,2>::set_depth_type     )
      .method("get_density_vector" , &R_DEPTH<2,2>::get_density_vector )
      .method("get_output_matrices", &R_DEPTH<2,2>::get_output_matrices)
      .method("get_ifd"            , &R_DEPTH<2,2>::get_ifd            )
      .method("get_storage"        , &R_DEPTH<2,2>::get_storage        )
      .method("init"               , &R_DEPTH<2,2>::init               )
      .method("solve"              , &R_DEPTH<2,2>::solve              )
      .method("predict"            , &R_DEPTH<2,2>::predict            );
}

// Rcpp modules definition
using cpp_surface_depth = R_DEPTH<2,3>;
RCPP_MODULE(cpp_surface_depth) {
    Rcpp::class_<R_DEPTH<2,3>>("cpp_surface_depth")
      .constructor<Rcpp::Environment, int>()
      .method("set_phi_matrix"     , &R_DEPTH<2,3>::set_phi_matrix     )
      .method("set_functional_data", &R_DEPTH<2,3>::set_functional_data)
      .method("set_locations"      , &R_DEPTH<2,3>::set_locations      )
      .method("set_depth_type"     , &R_DEPTH<2,3>::set_depth_type     )
      .method("get_density_vector" , &R_DEPTH<2,3>::get_density_vector )
      .method("get_output_matrices", &R_DEPTH<2,3>::get_output_matrices)
      .method("get_ifd"            , &R_DEPTH<2,3>::get_ifd            )
      .method("get_storage"        , &R_DEPTH<2,3>::get_storage        )
      .method("init"               , &R_DEPTH<2,3>::init               )
      .method("solve"              , &R_DEPTH<2,3>::solve              )
      .method("predict"            , &R_DEPTH<2,3>::predict            );
}

// Rcpp modules definition
using cpp_3d_depth = R_DEPTH<3,3>;
RCPP_MODULE(cpp_3d_depth) {
    Rcpp::class_<R_DEPTH<3,3>>("cpp_3d_depth")
      .constructor<Rcpp::Environment, int>()
      .method("set_phi_matrix"     , &R_DEPTH<3,3>::set_phi_matrix     )
      .method("set_functional_data", &R_DEPTH<3,3>::set_functional_data)
      .method("set_locations"      , &R_DEPTH<3,3>::set_locations      )
      .method("set_depth_type"     , &R_DEPTH<3,3>::set_depth_type     )
      .method("get_density_vector" , &R_DEPTH<3,3>::get_density_vector )
      .method("get_output_matrices", &R_DEPTH<3,3>::get_output_matrices)
      .method("get_ifd"            , &R_DEPTH<3,3>::get_ifd            )
      .method("get_storage"        , &R_DEPTH<3,3>::get_storage        )
      .method("init"               , &R_DEPTH<3,3>::init               )
      .method("solve"              , &R_DEPTH<3,3>::solve              )
      .method("predict"            , &R_DEPTH<3,3>::predict            );
}
