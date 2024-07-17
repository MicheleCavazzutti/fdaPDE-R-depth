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
#include "../inst/include/depth.h"

namespace fdapde {
  namespace r {

    // Rcpp modules definition
    //using cpp_network_depth = R_DEPTH<1,1>;
    //RCPP_MODULE(cpp_network_depth) {
    //    Rcpp::class_<R_DEPTH<1,1>>("cpp_linear_depth")
    //      .constructor<Rcpp::Environment>()
    //      .method("set_phi_function_evaluation"     , &R_DEPTH<1,1>::set_phi_function_evaluation)
    //      .method("set_functional_data"             , &R_DEPTH<1,1>::set_functional_data        )
    //      .method("set_locations"                   , &R_DEPTH<1,1>::set_locations              )
    //      .method("set_depth_types"                 , &R_DEPTH<1,1>::set_depth_types            )
    //	    .method("set_pred_depth_types"            , &R_DEPTH<1,1>::set_pred_depth_types       )
    //      .method("density_vector"                  , &R_DEPTH<1,1>::density_vector             )
    //      .method("ifd_fit"                         , &R_DEPTH<1,1>::ifd_fit                    )
    //      .method("ifd_pred"                        , &R_DEPTH<1,1>::ifd_pred                   )
    //	    .method("mhypo_fit"                       , &R_DEPTH<1,1>::mhypo_fit                  )
    //	    .method("mepi_fit"                        , &R_DEPTH<1,1>::mepi_fit                   )
    //	    .method("mhypo_pred"                      , &R_DEPTH<1,1>::mhypo_pred                 )
    // 	    .method("mepi_pred"                       , &R_DEPTH<1,1>::mepi_pred                  )
    //      .method("medians"                         , &R_DEPTH<1,1>::medians                    )
    //	    .method("medians_NA"                      , &R_DEPTH<1,1>::medians_NA                 )
    //      .method("first_quartile"                  , &R_DEPTH<1,1>::first_quartile             )
    //	    .method("third_quartile"                  , &R_DEPTH<1,1>::third_quartile             )
    // 	    .method("up_whisker"                      , &R_DEPTH<1,1>::up_whisker                 )
    //	    .method("low_whisker"                     , &R_DEPTH<1,1>::low_whisker                )
    // 	    .method("outliers"                        , &R_DEPTH<1,1>::outliers                   )
    //      .method("init"                            , &R_DEPTH<1,1>::init                       )
    //      .method("solve"                           , &R_DEPTH<1,1>::solve                      )
    //      .method("predict"                         , &R_DEPTH<1,1>::predict                    );
    //}

    // Rcpp modules definition
    //using cpp_network_depth = R_DEPTH<1,2>;
    //RCPP_MODULE(cpp_network_depth) {
    //    Rcpp::class_<R_DEPTH<1,2>>("cpp_network_depth")
    //      .constructor<Rcpp::Environment>()
    //      .method("set_phi_function_evaluation"     , &R_DEPTH<1,2>::set_phi_function_evaluation)
    //      .method("set_functional_data"             , &R_DEPTH<1,2>::set_functional_data        )
    //      .method("set_locations"                   , &R_DEPTH<1,2>::set_locations              )
    //      .method("set_depth_types"                 , &R_DEPTH<1,2>::set_depth_types            )
    //	    .method("set_pred_depth_types"            , &R_DEPTH<1,2>::set_pred_depth_types       )
    //      .method("density_vector"                  , &R_DEPTH<1,2>::density_vector             )
    //      .method("ifd_fit"                         , &R_DEPTH<1,2>::ifd_fit                    )
    //      .method("ifd_pred"                        , &R_DEPTH<1,2>::ifd_pred                   )
    //	    .method("mhypo_fit"                       , &R_DEPTH<1,2>::mhypo_fit                  )
    //	    .method("mepi_fit"                        , &R_DEPTH<1,2>::mepi_fit                   )
    //	    .method("mhypo_pred"                      , &R_DEPTH<1,2>::mhypo_pred                 )
    // 	    .method("mepi_pred"                       , &R_DEPTH<1,2>::mepi_pred                  )
    //      .method("medians"                         , &R_DEPTH<1,2>::medians                    )
    //	    .method("medians_NA"                      , &R_DEPTH<1,2>::medians_NA                 )
    //      .method("first_quartile"                  , &R_DEPTH<1,2>::first_quartile             )
    //	    .method("third_quartile"                  , &R_DEPTH<1,2>::third_quartile             )
    // 	    .method("up_whisker"                      , &R_DEPTH<1,2>::up_whisker                 )
    //	    .method("low_whisker"                     , &R_DEPTH<1,2>::low_whisker                )
    // 	    .method("outliers"                        , &R_DEPTH<1,2>::outliers                   )
    //      .method("init"                            , &R_DEPTH<1,2>::init                       )
    //      .method("solve"                           , &R_DEPTH<1,2>::solve                      )
    //      .method("predict"                         , &R_DEPTH<1,2>::predict                    );
    //}

    // Rcpp modules definition
    using cpp_2d_depth = R_DEPTH<2,2>;
    RCPP_MODULE(cpp_2d_depth) {
      Rcpp::class_<R_DEPTH<2,2>>("cpp_2d_depth")
	.constructor<Rcpp::Environment>()
	.method("set_phi_function_evaluation"     , &R_DEPTH<2,2>::set_phi_function_evaluation)
	.method("set_functional_data"             , &R_DEPTH<2,2>::set_functional_data        )
	.method("set_locations"                   , &R_DEPTH<2,2>::set_locations              )
	.method("set_depth_types"                 , &R_DEPTH<2,2>::set_depth_types            )
	.method("set_pred_depth_types"            , &R_DEPTH<2,2>::set_pred_depth_types       )
	.method("density_vector"                  , &R_DEPTH<2,2>::density_vector             )
	.method("ifd_fit"                         , &R_DEPTH<2,2>::ifd_fit                    )
	.method("ifd_pred"                        , &R_DEPTH<2,2>::ifd_pred                   )
	.method("mhypo_fit"                       , &R_DEPTH<2,2>::mhypo_fit                  )
	.method("mepi_fit"                        , &R_DEPTH<2,2>::mepi_fit                   )
	.method("mhypo_pred"                      , &R_DEPTH<2,2>::mhypo_pred                 )
	.method("mepi_pred"                       , &R_DEPTH<2,2>::mepi_pred                  )
	.method("medians"                         , &R_DEPTH<2,2>::medians                    )
	.method("medians_NA"                      , &R_DEPTH<2,2>::medians_NA                 )
	.method("first_quartile"                  , &R_DEPTH<2,2>::first_quartile             )
	.method("third_quartile"                  , &R_DEPTH<2,2>::third_quartile             )
	.method("up_whisker"                      , &R_DEPTH<2,2>::up_whisker                 )
	.method("low_whisker"                     , &R_DEPTH<2,2>::low_whisker                )
	.method("outliers"                        , &R_DEPTH<2,2>::outliers                   )
	.method("init"                            , &R_DEPTH<2,2>::init                       )
	.method("solve"                           , &R_DEPTH<2,2>::solve                      )
	.method("predict"                         , &R_DEPTH<2,2>::predict                    );
    }

    // Rcpp modules definition
    //using cpp_surface_depth = R_DEPTH<2,3>;
    //RCPP_MODULE(cpp_surface_depth) {
    //    Rcpp::class_<R_DEPTH<2,3>>("cpp_surface_depth")
    //      .constructor<Rcpp::Environment>()
    //      .method("set_phi_function_evaluation"     , &R_DEPTH<2,3>::set_phi_function_evaluation)
    //      .method("set_functional_data"             , &R_DEPTH<2,3>::set_functional_data        )
    //      .method("set_locations"                   , &R_DEPTH<2,3>::set_locations              )
    //      .method("set_depth_types"                 , &R_DEPTH<2,3>::set_depth_types            )
    // 	    .method("set_pred_depth_types"            , &R_DEPTH<2,3>::set_pred_depth_types       )
    //      .method("density_vector"                  , &R_DEPTH<2,3>::density_vector             )
    //      .method("ifd_fit"                         , &R_DEPTH<2,3>::ifd_fit                    )
    //      .method("ifd_pred"                        , &R_DEPTH<2,3>::ifd_pred                   )
    //	    .method("mhypo_fit"                       , &R_DEPTH<2,3>::mhypo_fit                  )
    //	    .method("mepi_fit"                        , &R_DEPTH<2,3>::mepi_fit                   )
    //	    .method("mhypo_pred"                      , &R_DEPTH<2,3>::mhypo_pred                 )
    // 	    .method("mepi_pred"                       , &R_DEPTH<2,3>::mepi_pred                  )
    //      .method("medians"                         , &R_DEPTH<2,3>::medians                    )
    //	    .method("medians_NA"                      , &R_DEPTH<2,3>::medians_NA                 )
    //      .method("first_quartile"                  , &R_DEPTH<2,3>::first_quartile             )
    //	    .method("third_quartile"                  , &R_DEPTH<2,3>::third_quartile             )
    // 	    .method("up_whisker"                      , &R_DEPTH<2,3>::up_whisker                 )
    //	    .method("low_whisker"                     , &R_DEPTH<2,3>::low_whisker                )
    // 	    .method("outliers"                        , &R_DEPTH<2,3>::outliers                   )
    //      .method("init"                            , &R_DEPTH<2,3>::init                       )
    //      .method("solve"                           , &R_DEPTH<2,3>::solve                      )
    //      .method("predict"                         , &R_DEPTH<2,3>::predict                    );
    //}

    // Rcpp modules definition
    using cpp_3d_depth = R_DEPTH<3,3>;
    RCPP_MODULE(cpp_3d_depth) {
        Rcpp::class_<R_DEPTH<3,3>>("cpp_3d_depth")
          .constructor<Rcpp::Environment>()
          .method("set_phi_function_evaluation"     , &R_DEPTH<3,3>::set_phi_function_evaluation)
          .method("set_functional_data"             , &R_DEPTH<3,3>::set_functional_data        ) 
          .method("set_locations"                   , &R_DEPTH<3,3>::set_locations              )
          .method("set_depth_types"                 , &R_DEPTH<3,3>::set_depth_types            )
          .method("set_pred_depth_types"            , &R_DEPTH<3,3>::set_pred_depth_types       )
          .method("density_vector"                  , &R_DEPTH<3,3>::density_vector             )
          .method("ifd_fit"                         , &R_DEPTH<3,3>::ifd_fit                    )
          .method("ifd_pred"                        , &R_DEPTH<3,3>::ifd_pred                   )
          .method("mhypo_fit"                       , &R_DEPTH<3,3>::mhypo_fit                  )
          .method("mepi_fit"                        , &R_DEPTH<3,3>::mepi_fit                   )
          .method("mhypo_pred"                      , &R_DEPTH<3,3>::mhypo_pred                 )
          .method("mepi_pred"                       , &R_DEPTH<3,3>::mepi_pred                  )
          .method("medians"                         , &R_DEPTH<3,3>::medians                    )
          .method("medians_NA"                      , &R_DEPTH<3,3>::medians_NA                 )
          .method("first_quartile"                  , &R_DEPTH<3,3>::first_quartile             )
          .method("third_quartile"                  , &R_DEPTH<3,3>::third_quartile             )
          .method("up_whisker"                      , &R_DEPTH<3,3>::up_whisker                 )
          .method("low_whisker"                     , &R_DEPTH<3,3>::low_whisker                )
          .method("outliers"                        , &R_DEPTH<3,3>::outliers                   )
          .method("init"                            , &R_DEPTH<3,3>::init                       )
          .method("solve"                           , &R_DEPTH<3,3>::solve                      )
          .method("predict"                         , &R_DEPTH<3,3>::predict                    );
    }

  }
}
