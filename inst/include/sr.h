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

#ifndef __R_SR_H__
#define __R_SR_H__

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include "fe_ls_elliptic.h"

namespace fdapde {
namespace r {

template <int LocalDim, int EmbedDim>
class sr_elliptic : public fe_ls_elliptic<LocalDim, EmbedDim, fdapde::SRPDE<internals::fe_ls_elliptic>> {
    using Base = fe_ls_elliptic<LocalDim, EmbedDim, fdapde::SRPDE<internals::fe_ls_elliptic>>;
   public:
    sr_elliptic() noexcept = default;
    sr_elliptic(
      const std::string& formula, const Rcpp::Environment& geoframe, const Rcpp::Nullable<Rcpp::List>& penalty) :
        Base(formula, geoframe, penalty) { }
};

}   // namespace r
}   // namespace fdapde

#endif   // __R_SR_H__
