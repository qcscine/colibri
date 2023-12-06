/**
 * @file Tensor.h
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group.\n See LICENSE.txt for details.
 */

#ifndef VIB_TENSOR_H
#define VIB_TENSOR_H

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <unordered_set>

template<int rank>
using tensor = std::unordered_map<std::array<int, rank>, double, boost::hash<std::array<int, rank>>>;

// tensorMC stores all coeffs for a specific mode-coupling order, regardless of
// the rank for example the coeff for $Q_3^5$ would be a one-mode term stored
// with {(3,5)}
template<int mcorder>
using tensorMC =
    std::unordered_map<std::array<std::pair<int, int>, mcorder>, double, boost::hash<std::array<std::pair<int, int>, mcorder>>>;

// tensorMCSP stores all terms for a specific mode-coupling order, regardless of
// the rank, in a more efficient manner for example the coeff for $Q_3^5 Q_7^2$
// would be a one-mode term stored with {3,7} with the value ({5, 2}, val)
template<int mcorder>
using tensorMCSP = std::unordered_map<std::array<int, mcorder>, std::vector<std::pair<std::array<int, mcorder>, double>>,
                                      boost::hash<std::array<int, mcorder>>>;

// setNbody stores all keys of n modes that are coupled in n-th order mode
// coupling
template<int mcorder>
using setNbody = std::unordered_set<std::array<int, mcorder>, boost::hash<std::array<int, mcorder>>>;

template<int mcorder>
using calculatedSPs =
    std::unordered_map<std::array<std::pair<int, int>, mcorder>, double, boost::hash<std::array<std::pair<int, int>, mcorder>>>;

#endif
