/*
 * SourceVariable.cpp
 *
 *  Created on: 13 Nov 2025
 *      Author: ntc132
 */

#include "SourceVariable.h"


double SourceVariable::getLuminosity(double time) {

	if (time <= times_.front()) return luminosities_.front();
    if (time >= times_.back())  return luminosities_.back();

    auto it = std::upper_bound(times_.begin(), times_.end(), time);
    size_t i = static_cast<size_t>((it - times_.begin()) - 1);

    double t0 = times_[i];
    double t1 = times_[i+1];
    double L0 = luminosities_[i];
    double L1 = luminosities_[i+1];

    double alpha = (time - t0) / (t1 - t0);
    double L = L0 + alpha * (L1 - L0);

    return L;
}

SourceVariable::~SourceVariable() {
	// TODO Auto-generated destructor stub
}

