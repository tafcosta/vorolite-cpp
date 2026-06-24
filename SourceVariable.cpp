/*
 * SourceVariable.cpp
 *
 *  Created on: 13 Nov 2025
 *      Author: Tiago Costa
 */

#include "SourceVariable.h"

SourceVariable::SourceVariable(std::vector<double> sourcePosition,
                               double lumTotal, std::string lightcurveFile)
    : Source(sourcePosition, lumTotal),
      lightcurveFile_(lightcurveFile)
{
    loadLightcurve_();
}


void SourceVariable::loadLightcurve_()
{
    std::ifstream in(lightcurveFile_);

    if (!in) {
        throw std::runtime_error("Cannot open light-curve file: " + lightcurveFile_);
    }

    std::cout << "Loading light-curve from file: " << lightcurveFile_ << std::endl;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Accept headers like "Time, Lion" and mixed comma/space separators.
        std::replace(line.begin(), line.end(), ',', ' ');
        std::stringstream ss(line);
        double t, L;
        if (!(ss >> t >> L)) {
            continue;
        }
 
        times_.push_back(t);
        luminosities_.push_back(L);
    }

    if (times_.empty() || luminosities_.empty()) {
        throw std::runtime_error(
            "Light-curve file contains no valid numeric (time, luminosity) pairs: " + lightcurveFile_);
    }

    if (times_.size() != luminosities_.size()) {
        throw std::runtime_error(
            "Light-curve parsing produced mismatched time/luminosity lengths: " + lightcurveFile_);
    }

    std::cout << "First time " << times_.front() << ", first luminosity " << luminosities_.front() << std::endl;
    std::cout << "Last time " << times_.back() << ", last luminosity " << luminosities_.back() << std::endl;
}

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
