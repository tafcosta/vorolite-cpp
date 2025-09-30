/*
 * Source.cpp
 *
 *  Created on: 30 Sept 2025
 *      Author: Tiago Costa
 */

#include "Source.h"

Source::Source(std::vector<double> sourcePosition, std::string lightcurveFile, double lumTotal): sourcePosition(sourcePosition), lightcurveFile_(std::move(lightcurveFile)), sourceLuminosity(lumTotal) {

}

void Source::loadLightcurve_() {
    std::ifstream in(lightcurveFile_);
    if (!in) {
        throw std::runtime_error("Failed to open lightcurve file: " + lightcurveFile_);
    }

    std::vector<std::pair<double,double>> rows;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue; // skip comments and blank lines

        double t, L;
        std::istringstream iss(line);
        if (iss >> t >> L) {
            times_.push_back(t);
            luminosities_.push_back(L);
        }
    }

    if (rows.size() < 2) {
        throw std::runtime_error("Lightcurve file has fewer than 2 valid rows: " + lightcurveFile_);
    }

    // sort by time
    std::sort(rows.begin(), rows.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    times_.clear();
    luminosities_.clear();
    for (auto& r : rows) {
        if (!times_.empty() && r.first == times_.back()) {
            luminosities_.back() = r.second; // overwrite duplicate time
        } else {
            times_.push_back(r.first);
            luminosities_.push_back(r.second);
        }
    }
}

double Source::getLuminosity(double time) {

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

Source::~Source() {
	// TODO Auto-generated destructor stub
}

