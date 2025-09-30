/*
 * Source.h
 *
 *  Created on: 30 Sept 2025
 *      Author: Tiago Costa
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "common_includes.h"

class Source {
public:
	Source(std::vector<double> sourcePosition, std::string lightcurveFile, double lumTotal);
	virtual ~Source();

	double getLuminosity(double time);

	std::vector<double> sourcePosition;
	double sourceLuminosity;

private:

    std::string lightcurveFile_;

    std::vector<double> times_;
    std::vector<double> luminosities_;

    void loadLightcurve_();


};

#endif /* SOURCE_H_ */
