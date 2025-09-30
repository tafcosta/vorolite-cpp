/*
 * Source.h
 *
 *  Created on: 30 Sept 2025
 *      Author: ntc132
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
	std::string& lightcurve;
	double sourceLuminosity;

};

#endif /* SOURCE_H_ */
