/*
 * Source.cpp
 *
 *  Created on: 30 Sept 2025
 *      Author: ntc132
 */

#include "Source.h"

Source::Source(std::vector<double> sourcePosition, std::string lightcurveFile, double lumTotal): sourcePosition(sourcePosition), lightcurve(lightcurveFile), sourceLuminosity(lumTotal) {

}

double Source::getLuminosity(double time){
	return sourceLuminosity;
}

Source::~Source() {
	// TODO Auto-generated destructor stub
}

