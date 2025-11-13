/*
 * Source.cpp
 *
 *  Created on: 30 Sept 2025
 *      Author: Tiago Costa
 */

#include "Source.h"

Source::Source(std::vector<double> sourcePosition, double lumTotal): sourcePosition(sourcePosition), sourceLuminosity(lumTotal) {

}

double Source::getLuminosity(double time) {

return sourceLuminosity;
}

Source::~Source() {
	// TODO Auto-generated destructor stub
}

