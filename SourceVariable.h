/*
 * SourceVariable.h
 *
 *  Created on: 13 Nov 2025
 *      Author: ntc132
 */

#ifndef SOURCEVARIABLE_H_
#define SOURCEVARIABLE_H_

#include "Source.h"

class SourceVariable: public Source {
public:
	SourceVariable();
	virtual ~SourceVariable();

	double getLuminosity(double time);

private:

    std::string lightcurveFile_;

    std::vector<double> times_;
    std::vector<double> luminosities_;

    void loadLightcurve_();
};

#endif /* SOURCEVARIABLE_H_ */
