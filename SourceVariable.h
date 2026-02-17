/*
 * SourceVariable.h
 *
 *  Created on: 13 Nov 2025
 *      Author: Tiago Costa
 */

#ifndef SOURCEVARIABLE_H_
#define SOURCEVARIABLE_H_


#include "Source.h"

class SourceVariable : public Source {
public:
	SourceVariable(std::vector<double> sourcePosition, double lumTotal, std::string lightcurveFile);
    ~SourceVariable() override;

    double getLuminosity(double time) override;

private:
    std::string lightcurveFile_;
    std::vector<double> times_;
    std::vector<double> luminosities_;

    void loadLightcurve_();
};


#endif /* SOURCEVARIABLE_H_ */
