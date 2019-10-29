#ifndef TIMESERIESPHYSICALPROPERTIES_H
#define TIMESERIESPHYSICALPROPERTIES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <math.h>

class TimeSeriesPhysicalProperties
{   
public:
    size_t TimePointNo;
    std::vector<std::vector<std::vector<double>>> TimeSeriesGrid;
    std::vector<double> TimeToApplyTimeSeriesGrid;

    bool applyToColumnarLayer;
    bool applyToPeripodialMembrane;
    bool applyToApicalLayer;
    bool applyToMidLayer;
    bool applyToBasalLayer;
    bool applyToPeripodialECM;
    bool applyToColumnarECM;

    TimeSeriesPhysicalProperties(size_t TimePointNoInput,const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput, const std::array<bool,7>& WhereToApplyInput);

    void ConsistencyCheck(const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput);
    bool timeStepConsistencyCheck(const double dt);
    bool ReadGrid(const std::vector<std::string>& nameListInput);
    void displayGrid(const std::vector<std::vector<std::vector<double>>>& TimeSeriesGrid);
    void DefineWhereToApply(const std::array<bool,7>& WhereToApplyInput);
};


#endif // TIMESERIESPHYSICALPROPERTIES_H
