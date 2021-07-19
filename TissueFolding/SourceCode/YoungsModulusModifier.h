#ifndef YOUNGSMODULUSMODIFIER_H
#define YOUNGSMODULUSMODIFIER_H

#include "TimeSeriesPhysicalProperties.h"
#include <iostream>
#include "ShapeBase.h"
#include <array>

class YoungsModulusModifier : public TimeSeriesPhysicalProperties
{
public:
    double CurrentTimeStep;
    double NextTimeStep;

    bool YoungsModulusMultiplierChangeRateIsLinear;

    size_t nGridX;
    size_t nGridY;
    size_t CurrentRateGridIndex;

    std::array<double,4> CurrElementXYGridIndices;

    std::vector<double> CurrentStepYoungsModulusGrid;
    std::vector<double> NextStepYoungsModulusGrid;
    std::vector<std::vector<std::vector<double>>> YoungsModulusChangeRateGrid;
    std::vector<std::vector<double>> currentMultiplierChangeRate;

    YoungsModulusModifier(size_t TimePointNoInput,const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput, const std::array<bool,7>& WhereToApplyInput, bool YoungsModulusMultiplierChangeRateIsLinear);

    void CalculateYoungsModulusChangeRate();                                                                                                      ///< this function calculates the Young's modulus change rate from the input time-series grid.
    bool calculateCurrentRateGridIndex(const double currTime);                                                                                    ///<this function gets the current time as input and calculates which grid of the timeseries grid to use.
    bool IsItApplicable(const ShapeBase* currElement);                                                                                            ///< this function gets the element and decides whether the change in the Young's modulus is applicable to it.
    std::array<double,4> GetElementPosition(ShapeBase* currElement);                                                                              ///< this function takes the size of the grid as input and calculates the element position on the grid.
    double GetCurrentMultiplierChangeRate(const size_t CurrentRateGridIndex, const std::array<double,4> CurrElementXYGridIndices);                ///< this function gets the current time index and the position of the element, looks through the timeseries multiplier rate grid and finds the multiplier change rate that is applicable to the element.
    void updateTimeSeriesStiffnessMultiplier(const double currTime, const double dt, ShapeBase* currElement);                                               ///< this function calculates the updated value of the element multiplier.
};


#endif // YOUNGSMODULUSMODIFIER_H
