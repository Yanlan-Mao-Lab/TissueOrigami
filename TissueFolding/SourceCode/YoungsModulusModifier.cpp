#include "YoungsModulusModifier.h"
#include "ShapeBase.h"
#include "Simulation.h"
#include <memory>

YoungsModulusModifier::YoungsModulusModifier(size_t TimePointNoInput,const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput, const std::array<bool,7>& WhereToApplyInput, const bool YoungsModulusMultiplierChangeRateIsLinearInput) : TimeSeriesPhysicalProperties(TimePointNoInput,nameListInput, WhenToApplyInput, WhereToApplyInput)
{
    YoungsModulusMultiplierChangeRateIsLinear = YoungsModulusMultiplierChangeRateIsLinearInput;

    CalculateYoungsModulusChangeRate();
}

void YoungsModulusModifier::CalculateYoungsModulusChangeRate()
{
    //this function takes the YoungsModulusMultiplierChangeRateIsLinear as input from the user which defines whetehr the user wants the new stiffness to be reached in a linear or exponential manner.
    nGridX = TimeSeriesGrid[0].size();
    nGridY = TimeSeriesGrid[0][0].size();
    std::vector<std::vector<std::vector<double>>> tempVec (TimePointNo-1,std::vector<std::vector<double>>(nGridX,std::vector<double>(nGridY,0)));
    YoungsModulusChangeRateGrid = tempVec;
    //YoungsModulusChangeRateGrid (TimePointNo-1,std::vector<std::vector<double>>(nX,std::vector<double>(nY,0)));
    // std::vector<std::vector<std::vector<double>>> YoungsModulusChangeRateGrid vect
    for (size_t t = 0; t<TimePointNo-1; ++t){
        double deltaT = TimeToApplyTimeSeriesGrid[t+1]-TimeToApplyTimeSeriesGrid[t];
        for (size_t i = 0; i<nGridX; ++i){
            for (size_t j = 0; j<nGridY; ++j){
                double rate=0.0;
                if (YoungsModulusMultiplierChangeRateIsLinear){
                    rate = ((TimeSeriesGrid[t+1][i][j]/TimeSeriesGrid[t][i][j])-1)/deltaT; ///< linear change of the stiffness
                }
                else{
                    rate = log(TimeSeriesGrid[t+1][i][j]/TimeSeriesGrid[t][i][j])/deltaT; ///< exponential change of the stiffness
                }
                YoungsModulusChangeRateGrid[t][i][j] = rate;
                //std::cout<<"Change rate grid is:"<<std::endl;
                //displayGrid(YoungsModulusChangeRateGrid);
            }
        }
    }
}

bool YoungsModulusModifier::calculateCurrentRateGridIndex(const double currTime){
    // this function gets the current time as input and calculates which grid of the timeseries grid to use.
    bool ModificationApplicableAtThisTimepoint = true;
    bool FoundCurrentRateGridIndex = false;
    CurrentRateGridIndex = 0;
    if (currTime>TimeToApplyTimeSeriesGrid[TimePointNo-1]){
        ModificationApplicableAtThisTimepoint = false;
    }
    else{
        for(size_t i =0; i<TimePointNo-1;++i){
            if (currTime>=TimeToApplyTimeSeriesGrid[i] && currTime<=TimeToApplyTimeSeriesGrid[i+1]){
                CurrentRateGridIndex=i;
                FoundCurrentRateGridIndex = true;
                break;
            }
        }
        if (!FoundCurrentRateGridIndex){
            std::cerr<<"could not find the application time in Young modulus modifier. "<<std::endl;
            ModificationApplicableAtThisTimepoint = false;
        }
    }
    return ModificationApplicableAtThisTimepoint;
}

bool YoungsModulusModifier::IsItApplicable(const ShapeBase* currElement){
    // this function decides whether the modifier is applicable to current element.
    bool IsModificationApplicable = false;

    if ((applyToColumnarLayer && (currElement->tissueType==0)) || (applyToPeripodialMembrane && (currElement->tissueType==1))){
        //Tissue type includes both tissue and ECM, so tissueType 0 is columnar cell layer and ECM, and tissueType 1 is peripodial layer and ECM.
        //This will apply the same map to the lateral ECM.
        if (currElement->isECMMimicing){
            if ((applyToColumnarECM && currElement->tissuePlacement==0) ||
                    (applyToColumnarECM && currElement->isECMMimimcingAtCircumference && !(currElement->tissuePlacement == 0))){              //This is the lateral ECM. But we check that the tissuePlacement is not 0 so that we don't take into account the one element that is placed where the lateral and basal ECM meet. That element is assumed to only be basal.
                IsModificationApplicable=true;
            }
        }
        else{
            if ((applyToApicalLayer && currElement->tissuePlacement==1) ||
                    (applyToMidLayer && currElement->tissuePlacement==2) ||
                    (applyToBasalLayer && (currElement->tissuePlacement==0 || currElement->atBasalBorderOfECM))){
                IsModificationApplicable=true;
            }
        } 
    }
    return IsModificationApplicable;
}


std::array<size_t,2> YoungsModulusModifier::GetElementPosition(ShapeBase* currElement){
    int IndexX = 0.0, IndexY = 0.0;
    double FracX = 1.0,  FracY = 1.0;
    currElement->getRelativePositionInTissueInGridIndex(nGridX,nGridY,IndexX, IndexY,FracX, FracY);
    std::array<size_t,2> PositionIndex{(size_t) IndexX,(size_t) IndexY};
    return PositionIndex;
}

double YoungsModulusModifier::GetCurrentMultiplierChangeRate(const size_t CurrentRateGridIndex, const std::array<size_t,2> CurrElementXYGridIndices){
    //I know th index of the rate grids within the time series that I will use
    // I need the x&y for the current element, then I can return current rate.
    size_t CurrElementXIndex = CurrElementXYGridIndices[1];
    size_t CurrElementYIndex = CurrElementXYGridIndices[2];
    return YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementXIndex][CurrElementYIndex];
}

void YoungsModulusModifier::UpdateStiffnessMultiplier(const double currTime, const double dt, ShapeBase* currElement){
    if (IsItApplicable(currElement)){
        bool ModificationApplicableAtThisTimepoint = calculateCurrentRateGridIndex(currTime);
        if (ModificationApplicableAtThisTimepoint){
            std::array<size_t,2> CurrElementXYGridIndices = GetElementPosition(currElement);
            double CurrElementMultiplierChangeRate = GetCurrentMultiplierChangeRate(CurrentRateGridIndex,CurrElementXYGridIndices);
            double CurrElementUpdatedMultiplier = 1.0;
            if (YoungsModulusMultiplierChangeRateIsLinear){
                CurrElementUpdatedMultiplier = (currElement -> StiffnessTimeSeriesMultiplier) + CurrElementMultiplierChangeRate * dt;
            }
            else{
                CurrElementUpdatedMultiplier = (currElement -> StiffnessTimeSeriesMultiplier) * exp(CurrElementMultiplierChangeRate * dt);
            }
            currElement->SetTimeseriesStiffnessMultiplier(CurrElementUpdatedMultiplier);
            currElement->updateElasticProperties();
        }
    }
}


void Simulation::ApplyAllYoungsModulusModifiers(){
    for (const auto& currentYoungsModulusModifier : AllYoungsModulusModifiers){
        for(const auto& itElement: Elements){
             currentYoungsModulusModifier->UpdateStiffnessMultiplier(currSimTimeSec, dt, itElement.get());
        }
    }
}



