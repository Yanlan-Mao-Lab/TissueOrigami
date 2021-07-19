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
    nGridX = TimeSeriesGrid[0][0].size();
    nGridY = TimeSeriesGrid[0].size();
    std::vector<std::vector<std::vector<double>>> tempVec (TimePointNo-1,std::vector<std::vector<double>>(nGridY,std::vector<double>(nGridX,0)));
    YoungsModulusChangeRateGrid = tempVec;
    //YoungsModulusChangeRateGrid (TimePointNo-1,std::vector<std::vector<double>>(nX,std::vector<double>(nY,0)));
    // std::vector<std::vector<std::vector<double>>> YoungsModulusChangeRateGrid vect
    for (size_t t = 0; t<TimePointNo-1; ++t){
        double deltaT = TimeToApplyTimeSeriesGrid[t+1]-TimeToApplyTimeSeriesGrid[t];
        for (size_t i = 0; i<nGridY; ++i){
            for (size_t j = 0; j<nGridX; ++j){
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
            if ((applyToApicalLayer && (currElement->tissuePlacement==1 || currElement->tissuePlacement==4)) ||
                    (applyToMidLayer && currElement->tissuePlacement==2) ||
                    (applyToBasalLayer && (currElement->tissuePlacement==0 || currElement->atBasalBorderOfECM))){
                IsModificationApplicable=true;
            }
        } 
    }
    return IsModificationApplicable;
}


std::array<double,4> YoungsModulusModifier::GetElementPosition(ShapeBase* currElement){
    int IndexX = 0.0, IndexY = 0.0;
    double FracX = 1.0,  FracY = 1.0;
    currElement->getRelativePositionInTissueInGridIndex(nGridX,nGridY,IndexX, IndexY,FracX, FracY);
    std::array<double,4> PositionIndex{(double) IndexX,(double) IndexY,(double) FracX, (double) FracY};
    return PositionIndex;
}

//double YoungsModulusModifier::GetCurrentMultiplierChangeRate(const size_t CurrentRateGridIndex, const std::array<size_t,2> CurrElementXYGridIndices){
//    //I know th index of the rate grids within the time series that I will use
//    // I need the x&y for the current element, then I can return current rate.
//    size_t CurrElementXIndex = CurrElementXYGridIndices[0];
//    size_t CurrElementYIndex = CurrElementXYGridIndices[1];
//    return YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex][CurrElementXIndex];
//}

double YoungsModulusModifier::GetCurrentMultiplierChangeRate(const size_t CurrentRateGridIndex, const std::array<double,4> CurrElementXYGridIndices){
    //I know th index of the rate grids within the time series that I will use
    // I need the x&y for the current element, then I can return current rate.
    double CurrElementXIndex = CurrElementXYGridIndices[0];
    double CurrElementYIndex = CurrElementXYGridIndices[1];
    double CurrElementFracX = CurrElementXYGridIndices[2];
    double CurrElementFracY = CurrElementXYGridIndices[3];
    //std::cout<<"Im in GetCurrentMultiplierChangeRate and curreElementXYGridIndices are"<<CurrentRateGridIndex<<" "<<CurrElementXYGridIndices[0]<<" "<<CurrElementXYGridIndices[1]<<" "<<CurrElementXYGridIndices[2]<<" "<<CurrElementXYGridIndices[3]<<std::endl;
    //std::cout<<"My YoungsModlulsChangeRateGrid is"<<YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex][CurrElementXIndex]<<std::endl;
    
    
    //taking growth data around 4 grid points
    double* YoungsModulusChangeRateNeighbours = new double[4];
    
    double outputYoungsModulusChangeRate;

    YoungsModulusChangeRateNeighbours[0]=YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex][CurrElementXIndex];
    YoungsModulusChangeRateNeighbours[1]=YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex][CurrElementXIndex+1];
    YoungsModulusChangeRateNeighbours[2]=YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex+1][CurrElementXIndex];
    YoungsModulusChangeRateNeighbours[3]=YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex+1][CurrElementXIndex+1];
    
    outputYoungsModulusChangeRate=(1-CurrElementFracX)*(1-CurrElementFracY)*YoungsModulusChangeRateNeighbours[0]
            +CurrElementFracX*(1-CurrElementFracY)*YoungsModulusChangeRateNeighbours[1]
            +(1-CurrElementFracX)*CurrElementFracY*YoungsModulusChangeRateNeighbours[2]
            +CurrElementFracX*CurrElementFracY*YoungsModulusChangeRateNeighbours[3];
    // return YoungsModulusChangeRateGrid[CurrentRateGridIndex][CurrElementYIndex][CurrElementXIndex];
    //std::cout<<"nGridX is "<<nGridX<<" and nGridY is "<<nGridY<<" and currElementFracX is "<<CurrElementFracX<<"and CurrElementFracY is "<<CurrElementFracY<<"CurrElementXIndex is "<<CurrElementXIndex<<"and CurrElementYIndex is "<<CurrElementYIndex<<std::endl;
    return outputYoungsModulusChangeRate;
}

void YoungsModulusModifier::updateTimeSeriesStiffnessMultiplier(const double currTime, const double dt, ShapeBase* currElement){
    if (IsItApplicable(currElement)){
        cout<<"applicable to element "<<currElement->Id<<endl;
        bool ModificationApplicableAtThisTimepoint = calculateCurrentRateGridIndex(currTime);
        if (ModificationApplicableAtThisTimepoint){
            std::array<double,4> CurrElementXYGridIndices = GetElementPosition(currElement);
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





