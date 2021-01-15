#include "TimeSeriesPhysicalProperties.h"

TimeSeriesPhysicalProperties::TimeSeriesPhysicalProperties(size_t TimePointNoInput, const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput, const std::array<bool,7>& WhereToApplyInput)
{
    TimePointNo=TimePointNoInput;
    TimeToApplyTimeSeriesGrid=WhenToApplyInput;

    ConsistencyCheck(nameListInput, WhenToApplyInput);
    ReadGrid(nameListInput);
    //displayGrid(TimeSeriesGrid);
    DefineWhereToApply(WhereToApplyInput);
}


void TimeSeriesPhysicalProperties::ConsistencyCheck(const std::vector<std::string>& nameListInput, const std::vector<double>& WhenToApplyInput){

    size_t nameListSize=nameListInput.size();
    size_t WhenToApplySize=WhenToApplyInput.size();
    if (nameListSize != TimePointNo || WhenToApplySize!=TimePointNo){
        std::cerr<<"The number of timepoints do not match the size of the inputs"<<std::endl;
    }
}

bool TimeSeriesPhysicalProperties::timeStepConsistencyCheck(const double dt){
    size_t WhenToApplySize=TimeToApplyTimeSeriesGrid.size();
    for (size_t i=0;i<WhenToApplySize;++i){
        if (fmod(TimeToApplyTimeSeriesGrid[i],dt) > 0){
            std::cerr<<"The timeseries "<<i+1<<" is not dividable to your timestep "<< dt<<std::endl;
            return false;
        }
    }
    return true;
}
bool TimeSeriesPhysicalProperties::ReadGrid(const std::vector<std::string>& nameListInput){

    //open physical property file
    for (size_t nT =0; nT<TimePointNo; ++nT){
        size_t nX;         //number of rows of the physical property grid
        size_t nY;         //number of columns of the physical property grid

        double PhysicalProperty;

        std:: string filename = nameListInput[nT];
        std::cout<<"physical property filename is "<<filename<<std::endl;
        std::ifstream PhysicalPropertyFile;
        PhysicalPropertyFile.open(filename,std::ios::in);

        //give an error if file doesn't exist or cannot be opened
        if (!(PhysicalPropertyFile.good() && PhysicalPropertyFile.is_open())){
            std::cerr<<"could not open physical property file : "<<filename<<std::endl;
            return false;
        }
        //std::cout<<"filename is : "<<filename<<std::endl;
        //read first row of the file which defines number of columns (nX) and number of rows (nY)
        PhysicalPropertyFile >> nX;
        PhysicalPropertyFile >> nY;

        if (nT==0){
            std::vector<std::vector<std::vector<double>>> tempVecTimeSeries (TimePointNo,std::vector<std::vector<double>>(nY,std::vector<double>(nX,0)));
            TimeSeriesGrid=tempVecTimeSeries;
        }
        std::cout<<"The TimeSeriesPhysicalPropertiesGrid before going into loop is:"<<std::endl;
        displayGrid(TimeSeriesGrid);

        for (size_t i= 0; i<nY; ++i){
            for(size_t j=0; j<nX; ++j){
                PhysicalPropertyFile>>PhysicalProperty;

                TimeSeriesGrid[nT][i][j]=PhysicalProperty;
            }
        }

//        //This section creates the timeseries grid.
//        std::vector<std::vector<double>> oneGrid;

//        for (size_t i= 0; i<nY; ++i){
//            std::vector<double> oneRow(nX,0);
//            for(size_t j=0; j<nX; ++j){
//                PhysicalPropertyFile>>PhysicalProperty;
//                oneRow[j]=PhysicalProperty;
//            }
//            oneGrid.push_back(oneRow);
//        }
//        TimeSeriesGrid.push_back(oneGrid);

        PhysicalPropertyFile.close();
    }
    std::cout<<"The TimeSeriesPhysicalPropertiesGrid is:"<<std::endl;
    displayGrid(TimeSeriesGrid);
    return true;
}


//This section displays the Young's modulus timeseries grid in the output tmp file.
//& is passing with address, so that you are not giving the whole grid to the function, but only its adress.
//Once you have its address, you can change it, which is dangerous if this is not your intention.
//a display function has no reason to change it.
//const makes sure it is not changeable (constant)
void TimeSeriesPhysicalProperties::displayGrid(const std::vector<std::vector<std::vector<double>>>& TimeSeriesGrid){
    //for (size_t k= 0; k<nT; ++k){
    for (auto gridOfCurrentTimePoint : TimeSeriesGrid){
        //std::cout<<"Inside displayGrid, the timepoint is : "<<k<<" and the modulus grid is:"<<std::endl;
        for (auto oneRowOfCurrentGrid : gridOfCurrentTimePoint){
            for (auto oneXYValue : oneRowOfCurrentGrid){
            //for (size_t i= 0; i<nX; ++i){
            //for (size_t j= 0; j<nY; ++j){
              std::cout<<oneXYValue<<"  ";
                //std::cout<<TimeSeriesGrid[k][i][j]<<"  ";
            }
            std::cout<<std::endl; //end the line between rows.
        }
        std::cout<<std::endl; //end the line between rows.
    }
}


void TimeSeriesPhysicalProperties::DefineWhereToApply(const std::array<bool,7>& WhereToApplyInput){
    applyToColumnarLayer = WhereToApplyInput[0];
    applyToPeripodialMembrane = WhereToApplyInput[1];
    applyToApicalLayer = WhereToApplyInput[2];
    applyToMidLayer = WhereToApplyInput[3];
    applyToBasalLayer = WhereToApplyInput[4];
    applyToPeripodialECM = WhereToApplyInput[5];
    applyToColumnarECM = WhereToApplyInput[6];
}

