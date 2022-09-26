# include "PropertyDisplayLayout.h"

# include <fstream>

DefaultLayout::DefaultLayout() {
    // place the placeholder label in the centre of the display
    addWidget(placeholder_label, 0, 0, n_rows, n_cols, AL_CENTRE);
}

PropertyDisplayLayout::PropertyDisplayLayout() {
    // defualt element starts as a nullptr
    current_element = nullptr;

    // setup the default display
    defaultBox = new QGroupBox; 
    defaultDisplay = new DefaultLayout;
    defaultBox->setLayout(defaultDisplay);
    addWidget(defaultBox);

    // setup the growth display
    growthBox = new QGroupBox; 
    growthDisplay = new MatrixLayout3by3;
    growthBox->setLayout(growthDisplay);
    addWidget(growthBox);
    connect(this, SIGNAL(clearEntries()), growthDisplay, SLOT(clearAllValues()));

    // setup the growth rate display
    growthRateBox = new QGroupBox; 
    growthRateDisplay = new VectorLayout3;
    growthRateBox->setLayout(growthRateDisplay);
    addWidget(growthRateBox);
    connect(this, SIGNAL(clearEntries()), growthRateDisplay, SLOT(clearAllComponents()));

    // setup the internal viscosity display
    intViscBox = new QGroupBox;
    intViscDisplay = new SingleBoxLayout;
    intViscBox->setLayout(intViscDisplay);
    addWidget(intViscBox);
    connect(this, SIGNAL(clearEntries()), intViscDisplay, SLOT(clearValue()));

    // setup the Young's modulus display
    youngModBox = new QGroupBox;
    youngModDisplay = new SingleBoxLayout;
    youngModBox->setLayout(youngModDisplay);
    addWidget(youngModBox);
    connect(this, SIGNAL(clearEntries()), youngModDisplay, SLOT(clearValue()));

    // setup the poisson ratio display
    poissonBox = new QGroupBox;
    poissonDisplay = new SingleBoxLayout;
    poissonBox->setLayout(poissonDisplay);
    addWidget(poissonBox);
    connect(this, SIGNAL(clearEntries()), poissonDisplay, SLOT(clearValue()));

    // setup the volume grwoth display
    volumeGrowthBox = new QGroupBox;
    volumeGrowthDisplay = new SingleBoxLayout;
    volumeGrowthBox->setLayout(volumeGrowthDisplay);
    addWidget(volumeGrowthBox);
    connect(this, SIGNAL(clearEntries()), volumeGrowthDisplay, SLOT(clearValue()));

    // setup the shape and size display
    shapeAndSizeBox = new QGroupBox;
    shapeAndSizeDisplay = new SingleBoxLayout;
    shapeAndSizeBox->setLayout(shapeAndSizeDisplay);
    addWidget(shapeAndSizeBox);
    connect(this, SIGNAL(clearEntries()), shapeAndSizeDisplay, SLOT(clearValue()));

    // setup the shape change rate (z) display
    shapeChangeRateZBox = new QGroupBox;
    shapeChangeRateZDisplay = new SingleBoxLayout;
    shapeChangeRateZBox->setLayout(shapeChangeRateZDisplay);
    addWidget(shapeChangeRateZBox);
    connect(this, SIGNAL(clearEntries()), shapeChangeRateZDisplay, SLOT(clearValue()));

    // setup the volume strain (Fe) display
    volStrainFeBox = new QGroupBox;
    volStrainFeDisplay = new SingleBoxLayout;
    volStrainFeBox->setLayout(volStrainFeDisplay);
    addWidget(volStrainFeBox);
    connect(this, SIGNAL(clearEntries()), volStrainFeDisplay, SLOT(clearValue()));

    // setup the strain (Dv, AP, AB) display
    strainBox = new QGroupBox;
    strainDisplay = new VectorLayout3(strainLabels);
    strainBox->setLayout(strainDisplay);
    addWidget(strainBox);
    connect(this, SIGNAL(clearEntries()), strainDisplay, SLOT(clearAllComponents()));

    // setup the strain (xy, xz, yz) display
    shearBox = new QGroupBox;
    shearDisplay = new VectorLayout3(shearLabels);
    shearBox->setLayout(shearDisplay);
    addWidget(shearBox);
    connect(this, SIGNAL(clearEntries()), shearDisplay, SLOT(clearAllComponents()));

    // upon initalisation, set to display the default display
    setCurrentWidget(defaultBox);
    // and set all the boxes to their default state by forcing a read
    // since we initilse the current element to nullptr, this just clears all the entries
    readAndUpdateElementProperties();
}
void PropertyDisplayLayout::recieveNewElement(std::unique_ptr<ShapeBase> *new_element) {
    current_element = new_element;
    // force an update
    readAndUpdateElementProperties();
}

void PropertyDisplayLayout::readAndUpdateElementProperties()
{
    // if we do not have the nullptr, we need to update all of our displays
    if (current_element) {
        // update all managed displays
        growthDisplay->fillValues((*current_element)->getFg());

        growthRateDisplay->fillValues((*current_element)->getGrowthRate());

        intViscDisplay->setDisplayValue((*current_element)->getInternalViscosity());

        youngModDisplay->setDisplayValue((*current_element)->getYoungModulus());

        poissonDisplay->setDisplayValue((*current_element)->getPoissonRatio());

        // the formula used by the source code for this value is
        // GrownVolume / ReferenceShape->Volume
        // see ShapeBase.cpp, line 1296 (type=5 for this field)
        // NOTE: this is a ratio, not absolute by the looks of things!
        double GrownVolume = (*current_element)->GrownVolume / (*current_element)->ReferenceShape->Volume;
        volumeGrowthDisplay->setDisplayValue(GrownVolume);

        shapeAndSizeDisplay->setDisplayValue((*current_element)->calculateEmergentShapeOrientation());

        std::array<double, 6> shapeChangeRate = (*current_element)->getShapeChangeRate();
        // change rate in z is stored at index 2
        shapeChangeRateZDisplay->setDisplayValue((shapeChangeRate[2]));

        /* Strains are computed by passing pointers;
        getStrain(type, pointer_to_float) assigns the value of the type of strain to the address.
        type takes the values:
        0 - Fe strain
        1/2/3 - DV/AP/AB strain
        4/5/6(?) - xy/yz/xz strain
        */
        float strainFe, strainNonPlanarComps[3], strainPlanarComps[3];
        (*current_element)->getStrain(0, strainFe);
        for(int i=0; i<3; i++) {
            (*current_element)->getStrain(i+1, strainNonPlanarComps[i]);
            (*current_element)->getStrain(i+4, strainPlanarComps[i]);
        }
        volStrainFeDisplay->setDisplayValue(strainFe);
        strainDisplay->fillValues(strainNonPlanarComps);
        shearDisplay->fillValues(strainPlanarComps);
    }
    // otherwise, deselection has occurred and we should clear everything
    else {
        emit clearEntries();
    }
}

/*
void GrowthLayout::writeTo(string filename, bool write_header) {
    // append to the file rather than overwriting, as we will be writing in sequence
    ofstream file(filename, std::ios_base::app);
    // write header if requested
    if (write_header) {
        file << "== Growth ==\n";
    }
    for(int row=0; row<n_growth_rows; row++) {
        for(int col=0; col<n_growth_cols; col++) {
            int box_index = boxIndex(row, col);
            file << growth_component_values[box_index] << ",";
        }
        file << "\n";
    }
    // close after finishing
    file.close();
}
void GrowthRateLayout::writeTo(string filename, bool write_header) {
    // append to the file rather than overwriting, as we will be writing in sequence
    ofstream file(filename, std::ios_base::app);
    // write header if requested
    if (write_header) {
        file << "== Growth Rate ==\n";
    }
    for(int index=0; index<n_growthrate_components; index++) {
        file << growthrate_component_values[index] << ",";
    }
    file << "\n";
    // close after finishing
    file.close();
}

void PropertyDisplayLayout::writeToFile(const QString &filename) {
    string f_name = filename.toStdString();
    // write, in succession, the outputs
    // Growth
    growthDisplayLayout->writeTo(f_name);
    // GrowthRate
    growthRateDisplayLayout->writeTo(f_name);
}

*/