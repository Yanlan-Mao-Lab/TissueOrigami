# include "PropertyDisplayLayout.h"

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
    connect(this, SIGNAL(updateEntries()), growthDisplay, SLOT(fillValues((*current_element)->getFg())));

    // setup the growth rate display
    growthRateBox = new QGroupBox; 
    growthRateDisplay = new VectorLayout3;
    growthRateBox->setLayout(growthRateDisplay);
    addWidget(growthRateBox);
    connect(this, SIGNAL(clearEntries()), growthRateDisplay, SLOT(clearAllComponents()));
    connect(this, SIGNAL(updateEntries()), growthRateDisplay, SLOT(fillValues((*current_element)->getGrowthRate())));

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
    strainDV_AP_ABBox = new QGroupBox;
    strainDV_AP_ABDisplay = new VectorLayout3(DV_AP_AB_Labels);
    strainDV_AP_ABBox->setLayout(strainDV_AP_ABDisplay);
    addWidget(strainDV_AP_ABBox);
    connect(this, SIGNAL(clearEntries()), strainDV_AP_ABDisplay, SLOT(clearAllComponents()));

    // setup the strain (xy, xz, yz) display
    strainPlanarDirsBox = new QGroupBox;
    strainPlanarDirsDisplay = new VectorLayout3(planarDirsLabels);
    strainPlanarDirsBox->setLayout(strainPlanarDirsDisplay);
    addWidget(strainPlanarDirsBox);
    connect(this, SIGNAL(clearEntries()), strainPlanarDirsDisplay, SLOT(clearAllComponents()));

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
        // update growth display
        //growthDisplay->fillValues((*current_element)->getFg());
        //growthRateDisplay->fillValues((*current_element)->getGrowthRate());
        emit updateEntries();
    }
    // otherwise, deselection has occurred and we should clear everything
    else {
        emit clearEntries();
    }
}

