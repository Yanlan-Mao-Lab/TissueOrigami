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

QString PropertyDisplayLayout::makeHeader(const QString &section_name) {
    QString section_header_rev = section_header;
    std::reverse(section_header_rev.begin(), section_header_rev.end());
    return section_header + section_name + section_header_rev + "\n";
}
void PropertyDisplayLayout::writeToFile(const QString &filename, bool write_headers) {
    // write the values of the managed properties to an output file.
    // we have to do this sequentially as we can't allow processes to write to the file simultaneously, as the data will be muddled

    // open the file, in append mode if necessary (might have asked for node positions to be written)
    ofstream file(filename.toStdString(), std::ios_base::app);

    // if writing with headers...
    if (write_headers) {
        // write growth
        file << makeHeader("Growth").toStdString();
        gsl_matrix *growth = (*current_element)->getFg();
        for (int row = 0; row < growthDisplay->n_rows; row++)
        {
            for (int col = 0; col < growthDisplay->n_cols; col++)
            {
                int box_index = growthDisplay->boxIndex(row, col);
                file << gsl_matrix_get(growth, row, col) << ",";
            }
            file << "\n";
        }

        // write growth rate
        file << makeHeader("Volume Growth Rate (xyz)").toStdString();
        std::array<double, 3> growth_rate = (*current_element)->getGrowthRate();
        for (int index = 0; index < growthRateDisplay->n_comps; index++)
        {
            file << growth_rate[index] << ",";
        }
        file << "\n";

        // write scalar properties
        file << makeHeader("Internal Viscosity").toStdString();
        file << (*current_element)->getInternalViscosity() << "\n";
        file << makeHeader("Young's modulus").toStdString();
        file << (*current_element)->getYoungModulus() << "\n";
        file << makeHeader("Poisson Ratio").toStdString();
        file << (*current_element)->getPoissonRatio() << "\n";

        // the formula used by the source code for this value is
        // GrownVolume / ReferenceShape->Volume
        // see ShapeBase.cpp, line 1296 (type=5 for this field)
        // NOTE: this is a ratio, not absolute by the looks of things!
        double GrownVolume = (*current_element)->GrownVolume / (*current_element)->ReferenceShape->Volume;
        file << makeHeader("Volume Growth Ratio").toStdString();
        file << GrownVolume << "\n";

        file << makeHeader("Emergent Shape and Size").toStdString();

        file << (*current_element)->calculateEmergentShapeOrientation() << "\n";

        file << makeHeader("Shape Change Rate (z)").toStdString();
        std::array<double, 6> shapeChangeRate = (*current_element)->getShapeChangeRate();
        // change rate in z is stored at index 2
        file << shapeChangeRate[2] << "\n";

        // write the strain and the shear
        float strainFe, strainComps[3], shearComps[3];
        (*current_element)->getStrain(0, strainFe);
        for (int i = 0; i < 3; i++)
        {
            (*current_element)->getStrain(i + 1, strainComps[i]);
            (*current_element)->getStrain(i + 4, shearComps[i]);
        }
        file << makeHeader("Volumetric strain (via Fe)").toStdString();
        file << strainFe << "\n";
        file << makeHeader("Strain (DV, AP, AB)").toStdString();
        for (int index = 0; index < strainDisplay->n_comps; index++)
        {
            file << strainComps[index] << ",";
        }
        file << "\n";
        file << makeHeader("Shear (xy, yz, xz)").toStdString();
        for (int index = 0; index < shearDisplay->n_comps; index++)
        {
            file << shearComps[index] << ",";
        }
        file << "\n";
    }
    else {
        // write growth
        gsl_matrix *growth = (*current_element)->getFg();
        for (int row = 0; row < growthDisplay->n_rows; row++)
        {
            for (int col = 0; col < growthDisplay->n_cols; col++)
            {
                int box_index = growthDisplay->boxIndex(row, col);
                file << gsl_matrix_get(growth, row, col) << ",";
            }
            file << "\n";
        }

        // write growth rate
        std::array<double, 3> growth_rate = (*current_element)->getGrowthRate();
        for (int index = 0; index < growthRateDisplay->n_comps; index++)
        {
            file << growth_rate[index] << ",";
        }
        file << "\n";

        // write scalar properties
        file << (*current_element)->getInternalViscosity() << "\n";
        file << (*current_element)->getYoungModulus() << "\n";
        file << (*current_element)->getPoissonRatio() << "\n";

        // the formula used by the source code for this value is
        // GrownVolume / ReferenceShape->Volume
        // see ShapeBase.cpp, line 1296 (type=5 for this field)
        // NOTE: this is a ratio, not absolute by the looks of things!
        double GrownVolume = (*current_element)->GrownVolume / (*current_element)->ReferenceShape->Volume;
        file << GrownVolume << "\n";

        file << (*current_element)->calculateEmergentShapeOrientation() << "\n";

        std::array<double, 6> shapeChangeRate = (*current_element)->getShapeChangeRate();
        // change rate in z is stored at index 2
        file << shapeChangeRate[2] << "\n";

        // write the strain and the shear
        float strainFe, strainComps[3], shearComps[3];
        (*current_element)->getStrain(0, strainFe);
        for (int i = 0; i < 3; i++)
        {
            (*current_element)->getStrain(i + 1, strainComps[i]);
            (*current_element)->getStrain(i + 4, shearComps[i]);
        }
        file << strainFe << "\n";
        for (int index = 0; index < strainDisplay->n_comps; index++)
        {
            file << strainComps[index] << ",";
        }
        file << "\n";
        for (int index = 0; index < shearDisplay->n_comps; index++)
        {
            file << shearComps[index] << ",";
        }
        file << "\n";
    }
}
