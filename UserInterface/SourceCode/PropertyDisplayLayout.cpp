# include "PropertyDisplayLayout.h"

# include <fstream>

DefaultLayout::DefaultLayout() {
    // place the placeholder label in the centre of the display
    addWidget(placeholder_label, 0, 0, n_rows, n_cols, AL_CENTRE);
}

GrowthLayout::GrowthLayout() {
    // add each box corresponding to a component to the display, in a grid
    for(int row=0; row<n_growth_rows; row++) {
        for(int col=0; col<n_growth_cols; col++) {
            int box_index = boxIndex(row, col);
            growth_components[box_index].setAlignment(AL_CENTRE);
            addWidget(&growth_components[box_index], row, col, 1, 1, AL_CENTRE);
        }
    }
}

void GrowthLayout::setBoxValue(int row, int col, const QString &text) {
    setBoxValue(boxIndex(row, col), text);
}
void GrowthLayout::setBoxValue(int box_index, const QString &text) {
    growth_components[box_index].setText(text);
}
void GrowthLayout::clearAllValues() {
    for(int box_index=0; box_index<n_growth_components; box_index++) {
        growth_components[box_index].setText("-");
        growth_component_values[box_index] = 0.;
    }
}
void GrowthLayout::newElement(std::unique_ptr<ShapeBase> *element) {
    // extract growth matrix
    gsl_matrix *growth_matrix = (*element)->getFg();
    // store the new values, and write them to the display
    for(int row=0; row<n_growth_rows; row++) {
        for(int col=0; col<n_growth_cols; col++) {
            int box_index = boxIndex(row, col);
            growth_component_values[box_index] = gsl_matrix_get(growth_matrix, row, col);
            setBoxValue(row, col, QString::number(growth_component_values[box_index]));
        }
    }
}
int GrowthLayout::boxIndex(int row, int col) {
    if ((row<0) || (row>=n_growth_rows) || (col<0) || (col>=n_growth_cols)) {
        throw std::runtime_error("Growth matrix index out of bounds\n");
    }
    return row*n_growth_cols + col;
}
double GrowthLayout::getGrowthComponent(int row, int col) {
    return getGrowthComponent(boxIndex(row, col));
}
double GrowthLayout::getGrowthComponent(int box_index) {
    return growth_component_values[box_index];
}
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

GrowthRateLayout::GrowthRateLayout() {
    for(int index=0; index<n_growthrate_components; index++) {
        // add the box labels to the display
        addWidget(component_labels[index], index, 0, 1, 1, AL_RIGHT);
        // add each box corresponding to a component to the display, as a column vector
        growthrate_components[index].setAlignment(AL_CENTRE);
        addWidget(&growthrate_components[index], index, 1, 1, 2, AL_CENTRE);
    }
}
void GrowthRateLayout::setComponentValue(int index, const QString &text) {
    growthrate_components[index].setText(text);
}
void GrowthRateLayout::clearAllComponents() {
    for(int index=0; index<n_growthrate_components; index++) {
        growthrate_components[index].setText("-");
        growthrate_component_values[index] = 0.;
    }
}
void GrowthRateLayout::newElement(std::unique_ptr<ShapeBase> *element) {
    // extract growth rate
    std::array<double, 3> growth_rate = (*element)->getGrowthRate();
    // write growth rates
    for(int index=0; index<n_growthrate_components; index++) {
        growthrate_component_values[index] = growth_rate[index];
        setComponentValue(index, QString::number(growthrate_component_values[index]));
    }
}
double GrowthRateLayout::getGrowthRateComponent(int index) {
    return growthrate_component_values[index];
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

PropertyDisplayLayout::PropertyDisplayLayout() {
    // setup the default display
    defaultDisplay.setLayout(defaultDisplayLayout);
    addWidget(&defaultDisplay);

    // setup the growth display
    growthDisplay.setLayout(growthDisplayLayout);
    addWidget(&growthDisplay);

    // setup the growth rate display
    growthRateDisplay.setLayout(growthRateDisplayLayout);
    addWidget(&growthRateDisplay);

    // upon initalisation, set to display the default display
    setCurrentWidget(&defaultDisplay);
}
void PropertyDisplayLayout::writeNewElementProperties(std::unique_ptr<ShapeBase> *element) {
    growthDisplayLayout->newElement(element);
    growthRateDisplayLayout->newElement(element);
}
void PropertyDisplayLayout::writeToFile(const QString &filename) {
    string f_name = filename.toStdString();
    // write, in succession, the outputs
    // Growth
    growthDisplayLayout->writeTo(f_name);
    // GrowthRate
    growthRateDisplayLayout->writeTo(f_name);
}
void PropertyDisplayLayout::clearEntries() {
    // remove any values that may have been written to the displays
    growthDisplayLayout->clearAllValues();
    growthRateDisplayLayout->clearAllComponents();
}