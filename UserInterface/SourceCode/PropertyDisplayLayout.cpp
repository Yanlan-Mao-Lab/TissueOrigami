# include "PropertyDisplayLayout.h"

DefaultLayout::DefaultLayout() {
    // place the placeholder label in the centre of the display
    addWidget(&placeholder_label, 0, 0, n_rows, n_cols, AL_CENTRE);
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
    }
}

int GrowthLayout::boxIndex(int row, int col) {
    if ((row<0) || (row>=n_growth_rows) || (col<0) || (col>=n_growth_cols)) {
        throw std::runtime_error("Growth matrix index out of bounds\n");
    }
    return row*n_growth_cols + col;
}

GrowthRateLayout::GrowthRateLayout() {
    for(int index=0; index<n_growthrate_components; index++) {
        // add the box labels to the display
        addWidget(&component_labels[index], index, 0, 1, 1, AL_RIGHT);
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
    }
}

PropertyDisplayLayout::PropertyDisplayLayout() {
    // setup the default display
    defaultDisplay.setLayout(&defaultDisplayLayout);
    addWidget(&defaultDisplay);

    // setup the growth display
    growthDisplay.setLayout(&growthDisplayLayout);
    addWidget(&growthDisplay);

    // setup the growth rate display
    growthRateDisplay.setLayout(&growthRateDisplayLayout);
    addWidget(&growthRateDisplay);

    // upon initalisation, set to display the default display
    setCurrentWidget(&defaultDisplay);
}

void PropertyDisplayLayout::clearEntries() {
    // remove any values that may have been written to the displays
    growthDisplayLayout.clearAllValues();
    growthRateDisplayLayout.clearAllComponents();
}