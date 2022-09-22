# include "GUIBuildingBlocks.h"
# include "ElementPropertySelection.h"

ElementPropertySelection::ElementPropertySelection() {
    // add the pannel header, align left
    addWidget(&pannel_header, 0, 0, 1, 2, AL_LEFT);

    // add the dropdown menu, align centre
    addWidget(&selection_dropdown, 0, 2, 1, 2, AL_CENTRE);

    // add the widget that will display the selected element property
    display_box.setLayout(&element_property_display);
    addWidget(&display_box, 1, 1, max_disp_rows-1, n_cols-1, AL_CENTRE);

    // add the save button
    addWidget(&save_element_properties, 1, 0, max_disp_rows-1, 1, AL_CENTRE);

    // connect: choosing a different property should change the property display
    connect(&selection_dropdown, SIGNAL(currentIndexChanged(int)),
            &element_property_display, SLOT(setCurrentIndex(int)));
    // connect: clicking the save button triggers the save action
    connect(&save_element_properties, SIGNAL(clicked()),
            this, SLOT(clickedSaveElementProperties()));
}

void ElementPropertySelection::updatePropertyValues(std::unique_ptr<ShapeBase> *element)
{
    // the dropdown menu should be re-enabled
    setDropdownEnabled(true);
    // all properties in element_property_display now require a mass update
    element_property_display.writeNewElementProperties(element);
}

void ElementPropertySelection::clickedSaveElementProperties() {
    std::cout << "I clicked the button!" << endl;
}

void ElementPropertySelection::setDropdownEnabled(bool enabled) {
    // either enable or disable the dropdown menu
    selection_dropdown.setEnabled(enabled);
    // if disabling, clear the display entries
    if (!enabled) {
        element_property_display.clearEntries();
    }
    // enable/disable widgets that depend on the dropdown menu
    save_element_properties.setEnabled(enabled);
    // signal that enabling/disabling has happened
    emit dropdownIsNowEnabled(enabled);
}