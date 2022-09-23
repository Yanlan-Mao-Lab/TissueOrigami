# include "GUIBuildingBlocks.h"
# include "ElementPropertySelection.h"

ElementPropertySelection::ElementPropertySelection() {
    // add the pannel header, align left
    addWidget(pannel_header, 0, 0, 1, 2, AL_LEFT);

    // add the dropdown menu, align centre
    addWidget(selection_dropdown, 0, 2, 1, 2, AL_CENTRE);

    // add the widget that will display the selected element property
    display_box.setLayout(&element_property_display);
    addWidget(&display_box, 1, 0, max_disp_rows-1, n_cols, AL_CENTRE);

    // connect: choosing a different property should change the property display
    connect(selection_dropdown, SIGNAL(currentIndexChanged(int)), &element_property_display, SLOT(setCurrentIndex(int)));
}

void ElementPropertySelection::enableDropdownSelection()
{
    // the dropdown menu should be enabled
    selection_dropdown->setEnabled(true);
}
void ElementPropertySelection::disableDropdownSelection()
{
    // the dropdown menu should be disabled
    selection_dropdown->setEnabled(false);
    // the values saved should be cleared
    element_property_display.clearEntries();
}
void ElementPropertySelection::updatePropertyValues(std::unique_ptr<ShapeBase> *element)
{
    // the dropdown menu should be re-enabled
    selection_dropdown->setEnabled(true);
    // all properties in element_property_display now require a mass update
    element_property_display.writeNewElementProperties(element);
}
