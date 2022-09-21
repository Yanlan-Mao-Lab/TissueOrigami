# include "TissueFolding_GUI_elements.h"
# include "ElementPropertySelection.h"

ElementPropertySelection::ElementPropertySelection() {
    // add the pannel header, align left
    addWidget(&pannel_header, 0, 0, 1, 2, AL_LEFT);

    // add the dropdown menu, align centre
    addWidget(&selection_dropdown, 0, 2, 1, 2, AL_CENTRE);

    // add the widget that will display the selected element property
    display_box.setLayout(&element_property_display);
    addWidget(&display_box, 1, 0, max_disp_rows-1, n_cols, AL_CENTRE);

    // connect: choosing a different property should change the property display
    connect(&selection_dropdown, SIGNAL(currentIndexChanged(int)), &element_property_display, SLOT(setCurrentIndex(int)));
}

void ElementPropertySelection::emitDropdownUpdate(const QString &option)
{
    emit dropdownUpdate(option);
}

void ElementPropertySelection::enableDropdownSelection(bool enabled)
{
    if (enabled)
    {
        // if we enabled the dropdown menu, we should allow user selection again
        selection_dropdown.setEnabled(true);
        // we should also update the value that is saved in the dropdown box
        emit dropdownUpdate();
    }
    else
    {
        // the dropdown menu should be disabled
        selection_dropdown.setEnabled(false);
    }
}