# include "TissueFolding_GUI_elements.h"
# include "ElementPropertySelection.h"

ElementPropertySelection::ElementPropertySelection() {
    // add the pannel header, span all columns and align left
    addWidget(&pannel_header, 0, 0, 1, 2, AL_LEFT);

    // add the dropdown menu, span all columns and align centre
    addWidget(&select_element_property_dropdown, 0, 2, 1, 2, AL_CENTRE);

    // add the property display - this will dynamically change so we need to be careful!
    addWidget(&select_element_property_display, 1, 0, 1, n_cols, AL_CENTRE);

    // connect the option selection in the dropdown menu to the dropdownUpdate signal
    connect(&select_element_property_dropdown, SIGNAL(currentTextChanged(const QString &)), this, SLOT(emitDropdownUpdate(const QString &)));
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
        select_element_property_dropdown.setEnabled(true);
        // we should also update the value that is saved in the dropdown box
        emit dropdownUpdate();
    }
    else
    {
        // the dropdown menu should be disabled
        select_element_property_dropdown.setEnabled(false);
        // remove any text in the property value box
        select_element_property_display.setText("-");
    }
}