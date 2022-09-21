# ifndef ELEMENT_PROPERTY_SELECTION_BOX_H
# define ELEMENT_PROPERTY_SELECTION_BOX_H

# include "TissueFolding_GUI_elements.h"
# include <QtWidgets>

// there is likely a better way of doing this using maps and enums, to save on several else-if blocks later
// the possible element-properties that we can display to the user in the UI
const QStringList element_property_options = {"Growth",
                                              "Growth Rate"};

class ElementPropertySelection : public QGridLayout 
{
    Q_OBJECT
public:
    ElementPropertySelection();

    // header for pannel
    Header pannel_header = Header("Display element property:");

    // selection box for the element property to display
    DropdownMenu select_element_property_dropdown = DropdownMenu(element_property_options);

    // Box to display the value of the select_element_property (needs to change with the property though!!!!)
    ReadOnlyBox select_element_property_display = ReadOnlyBox("-");

    /**
     * @brief Enables the dropdown menu and signals that the property box should also be updated.
     *
     * Emits dropdownUpdate() if enabled is true.
     *
     * @param enabled Whether to enable (true) or disable (false) the dropdown options
     */
    void enableDropdownSelection(bool enabled);

public slots:
    // emit the dropdownUpdate signal when the user changes the option in the element_property_dropdown menu
    void emitDropdownUpdate(const QString &option);

signals:
    // emitted when the element property selection dropdown box has been enabled, so needs updating
    void dropdownUpdate();
    // emitted when the element property selection dropdown box has been changed, so needs updating
    void dropdownUpdate(const QString &option);

private:
    // number of columns in the grid layout
    int n_cols = 5;
};

# endif
