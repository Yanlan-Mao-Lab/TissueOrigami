# ifndef ELEMENT_PROPERTY_SELECTION_BOX_H
# define ELEMENT_PROPERTY_SELECTION_BOX_H

# include "TissueFolding_GUI_elements.h"
# include "PropertyDisplayLayout.h"
# include <QtWidgets>

class ElementPropertySelection : public QGridLayout 
{
    Q_OBJECT
public:
    ElementPropertySelection();

    // header for pannel
    Header pannel_header = Header("Display element property:");

    // selection box for the element property to display
    DropdownMenu selection_dropdown = DropdownMenu(element_property_options);

    // the box that displays the currently selected element property
    PropertyDisplayLayout element_property_display;
    QGroupBox display_box;

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
    // max number of rows that the property display will ever use
    int max_disp_rows = 5;
};

# endif
