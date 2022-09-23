# ifndef ELEMENT_PROPERTY_SELECTION_BOX_H
# define ELEMENT_PROPERTY_SELECTION_BOX_H

# include "ShapeBase.h"

# include "GUIBuildingBlocks.h"
# include "PropertyDisplayLayout.h"
# include <QtWidgets>

class ElementPropertySelection : public QGridLayout 
{
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
     * @brief Enables options to be selected from the property selection box
     * 
     */
    void enableDropdownSelection();
    /**
     * @brief Disables the property selection box
     * 
     */
    void disableDropdownSelection();
    /**
     * @brief Updates the values stored in the element_property_display to match those of the new element
     *
     * @param element The new element whose properties should be displayed
     */
    void updatePropertyValues(std::unique_ptr<ShapeBase> *element);

private:
    // number of columns in the grid layout
    int n_cols = 5;
    // max number of rows that the property display will ever use
    int max_disp_rows = 5;
};

# endif
