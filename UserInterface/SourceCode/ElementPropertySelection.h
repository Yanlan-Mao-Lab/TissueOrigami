# ifndef ELEMENT_PROPERTY_SELECTION_BOX_H
# define ELEMENT_PROPERTY_SELECTION_BOX_H

# include "ShapeBase.h"

# include "GUIBuildingBlocks.h"
# include "PropertyDisplayLayout.h"
# include <QtWidgets>

class ElementPropertySelection : public QGridLayout 
{
    Q_OBJECT
public:
    ElementPropertySelection();

    // header for pannel
    Header *pannel_header = new Header("Display element property:");

    // selection box for the element property to display
    DropdownMenu *selection_dropdown = new DropdownMenu(element_property_options);

    // the box that displays the currently selected element property
    PropertyDisplayLayout element_property_display;
    QGroupBox display_box;

    // the box that allows the user to request the element properties be saved
    Button *save_element_properties = new Button("Export\n properties of\n selected\n element");

public slots:
    /**
     * @brief Updates the values stored in the element_property_display to match those of the new element
     *
     * @param element The new element whose properties should be displayed. nullptr is interpretted as deselection.
     */
    void updatePropertyValues(std::unique_ptr<ShapeBase> *element);
    /**
     * @brief Enables or disables the dropdown menu.
     * 
     * Emits dropdownIsNowEnabled(bool) upon execution.
     * 
     * @param enabled Whether to enable (true) or disable (false) the dropdown menu
     */
    void setDropdownEnabled(bool enabled);

public slots:
    /**
     * @brief Exports the current element properties to a file.
     * 
     */
    void clickedSaveElementProperties();

signals:
    /**
     * @brief Emitted whenever the dropdown menu is enabled or disabled.
     * 
     * Signal should be connected to the setEnabled(bool) slot of all attributes that should be enabled IF AND ONLY IF the dropdown menu object is enabled.
     * 
     * @param enabled Whether the dropdown has been enabled (true) or disabled (false)
     */
    void dropdownIsNowEnabled(bool enabled);

private:
    int n_cols = 5;         // number of columns in the grid layout
    int max_disp_rows = 5;  // max number of rows that the property display will ever use

    QString default_name = "ElementExport.txt";                     // default name to export elements to
    QString save_dialogue_title = "Export element properties";      // dialogue box title when exporting
};

# endif
