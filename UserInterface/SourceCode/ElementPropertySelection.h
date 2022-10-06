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

    Header *pannel_header;  // header for pannel
    DropdownMenu *selection_dropdown;  // selection box for the element property to display
    PropertyDisplayLayout element_property_display;  // display for the currently selected element property
    QGroupBox display_box; // box wrapping the property display
    Button *save_element_properties; // box that allows the user to request the element properties be saved
    CheckBox *node_positions_on_export;  // checkbox that toggles the export of node positions with the element properties

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
    /**
     * @brief Emitted whenever element properties are exported, and the user has requested this include the associated node properties
     * 
     * @param filename File to append the node information to
     * @param element (Pointer to) element whose node-information should be written
     */
    void writeNodePositionsToFile(QString filename, std::unique_ptr<ShapeBase> *element);

private:
    int n_cols = 5;         // number of columns in the grid layout
    int max_disp_rows = 5;  // max number of rows that the property display will ever use

    QString default_name = "ElementExport.txt";                     // default name to export elements to
    QString save_dialogue_title = "Export element properties";      // dialogue box title when exporting
};

# endif
