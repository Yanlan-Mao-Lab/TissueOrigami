# ifndef PROPERTY_DISPLAY_LAYOUT_H
# define PROPERTY_DISPLAY_LAYOUT_H

# include "GUIBuildingBlocks.h"
# include "ShapeBase.h"

# include <QtWidgets>

// there is likely a better way of doing this using maps and enums, to save on several else-if blocks later
// the possible element-properties that we can display to the user in the UI
const QStringList element_property_options = {"<Select property>",
                                              "Growth",
                                              "Growth Rate"};

/**
 * @brief The layout of the space set aside for displaying element properties, when no particular property has been chosen.
 * 
 */
class DefaultLayout : public QGridLayout 
{
public:
    DefaultLayout();

private:
    int n_rows = 1;     // Number of rows to occupy
    int n_cols = 1;     // Number of columns to occupy
    Label *placeholder_label = new Label("No property selected");    // Placeholder text to display
};

const int n_growth_rows = 3, n_growth_cols = 3;                 // number of rows/columns in the Growth matrix
const int n_growth_components = n_growth_rows * n_growth_cols;  // number of components in the Growth matrix

/**
 * @brief The layout of the space set aside for displaying element properties, when the Growth property has been selected.
 * 
 */
class GrowthLayout : public QGridLayout
{
public:
    GrowthLayout();

    /**
     * @brief Set the displayed text of the component at position (row, col)
     * 
     * @param row,col Position of component in the Growth matrix to set 
     * @param text The text to display
     */
    void setBoxValue(int row, int col, const QString &text);
    /**
     * @brief Set the displayed text of the component whose internal index is box_index
     *
     * @param box_index Internal index of the component to set, box_index = GrowthLayout::BoxIndex(row, col)
     * @param text The text to display
     */
    void setBoxValue(int box_index, const QString &text);
    /**
     * @brief Clears all values in the display boxes
     * 
     */
    void clearAllValues();
    /**
     * @brief When a new element is provided, update the displayed growth values
     *
     * @param element The new element
     */
    void newElement(std::unique_ptr<ShapeBase> *element);

    /**
     * @brief Retrieves the internal index for the component boxes, given a (row, col) index in the growth matrix.
     * 
     * box_index = row * n_growth_cols + col.
     * 
     * @param row,col Position of the component in the growth matrix
     * @return int Internal index for the box that displays this components value
     */
    int boxIndex(int row, int col);
private:
    ReadOnlyBox growth_components[n_growth_components]; // The boxes that display the components of the growth matrix
};

const int n_growthrate_components = 3;                  // Number of elements in the growth rate vector

/**
 * @brief The layout of the space set aside for displaying element properties, when the GrowthRate property has been selected
 * 
 */
class GrowthRateLayout : public QGridLayout
{
public:
    GrowthRateLayout();

    /**
     * @brief Sets the displayed text of the component with the given index
     * 
     * @param index The index of the component to set text for
     * @param text The text to set
     */
    void setComponentValue(int index, const QString &text);
    /**
     * @brief Clears all values displayed in the component boxes
     * 
     */
    void clearAllComponents();
    /**
     * @brief When a new element is provided, update the displayed growth rate values
     * 
     * @param element The new element
     */
    void newElement(std::unique_ptr<ShapeBase> *element);
private:
    ReadOnlyBox growthrate_components[n_growthrate_components];                                     // Boxes that display the values of the components
    Header *component_labels[n_growthrate_components] = { new Header("x"), new Header("y"), new Header("z") };   // Labels for the component boxes
};

/**
 * @brief Widget that can change between the various displays for each of the element properties that can be chosen for display.
 * 
 */
class PropertyDisplayLayout : public QStackedLayout 
{
    Q_OBJECT
public:
    PropertyDisplayLayout();

    void writeNewElementProperties(std::unique_ptr<ShapeBase> *element);

public slots:
    /**
     * @brief Clear all values in the various displays.
     * 
     * Recieved when the selection box is disabled, such as when an element is deselected.
     */
    void clearEntries();

private:
    QGroupBox defaultDisplay;                                       // Box containing the default layout
    DefaultLayout *defaultDisplayLayout = new DefaultLayout();      // Display for "<No property selected>"

    QGroupBox growthDisplay;                                        // Box containing the growth layout
    GrowthLayout *growthDisplayLayout = new GrowthLayout();         // Display for "Growth"

    QGroupBox growthRateDisplay;                                        // Box containing the growth rate layout
    GrowthRateLayout *growthRateDisplayLayout = new GrowthRateLayout(); // Display for "Growth Rate"
};

# endif
