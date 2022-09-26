# ifndef PROPERTY_DISPLAY_LAYOUT_H
# define PROPERTY_DISPLAY_LAYOUT_H

# include "GUIBuildingBlocks.h"
# include "ShapeBase.h"

# include <QtWidgets>

// there is likely a better way of doing this using maps and enums, to save on several else-if blocks later
// the possible element-properties that we can display to the user in the UI
// NOTE: EXTERNAL VISCOSITY IS NODE-BASED!!!
const QStringList element_property_options = {"<Select property>",          // default, no selection
                                              "Growth",                     // 3-by-3 matrix
                                              "Volume Growth Rate (xyz)",   // 3-vector (x,y,z)
                                              "Internal Viscosity",         // getInternalViscosity()
                                              "Young Modulus",              // getYoungModulus()
                                              "Poisson Ratio",              // getPoissonRatio()
                                              "Volume Growth",              // GrownVolume/ReferenceShape->Volume
                                              "Emergent Shape & Size",      // calculateEmergentShapeOrientation
                                              "ShapeChangeRate (z)",        // getShapeChangeRate(), index 2
                                              "Volumetric Strain (via Fe)", // scalar, element.getStrain(type=0)
                                              "Strain (DV, AP, AB)",        // display as 3-vector of scalars
                                              "Shear (xy, xz, yz)"};        // display as 3-vector of scalars

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

/**
 * @brief Widget that can change between the various displays for each of the element properties that can be chosen for display.
 * 
 */
class PropertyDisplayLayout : public QStackedLayout 
{
    Q_OBJECT
public:
    PropertyDisplayLayout();

public slots:
    /**
     * @brief Changes the element whose properties are currently being displayed
     * 
     * @param new_element Pointer to the new (selected) element. nullptr implies deselection has occurred.
     */
    void recieveNewElement(std::unique_ptr<ShapeBase> *new_element);
    /**
     * @brief Write the element properties to the file provided
     * 
     * @param filename Name of file to write to
     * @param write_headers Whether to break data up by writing section headers
     */
    void writeToFile(const QString &filename, bool write_headers=true);

signals:
    /**
     * @brief Clear all values in the managed displays.
     * 
     */
    void clearEntries();

private:
    std::unique_ptr<ShapeBase> *current_element;    // Pointer to the element properties we are displaying

    QGroupBox *defaultBox;                          // Box containing the default layout
    DefaultLayout *defaultDisplay;                  // Display for "<No property selected>"

    QGroupBox *growthBox;                           // Box containing the growth layout
    MatrixLayout3by3 *growthDisplay;                // Display for "Growth"

    QGroupBox *growthRateBox;                       // Box containing the growth rate layout
    VectorLayout3 *growthRateDisplay;               // Display for "Growth Rate"

    QGroupBox *intViscBox;                          // Box containing the internal viscosity layout
    SingleBoxLayout *intViscDisplay;                // Display for "Internal Viscosity"

    QGroupBox *youngModBox;                         // Box containing the Young's modulus layout
    SingleBoxLayout *youngModDisplay;               // Display for "Youngs modulus"

    QGroupBox *poissonBox;                          // Box containing the poisson ratio layout
    SingleBoxLayout *poissonDisplay;                // Display for "Poisson Ratio"

    QGroupBox *volumeGrowthBox;                     // Box containing the volume growth ratio layout
    SingleBoxLayout *volumeGrowthDisplay;           // Display for "Volume Growth"

    QGroupBox *shapeAndSizeBox;                     // Box containing the shape and size display
    SingleBoxLayout *shapeAndSizeDisplay;           // Display for "Emergant shape and size"

    QGroupBox *shapeChangeRateZBox;                 // Box containing the shape change rate (z) display
    SingleBoxLayout *shapeChangeRateZDisplay;       // Display for "Shape change rate (z)"

    QGroupBox *volStrainFeBox;                      // Box containing the volumetric strain (Fe) display
    SingleBoxLayout *volStrainFeDisplay;            // Display for volumetric strain (Fe)

    QGroupBox *strainBox;                           // Box containing the strain values display
    VectorLayout3 *strainDisplay;                   // 3-vector display for the DV/AP/AB strain
    QStringList strainLabels = {"DV", "AP", "AB"};  // component labels

    QGroupBox *shearBox;                            // Box containing the shear values display
    VectorLayout3 *shearDisplay;                    // 3-vector display for the shear in xy/yz/xz
    QStringList shearLabels = {"xy", "yz", "xz"};   // component labels

    /**
     * @brief Updates all displays in this panel by reading the properties of a new element
     * 
     */
    void readAndUpdateElementProperties();

    // strings used to denote sections in any written output file
    const QString section_header = "===== ";
    /**
     * @brief Creates the header of section_name to be written to an output file
     *
     * @param section_name The name of the section that is to be written next
     * @return QString The section header to be written to the output file
     */
    QString makeHeader(const QString &section_name);
};

# endif
