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
                                              "Internal Viscosity",         // scalar, Element.getInternalViscosity
                                              "Young Modulus",              // scalar, getYoungModulus
                                              "Poisson Ratio",              // scalar, Element.getPoissonRatio
                                              "Volume Growth",              // scalar, GrownVolume/ReferenceShape->Volume
                                              "Emergent Shape & Size",      // scalar, element.calculateEmergentShapeOrientation
                                              "ShapeChangeRate (z)",        // element.getShapeChangeRate(), then take index 2
                                              "Volumetric Strain (via Fe)", // scalar, element.getStrain(type=0)
                                              "Strain (DV, AP, AB)",        // each a scalar, getStrain(type1,2,3 respectively). Alternatively can be obtained from gsl_matrix_get(element.Strain,0/1/2,0)
                                              "Strian (xy, xz, yz)"};       // each a scalar, getStrain(type4,5,6 respectively). Alternatively can be obtained from gsl_matrix_get(element.Strain,3/4/5,0)

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
     * @param new_element Pointer to the new (selected) element
     */
    void recieveNewElement(std::unique_ptr<ShapeBase> *new_element);

signals:
    /**
     * @brief Clear all values in the managed displays.
     * 
     */
    void clearEntries();
    void updateEntries();

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

    QGroupBox *volumeGrowthBox;
    SingleBoxLayout *volumeGrowthDisplay;

    QGroupBox *shapeAndSizeBox;
    SingleBoxLayout *shapeAndSizeDisplay;

    QGroupBox *shapeChangeRateZBox;
    SingleBoxLayout *shapeChangeRateZDisplay;

    QGroupBox *volStrainFeBox;
    SingleBoxLayout *volStrainFeDisplay;

    QGroupBox *strainDV_AP_ABBox;
    VectorLayout3 *strainDV_AP_ABDisplay;
    QStringList DV_AP_AB_Labels = {"DV", "AP", "AB"};

    QGroupBox *strainPlanarDirsBox;
    VectorLayout3 *strainPlanarDirsDisplay;
    QStringList planarDirsLabels = {"xy", "xz", "yz"};

    void readAndUpdateElementProperties();
};

# endif
