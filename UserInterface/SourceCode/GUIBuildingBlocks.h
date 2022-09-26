# ifndef GUI_BUILDING_BLOCKS_H
# define GUI_BUILDING_BLOCKS_H

# include <string>
# include <QtWidgets>
# include <gsl/gsl_matrix.h>

// default standard font for the GUI
const QFont DEF_FONT("SansSerif", 10);
// default bold-typset font for the GUI
const QFont DEF_HEADER_FONT("SaneSerif", 10, QFont::Bold, true);

// convenient alias for left-alignment
const Qt::AlignmentFlag AL_LEFT = Qt::AlignmentFlag::AlignLeft;
// convenient alias for centre-alignment
const Qt::AlignmentFlag AL_CENTRE = Qt::AlignmentFlag::AlignCenter;
// convenient alias for right-alignment
const Qt::AlignmentFlag AL_RIGHT = Qt::AlignmentFlag::AlignRight;

// default fixed size for infoboxes
const int def_fixed_box_width = 70;

/**
 * @brief Standard style for labels in the Qt interface.
 *
 * Labels are used to qualify information in infoboxes, radio buttons, and other (less-important) GUI objects that need labels but are not deserving of a header.
 */
class Label : public QLabel
{
public:
    Label();
    /**
     * @brief Construct a new Label object
     *
     * @param text Text of this label
     */
    Label(std::string text);
};

/**
 * @brief Standard style for headers in the Qt interface.
 * 
 * Headers are used to mark important information, provide names for pannels in the Qt interface, or qualify other important information. 
 */
class Header : public QLabel
{
public:
    Header();
    /**
     * @brief Construct a new Header object
     *
     * @param text Text of this header
     * @param parent Parent QWidget
     */
    Header(std::string text);
    Header(QString &text);
};

/**
 * @brief Style for read-only boxes in the Qt interface.
 * 
 * These boxes cannot be manually edited, and update their content only upon recieving signals from other components of the interface. 
 */
class ReadOnlyBox : public QLineEdit
{
public:
    ReadOnlyBox();
    /**
     * @brief Construct a new Read Only Box object
     *
     * @param placeholder_text Text to display when no input has been provided to the box
     */
    ReadOnlyBox(std::string placeholder_text);

private:
    const QString default_text = "-";       // default text when box is initalised
};

/**
 * @brief Style for user-input/selection boxes in the Qt interface.
 * 
 * These boxes allow the user to input a value manually, and usually send out a signal to other components when this happens. 
 */
class SelectionBox : public QLineEdit 
{
public:
    SelectionBox();
    /**
     * @brief Construct a new Selection Box object
     *
     * @param placeholder_text Text to display when no input has been provided to the box
     */
    SelectionBox(std::string placeholder_text);

    /**
     * @brief Setup the selection validator for this box.
     * 
     * Valid inputs to this box are integers between 0 and n_max (inclusive).
     * 
     * @param n_max The maximum permitted input value 
     * @param parent The parent Qt object
     */
    void initialseValidator(int n_max, QObject *parent = nullptr);
};

/**
 * @brief Style for dropdown menu options in the Qt interface
 * 
 */
class DropdownMenu : public QComboBox
{
public:
    DropdownMenu();
    /**
     * @brief Construct a new Dropdown Menu object
     * 
     * @param options The dropdown options that will be automatically added to the selection
     * @param n_opts The number of options provided
     */
    DropdownMenu(std::string *options, int n_opts);
    /**
     * @brief Construct a new Dropdown Menu object
     *
     * @param options The dropdown options that will be automatically added to the selection
     */
    DropdownMenu(QStringList options);
};

/**
 * @brief The layout of the space set aside for displaying element properties as a 3-vector with headers
 *
 */
class VectorLayout3 : public QGridLayout
{
    Q_OBJECT
public:
    VectorLayout3();
    /**
     * @brief Construct a new Vector Layout 3 object, with component labels
     * 
     * @param comp_labels List of labels for the components
     */
    VectorLayout3(QStringList comp_labels);

    /**
     * @brief Sets the displayed text of the component with the given index
     *
     * @param index The index of the component to set text for
     * @param text The text to set
     */
    void setComponentValue(int index, const QString &text);

public slots:
    /**
     * @brief Clears all values displayed in the component boxes
     *
     */
    void clearAllComponents();
    /**
     * @brief Set the values of each of the component boxes.
     *
     * @param values Values to set in the matrix display boxes
     */
    void fillValues(double *values);
    void fillValues(float *values);
    /**
     * @brief Set the values of each of the component boxes.
     *
     * @param values Values to set in the matrix display boxes
     */
    void fillValues(std::array<double, 3> values);
    /**
     * @brief Set the values of each of the component boxes.
     *
     * @param values Values to set in the matrix display boxes
     */
    void fillValues(QStringList values);

private:
    const int n_comps = 3;                        // number of vector components

    ReadOnlyBox *vector_components[3]; // Boxes that display the values of the components

    Header *component_labels[3]; // Labels for the component boxes
};

/**
 * @brief The layout of the space set aside for displaying element properties in a 3-by-3 matrix
 *
 */
class MatrixLayout3by3 : public QGridLayout
{
    Q_OBJECT
public:
    MatrixLayout3by3();

    /**
     * @brief Set the displayed text of the component at position (row, col)
     *
     * @param row,col Position of component in the matrix to set
     * @param text The text to display
     */
    void setBoxValue(int row, int col, const QString &text);
    /**
     * @brief Set the displayed text of the component whose internal index is box_index
     *
     * @param box_index Internal index of the component to set, box_index = MatrixLayout::BoxIndex(row, col)
     * @param text The text to display
     */
    void setBoxValue(int box_index, const QString &text);
    /**
     * @brief Retrieves the internal index for the component boxes, given a (row, col) index in the matrix.
     *
     * box_index = row * n_cols + col.
     *
     * @param row,col Position of the component in the matrix
     * @return int Internal index for the box that displays this components value
     */
    int boxIndex(int row, int col);

public slots:
    /**
     * @brief Clears all values in the display boxes
     *
     */
    void clearAllValues();
    /**
     * @brief Set the values of the component boxes.
     * 
     * Array should be 9 components long and values[index] will be the value set to matrix_components[index].
     * 
     * @param values Values to set in the matrix display boxes
     */
    void fillValues(double *values);
    void fillValues(QStringList values);
    /**
     * @brief Set the values of the component boxes, by providing a new matrix of values
     * 
     * @param values Values to set in the matrix display boxes
     */
    void fillValues(gsl_matrix *values);

private:
    const int n_cols = 3, n_rows = 3, n_comps = 9;

    ReadOnlyBox *matrix_components[9]; // The boxes that display the components of the growth matrix
};

/**
 * @brief A layout containing a single, centred box
 * 
 */
class SingleBoxLayout : public QGridLayout
{
    Q_OBJECT
public:
    SingleBoxLayout();

public slots:
    /**
     * @brief Set the value displayed in the display box
     * 
     * @param value The value to display
     */
    void setDisplayValue(const QString &value);
    void setDisplayValue(double value);
    /**
     * @brief Clear the value in the display box
     * 
     */
    void clearValue();

private:
    QString default_text = "-"; // default displayed text

    ReadOnlyBox *display_box; // box that displays the single value
};

# endif
