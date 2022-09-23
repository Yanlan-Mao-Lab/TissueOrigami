# ifndef GUI_BUILDING_BLOCKS_H
# define GUI_BUILDING_BLOCKS_H

# include <string>
# include <QtWidgets>

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

/**
 * @brief Standard style for labels in the Qt interface.
 *
 * Labels are used to qualify information in infoboxes, radio buttons, and other (less-important) GUI objects that need labels but are not deserving of a header.
 */
class Label : public QLabel
{
    Q_OBJECT
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
    Q_OBJECT
public:
    Header();
    /**
     * @brief Construct a new Header object
     *
     * @param text Text of this header
     * @param parent Parent QWidget
     */
    Header(std::string text);
};

/**
 * @brief Style for read-only boxes in the Qt interface.
 * 
 * These boxes cannot be manually edited, and update their content only upon recieving signals from other components of the interface. 
 */
class ReadOnlyBox : public QLineEdit
{
    Q_OBJECT
public:
    ReadOnlyBox();
    /**
     * @brief Construct a new Read Only Box object
     *
     * @param placeholder_text Text to display when no input has been provided to the box
     */
    ReadOnlyBox(std::string placeholder_text);
    
private:
    // default fixed size for infoboxes
    const int def_fixed_box_width = 70;
};

/**
 * @brief Style for user-input/selection boxes in the Qt interface.
 * 
 * These boxes allow the user to input a value manually, and usually send out a signal to other components when this happens. 
 */
class SelectionBox : public QLineEdit 
{
    Q_OBJECT
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

private:
    // default fixed size for infoboxes
    const int def_fixed_box_width = 70;
};

/**
 * @brief Style for dropdown menu options in the Qt interface
 * 
 */
class DropdownMenu : public QComboBox
{
    Q_OBJECT
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
 * @brief Style for clickable buttons in the GUI interface
 * 
 */
class Button : public QPushButton 
{
    Q_OBJECT
public:
    Button();
    /**
     * @brief Construct a new Button object
     * 
     * @param text The text to display on the button
     * @param parent The parent widget
     */
    Button(const QString &text, QWidget *parent=nullptr);
private:
    const int default_button_width = 60;    // default width for clickable buttons
    const int default_button_height = 10;   // default height for clickable buttons
};

# endif