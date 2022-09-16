# ifndef TISSUEFOLDING_UI_COMPONENTS.H
# define TISSUEFOLDING_UI_COMPONENTS.H

# include <string>
# include <QtWidgets>

const QFont DEF_FONT("SansSerif", 10);
const QFont DEF_HEADER_FONT("SaneSerif", 10, QFont::Bold, true);

const Qt::AlignmentFlag AL_LEFT = Qt::AlignmentFlag::AlignLeft;
const Qt::AlignmentFlag AL_CENTRE = Qt::AlignmentFlag::AlignCenter;

class Label : public QLabel
{
public:
    Label();
    /**
     * @brief Construct a new Label object
     *
     * This object sets the style for labels that appear in the GUI.
     *
     * @param text Text of this label
     */
    Label(std::string text);
};

/**
 * @brief Label style for headers in the Qt interface
 * 
 */
class Header : public QLabel
{
public:
    Header();
    /**
     * @brief Construct a new Header object
     *
     * This object sets the style for header-labels that appear in the GUI.
     *
     * @param text Text of this header
     */
    Header(std::string text);
};

/**
 * @brief Style for read-only boxes in the Qt simulation
 * 
 */
class ReadOnlyBox : public QLineEdit
{
public:
    ReadOnlyBox();
    /**
     * @brief Construct a new Read Only Box object
     *
     * This object sets the style for read-only boxes in the GUI interface
     *
     * @param placeholder_text Text to display when no input has been provided to the box
     */
    ReadOnlyBox(std::string placeholder_text);
};

class SelectionBox : public QLineEdit 
{
public:
    SelectionBox();
    /**
     * @brief Construct a new Selection Box object
     *
     * This object sets the style for read-only boxes in the GUI interface
     *
     * @param placeholder_text Text to display when no input has been provided to the box
     */
    SelectionBox(std::string placeholder_text);
};

# endif
