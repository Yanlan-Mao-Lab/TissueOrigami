# include "TissueFolding_GUI_elements.h"

using namespace std;

Label::Label() : QLabel() {
    setFont(DEF_FONT);
}
Label::Label(string text) : QLabel(text.c_str()) {
    setFont(DEF_FONT);
}

Header::Header() : QLabel() {
    setFont(DEF_HEADER_FONT);
}
Header::Header(string text) : QLabel(text.c_str()) {
    setFont(DEF_HEADER_FONT);
}

ReadOnlyBox::ReadOnlyBox() : QLineEdit() {
    // text that fills the box when nothing is provided
    setPlaceholderText("-");
    // make read-only
    setReadOnly(true);
    // use our default font
    setFont(DEF_FONT);
}
ReadOnlyBox::ReadOnlyBox(string placeholder_text) : QLineEdit() {
    // text that fills the box when nothing is provided
    setPlaceholderText(placeholder_text.c_str());
    // make read-only
    setReadOnly(true);
    // use our default font
    setFont(DEF_FONT);
}

SelectionBox::SelectionBox() : QLineEdit() {
    // text that fills the box when nothing is provided
    setPlaceholderText("-");
    // use our default font
    setFont(DEF_FONT);
}
SelectionBox::SelectionBox(string placeholder_text) : QLineEdit() {
    // text that fills the box when nothing is provided
    setPlaceholderText(placeholder_text.c_str());
    // use our default font
    setFont(DEF_FONT);
}