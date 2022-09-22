# include "GUIBuildingBlocks.h"
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
    // background colour is white
    setStyleSheet("background-color: white");
    // box should be of fixed width
    setFixedWidth(def_fixed_box_width);
}
SelectionBox::SelectionBox(string placeholder_text) : QLineEdit() {
    // text that fills the box when nothing is provided
    setPlaceholderText(placeholder_text.c_str());
    // use our default font
    setFont(DEF_FONT);
}
void SelectionBox::initialseValidator(int n_max, QObject *parent) {
    setPlaceholderText( QString("# 0-%1").arg(n_max));
    setValidator( new QIntValidator(0, n_max, parent) );
}

DropdownMenu::DropdownMenu() : QComboBox() {
    // use our default font
    setFont(DEF_FONT);
    // start disabled (since an element needs to be selected before this option becomes available)
    setEnabled(false);
}
DropdownMenu::DropdownMenu(string *options, int n_opts) : QComboBox() {
    // use our default font
    setFont(DEF_FONT);
    // start disabled (since an element needs to be selected before this option becomes available)
    setEnabled(false);
    // add the options to the selection
    for(int i=0; i<n_opts; i++) {
        addItem(options[i].c_str());
    }
}
DropdownMenu::DropdownMenu(QStringList options) : QComboBox() {
    // use our default font
    setFont(DEF_FONT);
    // start disabled (since an element needs to be selected before this option becomes available)
    setEnabled(false);
    // add the options to the selection
    addItems(options);
}