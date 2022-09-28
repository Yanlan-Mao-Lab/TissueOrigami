# include "GUIBuildingBlocks.h"
using namespace std;

Label::Label(QWidget *parent) : QLabel(parent) {
    setFont(DEF_FONT);
}
Label::Label(string text, QWidget *parent) : QLabel(text.c_str(), parent) {
    setFont(DEF_FONT);
}
Label::Label(const QString &text, QWidget *parent) : QLabel(text, parent) {
    setFont(DEF_FONT);
}
Label::Label(const char *text, QWidget *parent) : QLabel(text, parent) {
    setFont(DEF_FONT);
}

Header::Header(QWidget *parent) : QLabel(parent) {
    setFont(DEF_HEADER_FONT);
}
Header::Header(string text, QWidget *parent) : QLabel(text.c_str(), parent) {
    setFont(DEF_HEADER_FONT);
}
Header::Header(const QString &text, QWidget *parent) : QLabel(text, parent) {
    setFont(DEF_HEADER_FONT);
}
Header::Header(const char *text, QWidget *parent) : QLabel(text, parent) {
    setFont(DEF_HEADER_FONT);
}

ReadOnlyBox::ReadOnlyBox(QWidget *parent) : QLineEdit(parent) {
    // text that fills the box when nothing is provided
    setPlaceholderText(default_text);
    // make read-only
    setReadOnly(true);
    // use our default font
    setFont(DEF_FONT);
}
ReadOnlyBox::ReadOnlyBox(string placeholder_text, QWidget *parent) : QLineEdit(parent) {
    // text that fills the box when nothing is provided
    setPlaceholderText(placeholder_text.c_str());
    // make read-only
    setReadOnly(true);
    // use our default font
    setFont(DEF_FONT);
}

SelectionBox::SelectionBox(QWidget *parent) : QLineEdit(parent) {
    // text that fills the box when nothing is provided
    setPlaceholderText("-");
    // use our default font
    setFont(DEF_FONT);
    // background colour is white
    setStyleSheet("background-color: white");
    // box should be of fixed width
    setFixedWidth(def_fixed_box_width);
}
SelectionBox::SelectionBox(string placeholder_text, QWidget *parent) : QLineEdit(parent) {
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

Button::Button(const QString &text, QWidget *parent) : QPushButton(text, parent) {
    // start disabled
    setEnabled(false);
}

VectorLayout3::VectorLayout3() {
    for(int index=0; index<n_comps; index++) {
        vector_components[index] = new ReadOnlyBox("-");
        // add each box corresponding to a component to the display, as a column vector
        vector_components[index]->setAlignment(AL_CENTRE);
        addWidget(vector_components[index], index, 0, 1, 3, AL_CENTRE);
    }
}
VectorLayout3::VectorLayout3(QStringList comp_labels) {
    for(int index=0; index<n_comps; index++) {
        // create the label and add to display
        component_labels[index] = new Header(comp_labels[index]);
        addWidget(component_labels[index], index, 0, 1, 1, AL_RIGHT);
        // add each box corresponding to a component to the display, as a column vector
        vector_components[index] = new ReadOnlyBox("-");
        vector_components[index]->setAlignment(AL_CENTRE);
        addWidget(vector_components[index], index, 1, 1, 2, AL_CENTRE);
    }
}
void VectorLayout3::setComponentValue(int index, const QString &text) {
    vector_components[index]->setText(text);
}
void VectorLayout3::clearAllComponents() {
    for(int index=0; index<n_comps; index++) {
        vector_components[index]->setText("-");
    }
}
void VectorLayout3::fillValues(double *values) {
    for(int index=0; index<n_comps; index++) {
        setComponentValue(index, QString::number(values[index]));
    }
}
void VectorLayout3::fillValues(float *values) {
    for(int index=0; index<n_comps; index++) {
        setComponentValue(index, QString::number(values[index]));
    }
}
void VectorLayout3::fillValues(std::array<double, 3> values) {
    for(int index=0; index<n_comps; index++) {
        setComponentValue(index, QString::number(values[index]));
    }
}
void VectorLayout3::fillValues(QStringList values) {
    for(int index=0; index<n_comps; index++) {
        setComponentValue(index, values[index]);
    }
}

MatrixLayout3by3::MatrixLayout3by3() {
    // add each box corresponding to a component to the display, in a grid
    for(int row=0; row<n_rows; row++) {
        for(int col=0; col<n_cols; col++) {
            int box_index = boxIndex(row, col);
            matrix_components[box_index] = new ReadOnlyBox("-");
            matrix_components[box_index]->setAlignment(AL_CENTRE);
            addWidget(matrix_components[box_index], row, col, 1, 1, AL_CENTRE);
        }
    }
}
void MatrixLayout3by3::setBoxValue(int row, int col, const QString &text) {
    setBoxValue(boxIndex(row, col), text);
}
void MatrixLayout3by3::setBoxValue(int box_index, const QString &text) {
    matrix_components[box_index]->setText(text);
}
void MatrixLayout3by3::clearAllValues() {
    for(int box_index=0; box_index<n_comps; box_index++) {
        matrix_components[box_index]->setText("-");
    }
}
void MatrixLayout3by3::fillValues(double *values) {
    for(int index=0; index<n_comps; index++) {
        QString value_as_string = QString::number(values[index]);
        setBoxValue(index, value_as_string);
    }
}
void MatrixLayout3by3::fillValues(QStringList values) {
    for(int index=0; index<n_comps; index++) {
        setBoxValue(index, values[index]);
    }   
}
void MatrixLayout3by3::fillValues(gsl_matrix *values) {
    for(int row=0; row<n_rows; row++) {
        for(int col=0; col<n_cols; col++) {
            int box_index = boxIndex(row, col);
            double component_value = gsl_matrix_get(values, row, col);
            setBoxValue(box_index, QString::number(component_value));
        }
    }
}
int MatrixLayout3by3::boxIndex(int row, int col) {
    if ((row<0) || (row>=n_rows) || (col<0) || (col>=n_cols)) {
        throw std::runtime_error("Growth matrix index out of bounds\n");
    }
    return row*n_cols + col;
}

SingleBoxLayout::SingleBoxLayout() : QGridLayout() {
    display_box = new ReadOnlyBox;
    display_box->setText(default_text);
    display_box->setAlignment(AL_CENTRE);
    addWidget(display_box);
}
void SingleBoxLayout::setDisplayValue(const QString &value) {
    display_box->setText(value);
}
void SingleBoxLayout::setDisplayValue(double value) {
    QString value_as_double = QString::number(value);
    setDisplayValue(value_as_double);
}
void SingleBoxLayout::clearValue() {
    display_box->setText(default_text);
}

#include<iostream>
CheckBox::CheckBox() : QCheckBox() {
    setFont(DEF_FONT);
    setCheckable(true);
    setChecked(false);
}
CheckBox::CheckBox(const QString &text, bool start_enabled) : QCheckBox(text) {
    setFont(DEF_FONT);
    setCheckable(true);
    setChecked(false);
    setEnabled(start_enabled);
}