# include "GUIBuildingBlocks.h"
# include "ElementBasicDisplay.h"
# include <string>

using namespace std;

void rowAndColOfBox(int box_index, int *node_num, NodeInfoHeader *col)
{
    if (box_index > n_coord_boxes || box_index < 0)
    {
        throw runtime_error("Error: requested out-of-bounds box_index");
    }
    // box_index = col * n_nodes_per_element + n_node_info_headers
    *node_num = box_index % n_nodes_per_element;
    *col = NodeInfoHeader(box_index / n_nodes_per_element);
}
int getInfoBoxIndex(int node_number, NodeInfoHeader header)
{
    // box_index = col*n_nodes_per_element + n_node_info_headers
    return ((int)header) * n_nodes_per_element + node_number;
}

ElementBasicDisplay::ElementBasicDisplay(QWidget *parent) : QGridLayout(parent)
{
    // add the "selected item properties" panel header to the grid
    // this will be the top-level parent of the QWidgets in this layout,
    // and answers to the parent of this layout
    selection_header = new Header(panel_header_text, parent);
    addWidget(selection_header, 0, 0, 1, 2, AL_LEFT);

    // add the element selection box
    element_selection_box = new SelectionBox(selection_header);
    addWidget(element_selection_box, 1, 4, 1, 1, AL_LEFT);
    // add the element selection box's label
    element_selection_label = new Label("Select <br> element", element_selection_box);
    addWidget(element_selection_label, 0, 4, 1, 1, AL_LEFT);

    // add node selection box
    node_selection_box = new SelectionBox(selection_header);
    addWidget(node_selection_box, 1, 3, 1, 1, AL_LEFT);
    // add node selection box label
    node_selection_label = new Label("Select <br> node:", node_selection_box);
    addWidget(node_selection_label, 0, 3, 1, 1, AL_LEFT);

    // add the "element name" label to the grid
    element_name_label = new Header(element_name_label_text, selection_header);
    addWidget(element_name_label, 1, 0, 1, 1, AL_LEFT);
    // add the "element name" box to the grid next to the label
    element_name_display = new ReadOnlyBox("No element selected", element_name_label);
    addWidget(element_name_display, 1, 1, 1, 2, AL_LEFT);

    // add the node information headers (horz) into the grid
    for (int i=0; i<n_node_info_headers; i++) {
        node_info_labels_horz[i] = new Header(horz_info_label_text[i], selection_header);
        addWidget(node_info_labels_horz[i], 2, i+1, 1, 1, AL_CENTRE);
    }

    // add the node information node numbers (vert) to the grid
    for (int i=0; i<n_nodes_per_element; i++) {
        node_info_numbers_vert[i] = new Header(vert_node_label_text + QString::number(i, 'f', 0), selection_header);
        addWidget(node_info_numbers_vert[i], 2+(i+1), 0, 1, 1, AL_LEFT);
    }

    // add the node information ID and coordinate boxes
    for (int i=0; i<n_coord_boxes; i++) {
        node_coord_boxes[i] = new ReadOnlyBox(selection_header);
        node_coord_boxes[i]->setFixedWidth(node_coord_box_width);
        // get the row and column indices in the grid for this box
        int row = 0; NodeInfoHeader col = NodeInfoHeader::ID;
        rowAndColOfBox(i, &row, &col);
        // add the widget to the grid
        // offset (row, col) by (3,1) to account for the selection boxes above,
        // and the node # labels to the left
        addWidget(node_coord_boxes[i], row+3, col+1, 1, 1, AL_LEFT);
    }
};

void ElementBasicDisplay::setNodeSelectionValidator(int max_node_index, QObject *parent)
{
    node_selection_box->initialseValidator(max_node_index, parent);
};
void ElementBasicDisplay::setElementSelectionValidator(int max_element_index, QObject *parent)
{
    element_selection_box->initialseValidator(max_element_index, parent);
};
void ElementBasicDisplay::setDisplayedElementName(const QString &text) {
    element_name_display->setText(text);
}

void ElementBasicDisplay::updateCoordBox(int box_number, QString text, bool set_enabled)
{
    node_coord_boxes[box_number]->setText(text);
    node_coord_boxes[box_number]->setEnabled(set_enabled);
}
void ElementBasicDisplay::updateCoordBox(int row, NodeInfoHeader col, QString text, bool set_enabled)
{
    int box_number = getInfoBoxIndex(row, col);
    updateCoordBox(box_number, text, set_enabled);
}

void ElementBasicDisplay::updateDisplayValues(std::unique_ptr<ShapeBase> *element) {
    // if we have not been passed a nullptr, we can safely update the information
    if (element) {
        // display the element's name
        QString new_name = QString::fromStdString((*element)->getName());
        setDisplayedElementName(new_name);
        // write the information to the infoboxes
        for(int row=0; row<n_nodes_per_element; row++) {
            // these are the box indices for this row
            int box_id_ind = getInfoBoxIndex(row, NodeInfoHeader::ID);
            int box_x_ind = getInfoBoxIndex(row, NodeInfoHeader::X);
            int box_y_ind = getInfoBoxIndex(row, NodeInfoHeader::Y);
            int box_z_ind = getInfoBoxIndex(row, NodeInfoHeader::Z);
            // extract ID of node
            int new_id = (*element)->NodeIds[row];
            // extract coordinates of node
            double x_pos = (*element)->Positions[row][0];
            double y_pos = (*element)->Positions[row][1];
            double z_pos = (*element)->Positions[row][2];
            // update information
            updateCoordBox(box_id_ind, QString::number(new_id, 'f', 0));
            updateCoordBox(box_x_ind, QString::number(x_pos, 'f', 2));
            updateCoordBox(box_y_ind, QString::number(y_pos, 'f', 2));
            updateCoordBox(box_z_ind, QString::number(z_pos, 'f', 2));
        }
    }
    // otherwise, interpret this as deselection and revert to "no display" state
    else {
        setDisplayedElementName("No element selected");
        for(int box_index=0; box_index<n_coord_boxes; box_index++) {
            updateCoordBox(box_index, "", false);
        }
    }
}

void ElementBasicDisplay::writeNodePositions(QString filename, std::unique_ptr<ShapeBase> *element) {
    // open the file, in append mode if necessary (might have asked for node positions to be written)
    ofstream file(filename.toStdString(), std::ios_base::app);
    file << "====== Node Information (id, x, y, z) =====\n";

    // now we just need to write the positions of the elements to the currently selected node
    // if we have not been passed a nullptr, we can safely update the information
    if (element) {
        // write the information to the infoboxes
        for (int row = 0; row < n_nodes_per_element; row++) {
            // these are the box indices for this row
            int box_id_ind = getInfoBoxIndex(row, NodeInfoHeader::ID);
            int box_x_ind = getInfoBoxIndex(row, NodeInfoHeader::X);
            int box_y_ind = getInfoBoxIndex(row, NodeInfoHeader::Y);
            int box_z_ind = getInfoBoxIndex(row, NodeInfoHeader::Z);
            // extract this node's information
            double node_id = (*element)->NodeIds[row];
            double x_pos = (*element)->Positions[row][0];
            double y_pos = (*element)->Positions[row][1];
            double z_pos = (*element)->Positions[row][2];
            // output information to file
            file << node_id << "," << x_pos << "," << y_pos << "," << z_pos << ",";
            file << "\n";
        }
    }
    // throw error if we recieved a bad element
    else {
        throw runtime_error("Error: got nullptr when attempting to write node positions of current element");
    }
    // close file now that we're done with it
    file.close();
}

void ElementBasicDisplay::resetElementSelection(int max_element_index, QWidget *parent) {
    // block signals whilst resetting
    element_selection_box->blockSignals(true);
    // reset the text in the selection box
    setElementSelectionValidator(max_element_index, parent);
    element_selection_box->setText("");
    // reopen to user input
    element_selection_box->blockSignals(false);
}

void ElementBasicDisplay::resetNodeSelection(int max_node_index, QWidget *parent) {
    // block signals whilst resetting
    node_selection_box->blockSignals(true);
    // reset the text in the selection box
    setNodeSelectionValidator(max_node_index, parent);
    node_selection_box->setText("");
    // reopen to user input
    node_selection_box->blockSignals(false);
}