# include "TissueFolding_GUI_elements.h"
# include "ElementPropertiesUI.h"

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

ElementPropertiesUI::ElementPropertiesUI() {
    // add the "selected item properties" label to the grid
    addWidget(&selection_header, 0, 0, 1, 2, AL_LEFT);

    // add the "element name" label to the grid
    addWidget(&element_name_label, 1, 0, 1, 1, AL_LEFT);

    // add the "element name" box to the grid next to the label
    addWidget(&element_name_display, 1, 1, 1, 2, AL_LEFT);

    // add the node information headers (horz) into the grid
    for (int i=0; i<n_node_info_headers; i++) {
        addWidget(&node_info_labels_horz[i], 2, i+1, 1, 1, AL_CENTRE);
    }

    // add the node information node numbers (vert) to the grid
    for (int i=0; i<n_nodes_per_element; i++) {
        addWidget(&node_info_numbers_vert[i], 2+(i+1), 0, 1, 1, AL_LEFT);
    }

    // add the node information ID and coordinate boxes
    for (int i=0; i<n_coord_boxes; i++) {
        node_coord_boxes[i].setFixedWidth(def_fixed_box_width);
        // get the row and column indices in the grid for this box
        int row = 0; NodeInfoHeader col = NodeInfoHeader::ID;
        rowAndColOfBox(i, &row, &col);
        // add the widget to the grid
        // offset (row, col) by (3,1) to account for the selection boxes above,
        // and the node #x labels to the left
        addWidget(&node_coord_boxes[i], row+3, col+1, 1, 1, AL_LEFT);
    }

    // add node and element display labels to the grid
    addWidget(&node_selection_label, 0, 3, 1, 1, AL_LEFT);
    addWidget(&element_selection_label, 0, 4, 1, 1, AL_LEFT);
    // now add the corresponding boxes
    addWidget(&node_selection_box, 1, 3, 1, 1, AL_LEFT);
    addWidget(&element_selection_box, 1, 4, 1, 1, AL_LEFT);
};

void ElementPropertiesUI::setNodeSelectionValidator(int max_node_index, QObject *parent) {
    node_selection_box.initialseValidator(max_node_index, parent);
};

void ElementPropertiesUI::setElementSelectionValidator(int max_element_index, QObject *parent) {
    element_selection_box.initialseValidator(max_element_index, parent);
};

void ElementPropertiesUI::updateCoordBox(int box_number, QString text, bool set_enabled) {
    node_coord_boxes[box_number].setText(text);
    node_coord_boxes[box_number].setEnabled(set_enabled);
}
void ElementPropertiesUI::updateCoordBox(int row, NodeInfoHeader col, QString text, bool set_enabled)
{
    int box_number = getInfoBoxIndex(row, col);
    updateCoordBox(box_number, text, set_enabled);
}