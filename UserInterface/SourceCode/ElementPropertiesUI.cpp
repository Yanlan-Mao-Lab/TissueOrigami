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
int getInfoBoxIndex(NodeInfoHeader header, int node_number)
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
        addWidget(&node_info_numbers_vert[i], 2+(i+1), 1, 1, AL_LEFT);
    }

    // add the node information ID and coordinate boxes
    for (int i=0; i<n_coord_boxes; i++) {
        node_coord_boxes[i].setFixedWidth(coord_box_fixed_width);
        // get the row and column indices in the grid for this box
        int *row; NodeInfoHeader *col;
        rowAndColOfBox(i, row, col);
        // add the widget to the grid
        addWidget(&node_coord_boxes[i], *row, *col, 1, 1, AL_LEFT);
    }

    // add the node selection label and box
    addWidget(&node_selection_label, 0, 3, 1, 1, AL_LEFT);

    //TSTK GOT TO HERE BEFORE ENOUGH QT FOR THE DAY. RESUME ON LINE 142 of MainWindow.cpp, the selection boxes are more complex and have our first signals and switches
    
    // add the element selection label and box
    addWidget(&element_selection_label, 0, 4, 1, 1, AL_LEFT);

};
