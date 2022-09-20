# ifndef ELEMENT_PROPERTIES_UI_BOX_H
# define ELEMENT_PROPERTIES_UI_BOX_H

# include "TissueFolding_GUI_elements.h"

enum NodeInfoHeader {
    ID = 0,
    X = 1,
    Y = 2,
    Z = 3
};

const int n_node_info_headers = 4;
const int n_nodes_per_element = 6;
const int n_coord_boxes = n_node_info_headers * n_nodes_per_element;

/**
 * @brief Get the node number and column header of an info box
 * 
 * box_index = col*n_nodes_per_element + n_node_info_headers
 * 
 * @param[in] box_index The index of the info box
 * @param[out] node_num The node number (row ind))
 * @param[out] col The column header (col ind)
 */
void rowAndColOfBox(int box_index, int *node_num, NodeInfoHeader *col);
/**
 * @brief Get the index of the info box corresponding to a header and node number
 * 
 * @param node_number The node number
 * @param header The column header
 * @return int The index of the info box
 */
int infoBoxIndex(int node_number, NodeInfoHeader header);

/**
 * @brief Replaces the right-hand element properties box in the GUI
 *
 * Designed to replace the role of SelectionDisplayGrid
 */
class ElementPropertiesUI : public QGridLayout
{
public:
    /**
     * @brief Constructs the element-selection pannel.
     *
     */
    ElementPropertiesUI();

    // headers and labels
    Header selection_header = Header("Selected Item Properties");
    Header element_name_label = Header("Element name:");
    Header node_info_labels_horz[n_node_info_headers] = {Header("id"),
                                                         Header("x"),
                                                         Header("y"),
                                                         Header("z")};
    Header node_info_numbers_vert[n_nodes_per_element] = {Header("Node 0"),
                                                          Header("Node 1"),
                                                          Header("Node 2"),
                                                          Header("Node 3"),
                                                          Header("Node 4"),
                                                          Header("Node 5")};
    
    Label node_selection_label = Label("Select <br> node:");
    Label element_selection_label = Label("Select <br> element");

    // read-only boxes
    ReadOnlyBox element_name_display = ReadOnlyBox("No element selected");

    // node_info_boxes displays nodal information about the selected element
    ReadOnlyBox node_coord_boxes[n_coord_boxes];

    // input boxes, for manual (non-click) node and/or element selection
    SelectionBox node_selection_box;
    SelectionBox element_selection_box;

    void createNodeSelectionValidator(int max_node_index, QObject *parent = nullptr);
    void createElementSelectionValidator(int max_element_index, QObject *parent = nullptr);
};

#endif