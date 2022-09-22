# ifndef ELEMENT_BASIC_DISPLAY_H
# define ELEMENT_BASIC_DISPLAY_H

# include "GUIBuildingBlocks.h"
# include <QtWidgets>

/**
 * @brief Translates the headers of the node information into the column index they appear in.
 * 
 */
enum NodeInfoHeader {
    ID = 0,
    X = 1,
    Y = 2,
    Z = 3
};

// the number of information headers to display in the element information pannel (currently 4: id, x, y, z)
const int n_node_info_headers = 4;
// the number of nodes associated to each element in the mesh (6)
const int n_nodes_per_element = 6;
// the number of information boxes that will be required to display all node information pertaining to a particular element
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
int getInfoBoxIndex(int node_number, NodeInfoHeader header);

/**
 * @brief ElementBasicDisplay handles the display of basic element information, and selection of nodes and elements by passing their ID.
 *
 * This panel renders:
 * - The "selected item properties" header
 * - The node selection and element selection boxes
 * - The printout of node information associated to a given element
 */
class ElementBasicDisplay : public QGridLayout
{
    Q_OBJECT
public:
    /**
     * @brief Constructs the element-selection pannel.
     *
     */
    ElementBasicDisplay();

    // The "Selected Item Properties" header
    Header selection_header = Header("Selected Item Properties");
    // The label for the content of element_name_display
    Header element_name_label = Header("Element name:");
    // The (horizontal) headers of the information about each node that we wish to fetch and display
    Header node_info_labels_horz[n_node_info_headers] = {Header("id"),
                                                         Header("x"),
                                                         Header("y"),
                                                         Header("z")};
    // The (vertical) naming convention for the nodes associated to an element
    Header node_info_numbers_vert[n_nodes_per_element] = {Header("Node 0"),
                                                          Header("Node 1"),
                                                          Header("Node 2"),
                                                          Header("Node 3"),
                                                          Header("Node 4"),
                                                          Header("Node 5")};
    // The label for the node_selection_box
    Label node_selection_label = Label("Select <br> node:");
    // The label for the element_selection box
    Label element_selection_label = Label("Select <br> element");

    // Box to display the internal name of the selected element
    ReadOnlyBox element_name_display = ReadOnlyBox("No element selected");
    // Boxes that will display the information about the nodes associated to the selected element;
    // uses the convention box_index = col*n_nodes_per_element + n_node_info_headers
    ReadOnlyBox node_coord_boxes[n_coord_boxes];

    // input box for manual element selection by requesting a node that forms this element
    SelectionBox node_selection_box;
    // input box for manual element selection by element ID
    SelectionBox element_selection_box;

    /**
     * @brief Sets (or resets) the validator object for the node_selection_box
     * 
     * Inputs to the node selection must be a valid node ID, between 0 and max_node_index inclusive.
     * 
     * @param max_node_index The maximum valid node index
     * @param parent The parent object
     */
    void setNodeSelectionValidator(int max_node_index, QObject *parent = nullptr);
    /**
     * @brief Sets (or resets) the validator object for the element_selection_box
     *
     * Inputs to the element selection must be a valid element ID, between 0 and max_element_index inclusive.
     *
     * @param max_node_index The maximum valid element index
     * @param parent The parent object
     */
    void setElementSelectionValidator(int max_element_index, QObject *parent = nullptr);
    /**
     * @brief Updates the displayed content of a coordinate box; old content is overwritten.
     * 
     * @param box_number The box_index of the coordinate information box to update
     * @param text The (new) text to be displayed in the box
     * @param set_enabled Whether to render the box as enabled or disabled
     */
    void updateCoordBox(int box_number, QString text = "", bool set_enabled = true);
    /**
     * @brief Updates the displayed content of a coordinate box; old content is overwritten.
     *
     * @param row,col The row and column (header) in the pannel display of the box
     * @param text The (new) text to be displayed in the box
     * @param set_enabled Whether to render the box as enabled or disabled
     */
    void updateCoordBox(int row, NodeInfoHeader col, QString text="", bool set_enabled=true);

    /**
     * @brief Enables the dropdown menu and signals that the property box should also be updated.
     * 
     * Emits dropdownUpdate() if enabled is true.
     * 
     * @param enabled Whether to enable (true) or disable (false) the dropdown options
     */
    void enableDropdownSelection(bool enabled);

private:
    const int node_coord_box_width = 70;    // node coord box width
};

#endif