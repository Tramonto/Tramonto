/*
  Prototypes from dft_mesh_lib.c
*/
int round_to_int(double x);
int offset_to_node(int *inode_ijk, int *offset_node,
                   int *reflect_flag);
void ijk_to_ijk_box(int *ijk, int *ijk_box);
void ijk_to_ijk_box_no_bound(int *ijk, int *ijk_box);
void ijk_box_to_ijk(int *ijk_box, int *ijk);
int offset_to_node_box(int *ijk_box, int *offset,
		       int *reflect_flag);
int element_to_node(int ielement);
void element_to_nodes(int ielement,int *nodes);
int node_to_elem(int inode_all, int local_node, int *reflect_flag);
int node_to_elem_return_dim(int inode_all, int local_node, int *reflect_flag,
			    int *idim_return, int *iside, int *periodic_flag);
int node_to_elem_no_reflect(int inode_all, int local_node);
int node_to_elem_v2(int inode_all, int local_node);
void node_to_position(int inode, double *NodePos);
int position_to_node(double *NodePos);
int map_0th_plane(int i, int Nplane);
void node_to_ijk(int node, int *ijk);
void node_box_to_ijk_box(int node_box, int *ijk_box);
int ijk_to_node(int *ijk);
int ijk_box_to_node_box(int *ijk_box);
int node_box_to_node(int inode_box);
int node_to_node_box(int inode);
int node_to_node_box_no_bound(int inode);
int unk_box_to_unk(int i_box);
int unk_to_unk_box(int i);
int el_box_to_el(int iel_box);
int element_box_to_node_box(int iel_box);
int el_to_el_box(int iel);
int node_box_to_elem_box(int inode_box, int local_node);
int node_box_to_elem_box_reflect(int inode_box, int local_node, int *reflect_flag);
