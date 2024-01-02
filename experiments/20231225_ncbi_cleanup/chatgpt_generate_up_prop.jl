struct TreeNode
    label::String
    group::Union{Nothing,String}
    children::Vector{TreeNode}
end

function is_leaf(node::TreeNode)
    return isempty(node.children)
end

function set_group(node::TreeNode, group::String)
    node.group = group
end

function propagate_groups_up(node::TreeNode)
    if is_leaf(node)
        return
    end

    # Check if both children have groups assigned
    children_groups = Set(child.group for child in node.children if child.group !== nothing)

    if length(children_groups) == 1
        # Both children belong to the same group, set the current node to that group
        set_group(node, first(children_groups))
    end

    # Propagate the group information up to the parent
    parent = node
    while parent !== nothing
        # Check if all children of the parent have groups assigned
        parent_children_groups = Set(child.group for child in parent.children if child.group !== nothing)

        if length(parent_children_groups) == 1
            # All children of the parent belong to the same group, set the parent to that group
            set_group(parent, first(parent_children_groups))
        else
            # Stop propagating if not all children have the same group
            break
        end

        # Move up to the parent's parent
        parent = parent !== nothing ? parent.parent : nothing
    end
end

# Example usage
leaf1 = TreeNode("A1", "Group1", [])
leaf2 = TreeNode("A2", "Group1", [])
leaf3 = TreeNode("B1", "Group2", [])
leaf4 = TreeNode("B2", "Group2", [])
leaf5 = TreeNode("C1", "Group3", [])
leaf6 = TreeNode("C2", "Group3", [])

internal_node_A = TreeNode("A", nothing, [leaf1, leaf2])
internal_node_B = TreeNode("B", nothing, [leaf3, leaf4])
internal_node_C = TreeNode("C", nothing, [leaf5, leaf6])

root = TreeNode("Root", nothing, [internal_node_A, internal_node_B, internal_node_C])

# Propagate groups starting from the leaves
for leaf in [leaf1, leaf2, leaf3, leaf4, leaf5, leaf6]
    propagate_groups_up(leaf)
end


