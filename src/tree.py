from dendropy import Tree,Node
class CustomNode(Node):
    
    def __init__(self,id=0,**kwargs):
        super().__init__(**kwargs)
        self.id = id

    def new_node(self,**kwargs):
        n = self.__class__(0,**kwargs)
        return self.add_child(node=n)

class CustomTree(Tree):
    counter = 0
    def node_factory(cls, **kwargs):
        cls.counter += 1
        return CustomNode(cls.counter,**kwargs)

    def __init__(self,*args,**kwargs):
        counter = 0
        super().__init__(*args,**kwargs)

    def reindex(self):
        for i,v in enumerate(self.nodes()):
            v.id = i+self.counter


"""

So when we do spr we do the following:
1. copy_tree = t1.extract_tree(node_factory=lambda x:return CustomNode(x))
2. spr(tree,prune_node,attach_node):
    
2. new_node = Custom
    p_n = tree.find_node(lambda x:x.id==prune_node.id)
    tree.prune_subtree(p_n)
    attach_above_node(tree,p_n,attach_node)
"""
