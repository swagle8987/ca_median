from dendropy import Tree,Node,Taxon
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

    def regenerate_taxon(self):
        for v in self.leaf_nodes():
            v.taxon = Taxon()
            v.taxon.label = v.label
        self.reconstruct_taxon_namespace()

