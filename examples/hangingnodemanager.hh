#ifndef HANGINGNODEMANAGER_HH
#define HANGINGNODEMANAGER_HH

#include<dune/grid/common/grid.hh>
#include<dune/common/float_cmp.hh>

template<class Grid>
class HangingNodeManager
{
private:
  enum{ verbosity = 3 };
  // This will hold the information on whether a vertex is hanging. Is
  // is accessed via the local index of the vertex in the reference
  // element. 
  std::vector<bool> is_hanging;

  // Codim 0 Mapper
  template<int dim>
  struct Codim0Layout {
    bool contains (Dune::GeometryType gt) const {
      if(gt.dim() == dim)
	return true;
      return false;
    }
  };

  // Codim dim Mapper
  template<int dim>
  struct CodimDimLayout {
    bool contains (Dune::GeometryType gt) const {
      if(gt.dim() == 0)
	return true;
      return false;
    }
  };
 

  class NodeInfo
  {
  public:
    // Minimum level of elements containing this node
    unsigned short minimum_level;

    // Maximum level of elements containing this node
    unsigned short maximum_level;
    
    // Minimum level of elements touching this node
    unsigned short minimum_touching_level;

    NodeInfo() : minimum_level(std::numeric_limits<unsigned short>::max()), 
		 maximum_level(0),
		 minimum_touching_level(std::numeric_limits<unsigned short>::max())
    {}
  };

    std::vector<NodeInfo> node_info;

public:
  typedef typename Grid::LeafGridView GridView;
  enum {dim = GridView::dimension};
  typedef typename GridView::template Codim<0>::EntityPointer CellEntityPointer;
  typedef typename GridView::template Codim<dim>::EntityPointer VertexEntityPointer;
  typedef typename GridView::template Codim<0>::Iterator Iterator;
  typedef typename GridView::IntersectionIterator IntersectionIterator;
  typedef typename GridView::Grid::ctype ctype;
  typedef typename Dune::FieldVector<ctype,dim> Point;
  
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,Codim0Layout> CellMapper;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,CodimDimLayout> VertexMapper;

  Grid & grid;
  CellMapper cell_mapper;
  VertexMapper vertex_mapper;

public:

  void analyzeView()
  {
    cell_mapper.update();
    vertex_mapper.update();

    node_info = std::vector<NodeInfo>(vertex_mapper.size());

    const GridView & gv = grid.leafView();

    Iterator it = gv.template begin<0>();
    Iterator eit = gv.template end<0>();

    for(;it!=eit;++it){

      const Dune::GenericReferenceElement<double,dim> & 
	reference_element = 
	Dune::GenericReferenceElements<double,dim>::general(it->geometry().type()); 
      
      // level of this element
      const unsigned short level = it->level();
      
      // number of vertices in this element
      const size_t v_size = reference_element.size(dim);

      // update minimum_level and maximum_level for vertices in this
      // cell
      for(size_t i=0; i<v_size; ++i){
	const VertexEntityPointer vertex = it->template subEntity<dim>(i);
	const size_t v_globalindex = vertex_mapper.map( *vertex );
	unsigned short & min = node_info[v_globalindex].minimum_level;
	unsigned short & max = node_info[v_globalindex].maximum_level;
	if (level < min) min = level;
	if (level > max) max = level;
      }

      // Now we still have to update minimum_touching_level for this
      // cell
      
      IntersectionIterator fit = gv.ibegin(*it);
      IntersectionIterator efit = gv.iend(*it);

      for(;fit!=efit;++fit){
	if(!(*fit).boundary()){

	  const int eLocalIndex =  fit->indexInInside();
	  const int fLocalIndex =  fit->indexInOutside();

	  const int e_level = fit->inside()->level();
	  const int f_level = fit->outside()->level();
	  
	  // a conforming face has no hanging nodes
	  if(fit->conforming())
	    continue;

	  // so far no support for initially non conforming grids
	  assert(e_level != f_level);

	  // this check needs to be performed on the element containing
	  // the hanging node only
	  if(e_level < f_level)
	    continue;

	  // numbero of vertices in face
	  const int e_v_size = reference_element.size(eLocalIndex,1,dim);
            
	  for(int i=0; i<e_v_size;++i){
	    const int e_v_index = reference_element.subEntity(eLocalIndex,1,i,dim);
	    const VertexEntityPointer vertex = it->template subEntity<dim>(e_v_index);
	    const size_t v_globalindex = vertex_mapper.map( *vertex );
	    unsigned short & min = node_info[v_globalindex].minimum_touching_level;
	    if( f_level < min) min = f_level;
	  }
	}
      }
    }
  }

  HangingNodeManager(Grid & _grid)
    : grid(_grid),
      cell_mapper(grid.leafView()),
      vertex_mapper(grid.leafView())
  { analyzeView(); }

  const std::vector<bool> & hangingNodes(CellEntityPointer & e)
  {
      const Dune::GenericReferenceElement<double,dim> & 
	reference_element = 
	Dune::GenericReferenceElements<double,dim>::general(e->geometry().type()); 

      // number of vertices in this element
      const size_t v_size = reference_element.size(dim);

      // make sure the return array is big enough
      is_hanging.resize(v_size);

      // update minimum_level and maximum_level for vertices in this
      // cell
      for(size_t i=0; i<v_size; ++i){
	const VertexEntityPointer & vertex = e->template subEntity<dim>(i);
	const size_t v_globalindex = vertex_mapper.map( *vertex );

	// here we make use of the fact that a node is hanging if and
	// only if it touches a cell of a level smaller than the
	// smallest level of all element containing the node
	const NodeInfo & v_info = node_info[v_globalindex];
	if(v_info.minimum_touching_level < v_info.minimum_level){
	  is_hanging[i] = true;
#ifndef NDEBUG
	  if(verbosity){
	    const Point & local  = reference_element.position(i,dim);
	    const Point global = e->geometry().global(local);
	    cout << "Found hanging node with id " << v_globalindex << " at " << global << endl;
	  }
#endif
	}
	else
	  is_hanging[i] = false;
      }
  }

  void adaptToIsolatedHangingNodes()
  {
    bool reiterate;
    do{
      reiterate = false;
      
      const GridView & gv = grid.leafView();

      Iterator it = gv.template begin<0>();
      Iterator eit = gv.template end<0>();

      for(;it!=eit;++it){

	const Dune::GenericReferenceElement<double,dim> & 
	  reference_element = 
	  Dune::GenericReferenceElements<double,dim>::general(it->geometry().type()); 

	const unsigned short level = it->level();
	// number of vertices in this element
	const size_t v_size = reference_element.size(dim);

	// make sure the return array is big enough
	is_hanging.resize(v_size);

	// update minimum_level and maximum_level for vertices in this
	// cell
	for(size_t i=0; i<v_size; ++i){
	  const VertexEntityPointer & vertex = it->template subEntity<dim>(i);
	  const size_t v_globalindex = vertex_mapper.map( *vertex );
	  const NodeInfo & v_info = node_info[v_globalindex];

	  const unsigned short level_diff = v_info.maximum_level - level;
	  if( level_diff > 1){
	    grid.mark(1, *it);
	    reiterate = true;
	    if(verbosity){
	      cout << "Refining element nr " << element_mapper.map(*it) 
		   << " to isolate hanging nodes " << endl;
	    }
	    break;
	  }
	} // i

      } // it

      if(reiterate){
	grid.preAdapt();
	grid.adapt();
	grid.postAdapt();
	analyzeView();
      }

    }while(reiterate);
  }

};

#endif 
