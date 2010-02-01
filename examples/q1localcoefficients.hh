#include<dune/localfunctions/common/localkey.hh>
class Q1LocalCoefficients 
{
public:
  Q1LocalCoefficients () : li(4)  {
	for (int i=0; i<4; i++) li[i] = Dune::LocalKey(i,2,0);
  }
  
  //! number of coefficients
  int size () const { return 4; }
  
  //! map index i to local key
  const Dune::LocalKey& localKey (int i) const  {
	return li[i];
  } 
  
private:
  std::vector<Dune::LocalKey> li;
};
