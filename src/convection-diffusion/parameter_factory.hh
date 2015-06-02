// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Implementation of the 'factory method pattern' to create parameter classes on the fly
*/
#ifndef DUNE_PARAMETER_FACTORY_HH
#define DUNE_PARAMETER_FACTORY_HH
//
#include <string>
#include <iostream>
#include <list>
#include <map>
//
#include "parameter_base.hh"
#include "parameterA.hh"
#include "parameterB.hh"
#include "parameterC.hh"
#include "parameterD.hh"
#include "parameterE.hh"
#include "parameterF.hh"

//
// Generic factory that orchestrates product creation:
// The concrete product does not have to be known by the factory.
// We allow only one factory for one abstract parameter.
//
template<typename AbstractParameter,
     typename IdentifierType,
     typename ParameterCreator = AbstractParameter*(*)(const typename AbstractParameter::Traits::GridViewType&)>
class ParameterFactory
{

private:
  ParameterFactory(){
  };
  ParameterFactory(const ParameterFactory&);
  void operator=(const ParameterFactory &);

private:
  typedef std::map<IdentifierType,ParameterCreator> AssocMap;
  AssocMap _associations;

  typedef typename AbstractParameter::Traits::GridViewType GV;

public:
  static ParameterFactory& getInstance(){
    static ParameterFactory instance;
    return instance;
  }

  void registerParameter( const IdentifierType& docId, ParameterCreator createFN ){
    _associations[docId] = createFN;
  };

  bool unregisterParameter( const IdentifierType& docId ){
    return _associations.erase( docId ) == 1;
  };

  AbstractParameter* createParameter( const GV& gv, const IdentifierType& docId ){
    typename AssocMap::const_iterator it = _associations.find( docId );
    if( it == _associations.end() )
      return nullptr;
    else
      return (it->second)(gv);
  };



  template<typename GVType, typename RFType, typename IdType>
  static void registerAll(const GVType& gv){
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'A', createParameterA<GVType,RFType> );
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'B', createParameterB<GVType,RFType> );
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'C', createParameterC<GVType,RFType> );
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'D', createParameterD<GVType,RFType> );
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'E', createParameterE<GVType,RFType> );
    ParameterFactory<ParameterBase<GVType,RFType>,IdType>::getInstance().registerParameter( 'F', createParameterF<GVType,RFType> );
  }


};

#endif // DUNE_PARAMETER_FACTORY_HH
