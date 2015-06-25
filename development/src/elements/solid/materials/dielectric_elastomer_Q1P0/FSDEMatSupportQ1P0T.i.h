
namespace Tahoe {

  //
  //
  //
  inline FSDEMatSupportQ1P0T::FSDEMatSupportQ1P0T(int ndof, int nip) :
    FSMatSupportT(ndof, nip), fE_List(0), fFSDielectricElastomerQ1P0(0)
  {

  }

  //
  //
  //
  inline const FSDielectricElastomerQ1P0T*
  FSDEMatSupportQ1P0T::FSDielectricElastomerQ1P0() const
  {
    return fFSDielectricElastomerQ1P0;
  }

  inline const dArrayT&
  FSDEMatSupportQ1P0T::ElectricField() const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[CurrIP()];
  }

  inline const dArrayT&
  FSDEMatSupportQ1P0T::ElectricField(int ip) const
  {
    if (fE_List == 0) {
      throw ExceptionT::kGeneralFail;
    }

    return (*fE_List)[ip];
  }

  inline void FSDEMatSupportQ1P0T::SetElectricField(
      const ArrayT<dArrayT>* E_List)
  {
    fE_List = E_List;
  }

  //
  // Return pointer to specified local array
  //
  inline const LocalArrayT*
  FSDEMatSupportQ1P0T::LocalArray(LocalArrayT::TypeT t) const
  {
    const LocalArrayT* pla = 0;

    switch (t) {

    case LocalArrayT::kESP:
      pla = fESP;
      break;

    default:
      //
      // Inherited
      //
      pla = FSMatSupportT::LocalArray(t);

    }

    return pla;

  }

  //
  // Set pointer to local array
  //
  inline void FSDEMatSupportQ1P0T::SetLocalArray(const LocalArrayT& array)
  {
    switch (array.Type()) {

    case LocalArrayT::kESP:
      fESP = &array;
      break;

    default:
      //
      // Inherited
      //
      FSMatSupportT::SetLocalArray(array);

    }

  }

  //
  // Nodal electric scalar potentials
  //
  inline const LocalArrayT*
  FSDEMatSupportQ1P0T::ScalarPotentials() const
  {
    return fESP;
  }

  //
  //
  //
  inline void FSDEMatSupportQ1P0T::SetScalarPotentials(
      const LocalArrayT& scalarPotentials)
  {
    fESP = &scalarPotentials;
  }


} //namespace Tahoe
