// application
#include <application/exceptions.hpp>

// canoe
#include "variable.hpp"

void Variable::ConvertToPrimitive() {
  if (mytype_ == Type::MassFrac) {
    return;
  }

  if (mytype_ == Type::MassConc) {
    conservedToPrimitive();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToPrimitive();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToPrimitive();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  mytype_ = Type::MassFrac;
}

void Variable::ConvertToConserved() {
  if (mytype_ == Type::MassConc) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    primitiveToConserved();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToConserved();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToConserved();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  mytype_ = Type::MassConc;
}

void Variable::ConvertToMoleFraction() {
  if (mytype_ == Type::MoleFrac) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    primitiveToMoleFraction();
  } else if (mytype_ == Type::MassConc) {
    conservedToMoleFraction();
  } else if (mytype_ == Type::MoleConc) {
    moleConcentrationToMoleFraction();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  mytype_ = Type::MoleFrac;
}

void Variable::ConvertToMoleConcentration() {
  if (mytype_ == Type::MoleConc) {
    return;
  }

  if (mytype_ == Type::MassFrac) {
    primitiveToMoleConcentration();
  } else if (mytype_ == Type::MassConc) {
    conservedToMoleConcentration();
  } else if (mytype_ == Type::MoleFrac) {
    moleFractionToMoleConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }

  mytype_ = Type::MoleConc;
}
