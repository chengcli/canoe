// application
#include <application/application.hpp>

// canoe
#include "variable.hpp"

void Variable::ConvertTo(Variable::Type type) {
  if (type == mytype_) {
    return;
  }

  if (type == Type::MassFrac) {
    convertToPrimitive();
  } else if (type == Type::MassConc) {
    convertToConserved();
  } else if (type == Type::MoleFrac) {
    convertToMoleFraction();
  } else if (type == Type::MoleConc) {
    convertToMoleConcentration();
  } else {
    throw RuntimeError("Variable", "Unknown variable type");
  }
}

void Variable::convertToPrimitive() {
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

void Variable::convertToConserved() {
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

void Variable::convertToMoleFraction() {
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

void Variable::convertToMoleConcentration() {
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
