/*!
 * \file CTransLMVariable.cpp
 * \brief Definition of the solution fields.
 * \author A. Aranake, S. Kang
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../../include/variables/CTransIntermittencyVariable.hpp"

CTransIntermittencyVariable::CTransIntermittencyVariable(su2double Intermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CTurbVariable(npoint, ndim, nvar, config) {

  for(unsigned long iPoint=0; iPoint<nPoint; ++iPoint)
  {
    Solution(iPoint,0) = Intermittency;
  }

  Solution_Old = Solution;

  /*--- Setting CTransLMVariable of intermittency---*/
  IntermittencyVari.resize(nPoint) = Intermittency;
  zeta.resize(nPoint) = 0.0;
  sgn.resize(nPoint) = 0.0;
  extraTau.resize(nPoint) = 0.0;
  deflection.resize(nPoint) = 0.0;
  Mrel.resize(nPoint) = 0.0;
  TempVar1.resize(nPoint) = 0.0;
  TempVar2.resize(nPoint) = 0.0;
  TempVar3.resize(nPoint) = 0.0;
  TempVar4.resize(nPoint) = 0.0;
  TempVar5.resize(nPoint) = 0.0;
  TempVar6.resize(nPoint) = 0.0;

}

void CTransIntermittencyVariable::SetIntermittency_Fu_Func(unsigned long iPoint, su2double val_zeta, 
          su2double val_sgn, su2double val_extraTau, su2double val_deflection ) {

  zeta(iPoint) = val_zeta;
  sgn(iPoint) = val_sgn;
  extraTau(iPoint) = val_extraTau;
  deflection(iPoint) = val_deflection;
}

void CTransIntermittencyVariable::SetIntermittency_Wonder_Func(unsigned long iPoint, su2double var1, su2double var2 
          , su2double var3, su2double var4, su2double var5, su2double var6) {

  TempVar1(iPoint) = var1;
  TempVar2(iPoint) = var2;
  TempVar3(iPoint) = var3;
  TempVar4(iPoint) = var4;
  TempVar5(iPoint) = var5;
  TempVar6(iPoint) = var6;
  
}


void CTransIntermittencyVariable::SetIntermittency_Zhou_Func(unsigned long iPoint, su2double val_zeta, 
          su2double val_Mrel, su2double val_extraTau) {

  zeta(iPoint) = val_zeta;
  Mrel(iPoint) = val_Mrel;
  extraTau(iPoint) = val_extraTau;
}