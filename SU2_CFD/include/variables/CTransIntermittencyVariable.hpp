/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "CTurbVariable.hpp"

/*!
 * \class CTransLMVariable
 * \brief Transition model variables.
 * \ingroup Turbulence_Model
 * \author A. Bueno, S. Kang.
 */

class CTransIntermittencyVariable final : public CTurbVariable {
protected:
  VectorType IntermittencyVari, zeta, sgn, extraTau, deflection, Mrel;
  VectorType TempVar1, TempVar2, TempVar3, TempVar4, TempVar5, TempVar6;
  
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] Intermittency - intermittency(gamma) (initialization value).
   * \param[in] ReThetaT - momentum thickness Reynolds number(ReThetaT)(initialization value).
   * \param[in] gammaSep - separation intermittency(gamma) (initialization value).
   * \param[in] gammaEff - effective intermittency(gamma) (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransIntermittencyVariable(su2double Intermittency, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransIntermittencyVariable() override = default;

  void SetIntermittency_Fu_Func(unsigned long iPoint, su2double zeta, su2double sgn, su2double extraTau, su2double deflection ) override;

  void SetIntermittency_Wonder_Func(unsigned long iPoint, su2double tempVar1, su2double tempVar2, su2double tempVar3, su2double tempVar4, su2double tempVar5, su2double tempVar6 ) override;

  void SetIntermittency_Zhou_Func(unsigned long iPoint, su2double zeta, su2double Mrel,su2double extraTau) override;

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Fu_Func_zeta(unsigned long iPoint) const override { return zeta(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Fu_Func_sgn(unsigned long iPoint) const override { return sgn(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Fu_Func_extraTau(unsigned long iPoint) const override { return extraTau(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Fu_Func_deflection(unsigned long iPoint) const override { return deflection(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Zhou_Func_zeta(unsigned long iPoint) const override { return zeta(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Zhou_Func_Mrel(unsigned long iPoint) const override { return Mrel(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermit_Zhou_Func_extraTau(unsigned long iPoint) const override { return extraTau(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var1(unsigned long iPoint) const override { return TempVar1(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var2(unsigned long iPoint) const override { return TempVar2(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var3(unsigned long iPoint) const override { return TempVar3(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var4(unsigned long iPoint) const override { return TempVar4(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var5(unsigned long iPoint) const override { return TempVar5(iPoint); }

  /*!
   * \brief Value of Wonder Variable.
   */
  inline su2double GetIntermit_Wonder_Func_var6(unsigned long iPoint) const override { return TempVar6(iPoint); }



};
