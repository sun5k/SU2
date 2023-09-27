/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
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
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"
#include "./trans_correlations.hpp"

/*!
 * \class CSourcePieceWise_TranLM
 * \brief Class for integrating the source terms of the LM transition model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  const LM_ParsedOptions options;

  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double c_CF = 0.6;
  const su2double sigmat = 2.0;

  TURB_FAMILY TurbFamily;
  su2double hRoughness;

  su2double IntermittencySep = 1.0;
  su2double IntermittencyEff = 1.0;

  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];  // Static storage for the Jacobian (which needs to be pointer for return type).

  TransLMCorrelations TransCorrelations;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 2, config), idx(val_nDim, config->GetnSpecies()), options(config->GetLMParsedOptions()){
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;

    TurbFamily = TurbModelFamily(config->GetKind_Turb_Model());

    hRoughness = config->GethRoughness();

    TransCorrelations.SetOptions(options);

  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma, and TransVar[0] = ReThetaT ---*/
    /*--- dU/dx = PrimVar_Grad[1][0] ---*/
    AD::StartPreacc();
    AD::SetPreaccIn(StrainMag_i);
    AD::SetPreaccIn(ScalarVar_i, nVar);
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(TransVar_i, nVar);
    AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(Volume);
    AD::SetPreaccIn(dist_i);
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);

    su2double VorticityMag =
        sqrt(Vorticity_i[0] * Vorticity_i[0] + Vorticity_i[1] * Vorticity_i[1] + Vorticity_i[2] * Vorticity_i[2]);

    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1 + idx.Velocity()];
    const su2double vel_w = (nDim == 3) ? V_i[2 + idx.Velocity()] : 0.0;

    const su2double Velocity_Mag = sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w);

    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
    Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

    Residual[0] = 0.0;
    Residual[1] = 0.0;
    Jacobian_i[0][0] = 0.0;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = 0.0;

    if (dist_i > 1e-10) {
      su2double Tu = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) Tu = max(100.0 * sqrt(2.0 * ScalarVar_i[0] / 3.0) / Velocity_Mag, 0.027);
      if (TurbFamily == TURB_FAMILY::SA) Tu = config->GetTurbulenceIntensity_FreeStream() * 100;

      /*--- Corr_RetC correlation*/
      const su2double Corr_Rec = TransCorrelations.ReThetaC_Correlations(Tu, TransVar_i[1]);

      /*--- F_length correlation*/
      const su2double Corr_F_length = TransCorrelations.FLength_Correlations(Tu, TransVar_i[1]);

      /*--- F_length ---*/
      su2double F_length = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double r_omega = Density_i * dist_i * dist_i * ScalarVar_i[1] / Laminar_Viscosity_i;
        const su2double f_sub = exp(-pow(r_omega / 200.0, 2));
        F_length = Corr_F_length * (1. - f_sub) + 40.0 * f_sub;
      }
      if (TurbFamily == TURB_FAMILY::SA) F_length = Corr_F_length;

      /*--- F_onset ---*/
      su2double R_t = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) R_t = Density_i * ScalarVar_i[0] / Laminar_Viscosity_i / ScalarVar_i[1];
      if (TurbFamily == TURB_FAMILY::SA) R_t = Eddy_Viscosity_i / Laminar_Viscosity_i;

      const su2double Re_v = Density_i * dist_i * dist_i * StrainMag_i / Laminar_Viscosity_i;
      const su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
      su2double F_onset2 = 1.0;
      su2double F_onset3 = 1.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
        F_onset3 = max(1.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      if (TurbFamily == TURB_FAMILY::SA) {
        F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 4.0);
        F_onset3 = max(2.0 - pow(R_t / 2.5, 3.0), 0.0);
      }
      const su2double F_onset = max(F_onset2 - F_onset3, 0.0);

      /*-- Gradient of velocity magnitude ---*/

      su2double dU_dx = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][0] + 2. * vel_v * PrimVar_Grad_i[2][0]);
      if (nDim == 3) dU_dx += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][0]);

      su2double dU_dy = 0.5 / Velocity_Mag * (2. * vel_u * PrimVar_Grad_i[1][1] + 2. * vel_v * PrimVar_Grad_i[2][1]);
      if (nDim == 3) dU_dy += 0.5 / Velocity_Mag * (2. * vel_w * PrimVar_Grad_i[3][1]);

      su2double dU_dz = 0.0;
      if (nDim == 3)
        dU_dz =
            0.5 / Velocity_Mag *
            (2. * vel_u * PrimVar_Grad_i[1][2] + 2. * vel_v * PrimVar_Grad_i[2][2] + 2. * vel_w * PrimVar_Grad_i[3][2]);

      su2double du_ds = vel_u / Velocity_Mag * dU_dx + vel_v / Velocity_Mag * dU_dy;
      if (nDim == 3) du_ds += vel_w / Velocity_Mag * dU_dz;

      /*-- Calculate blending function f_theta --*/
      su2double time_scale = 500.0 * Laminar_Viscosity_i / Density_i / Velocity_Mag / Velocity_Mag;
      if (options.LM2015)
        time_scale = min(time_scale,
                         Density_i * LocalGridLength_i * LocalGridLength_i / (Laminar_Viscosity_i + Eddy_Viscosity_i));
      const su2double theta_bl = TransVar_i[1] * Laminar_Viscosity_i / Density_i / Velocity_Mag;
      const su2double delta_bl = 7.5 * theta_bl;
      const su2double delta = 50.0 * VorticityMag * dist_i / Velocity_Mag * delta_bl + 1e-20;

      su2double f_wake = 0.0;
      if (TurbFamily == TURB_FAMILY::KW) {
        const su2double re_omega = Density_i * ScalarVar_i[1] * dist_i * dist_i / Laminar_Viscosity_i;
        f_wake = exp(-pow(re_omega / (1.0e+05), 2));
      }
      if (TurbFamily == TURB_FAMILY::SA) f_wake = 1.0;

      const su2double var1 = (TransVar_i[0] - 1.0 / c_e2) / (1.0 - 1.0 / c_e2);
      const su2double var2 = 1.0 - pow(var1, 2.0);
      const su2double f_theta = min(max(f_wake * exp(-pow(dist_i / delta, 4)), var2), 1.0);
      const su2double f_turb = exp(-pow(R_t / 4, 4));

      su2double f_theta_2 = 0.0;
      if (options.LM2015)
        f_theta_2 = min(f_wake * exp(-pow(dist_i / delta, 4.0)), 1.0);

      /*--- Corr_Ret correlation*/
      const su2double Corr_Ret_lim = 20.0;
      su2double f_lambda = 1.0;

      su2double Retheta_Error = 200.0, Retheta_old = 0.0;
      su2double lambda = 0.0;
      su2double Corr_Ret = 20.0;

      for (int iter = 0; iter < 100; iter++) {
        su2double theta = Corr_Ret * Laminar_Viscosity_i / Density_i / Velocity_Mag;
        lambda = Density_i * theta * theta / Laminar_Viscosity_i * du_ds;
        lambda = min(max(-0.1, lambda), 0.1);

        if (lambda <= 0.0) {
          f_lambda = 1. - (-12.986 * lambda - 123.66 * lambda * lambda - 405.689 * lambda * lambda * lambda) *
                              exp(-pow(Tu / 1.5, 1.5));
        } else {
          f_lambda = 1. + 0.275 * (1. - exp(-35. * lambda)) * exp(-Tu / 0.5);
        }

        if (Tu <= 1.3) {
          Corr_Ret = f_lambda * (1173.51 - 589.428 * Tu + 0.2196 / Tu / Tu);
        } else {
          Corr_Ret = 331.5 * f_lambda * pow(Tu - 0.5658, -0.671);
        }
        Corr_Ret = max(Corr_Ret, Corr_Ret_lim);

        Retheta_Error = fabs(Retheta_old - Corr_Ret) / Retheta_old;

        if (Retheta_Error < 0.0000001) {
          break;
        }

        Retheta_old = Corr_Ret;
      }

      /*-- Corr_RetT_SCF Correlations--*/
      su2double ReThetat_SCF = 0.0;
      if (options.LM2015) {
        su2double VelocityNormalized[3];
        VelocityNormalized[0] = vel_u / Velocity_Mag;
        VelocityNormalized[1] = vel_v / Velocity_Mag;
        if (nDim == 3) VelocityNormalized[2] = vel_w / Velocity_Mag;

        su2double StreamwiseVort = 0.0;
        for (auto iDim = 0u; iDim < nDim; iDim++) {
          StreamwiseVort += VelocityNormalized[iDim] * Vorticity_i[iDim];
        }
        StreamwiseVort = abs(StreamwiseVort);

        const su2double H_CF = StreamwiseVort * dist_i / Velocity_Mag;
        const su2double DeltaH_CF = H_CF * (1.0 + min(Eddy_Viscosity_i / Laminar_Viscosity_i, 0.4));
        const su2double DeltaH_CF_Minus = max(-1.0 * (0.1066 - DeltaH_CF), 0.0);
        const su2double DeltaH_CF_Plus = max(0.1066 - DeltaH_CF, 0.0);
        const su2double fDeltaH_CF_Minus = 75.0 * tanh(DeltaH_CF_Minus / 0.0125);
        const su2double fDeltaH_CF_Plus = 6200 * DeltaH_CF_Plus + 50000 * DeltaH_CF_Plus * DeltaH_CF_Plus;

        const su2double toll = 1e-5;
        su2double error = toll + 1.0;
        su2double thetat_SCF = 0.0;
        su2double rethetat_SCF_old = 20.0;
        const int nMax = 100;

        int iter;
        for (iter = 0; iter < nMax && error > toll; iter++) {
          thetat_SCF = rethetat_SCF_old * Laminar_Viscosity_i / (Density_i * (Velocity_Mag / 0.82));
          thetat_SCF = max(1e-20, thetat_SCF);

          ReThetat_SCF = -35.088 * log(hRoughness / thetat_SCF) + 319.51 + fDeltaH_CF_Plus - fDeltaH_CF_Minus;

          error = abs(ReThetat_SCF - rethetat_SCF_old) / rethetat_SCF_old;

          rethetat_SCF_old = ReThetat_SCF;
        }
      }

      /*-- production term of Intermeittency(Gamma) --*/
      const su2double Pg =
          F_length * c_a1 * Density_i * StrainMag_i * sqrt(F_onset * TransVar_i[0]) * (1.0 - c_e1 * TransVar_i[0]);

      /*-- destruction term of Intermeittency(Gamma) --*/
      const su2double Dg = c_a2 * Density_i * VorticityMag * TransVar_i[0] * f_turb * (c_e2 * TransVar_i[0] - 1.0);

      /*-- production term of ReThetaT --*/
      const su2double PRethetat = c_theta * Density_i / time_scale * (Corr_Ret - TransVar_i[1]) * (1.0 - f_theta);

      /*-- destruction term of ReThetaT --*/
      // It should not be with the minus sign but I put for consistency
      su2double DRethetat = 0.0;
      if (options.LM2015)
        DRethetat = -c_theta * (Density_i / time_scale) * c_CF * min(ReThetat_SCF - TransVar_i[1], 0.0) * f_theta_2;

      /*--- Source ---*/
      Residual[0] += (Pg - Dg) * Volume;
      Residual[1] += (PRethetat - DRethetat) * Volume;

      /*--- Implicit part ---*/
      Jacobian_i[0][0] = (F_length * c_a1 * StrainMag_i * sqrt(F_onset) *
                              (0.5 * pow(TransVar_i[0], -0.5) - 1.5 * c_e1 * pow(TransVar_i[0], 0.5)) -
                          c_a2 * VorticityMag * f_turb * (2.0 * c_e2 * TransVar_i[0] - 1.0)) *
                         Volume;
      Jacobian_i[0][1] = 0.0;
      Jacobian_i[1][0] = 0.0;
      Jacobian_i[1][1] = -c_theta / time_scale * (1.0 - f_theta) * Volume;
      if (options.LM2015 && ReThetat_SCF - TransVar_i[1] < 0)
        Jacobian_i[1][1] += (c_theta / time_scale) * c_CF * f_theta_2 * Volume;
    }

    AD::SetPreaccOut(Residual, nVar);
    AD::EndPreacc();

    return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};


/*!
 * \class CIntermittencyVariables
 * \ingroup SourceDiscr
 * \brief Structure with Intermittency common auxiliary functions and constants.
 */
struct CIntermittencyVariables {
  /*--- List of constants ---*/
  const su2double cv1_3 = pow(7.1, 3);
  
  /*--- List of General variables ---*/
  su2double density, laminar_viscosity;
  /*--- List of Fu2013 model ---*/
  su2double Velocity_Mag = 0.0, vel_u = 0.0, vel_v = 0.0, vel_w = 0.0,
            VorticityMag = 0.0, StrainMag = 0.0, dist = 0.0, cordiX = 0.0, cordiY = 0.0,
            Eu = 0.0, norm_Eu = 0.0, norm_k = 0.0, tke = 0.0, omega = 0.0, intermittency = 0.0,
            eddyViscousity = 0.0;
  
  su2double rho_eL = 0.0, U_eL = 0.0, a_eL = 0.0, T_eL = 0.0, Ma_eL = 0.0, He = 0.0;

  bool Liu2022_WallType = false;
            
};


/*!
 * \class CSourceBase_Intermittency
 * \ingroup SourceDiscr
 * \brief Class for integrating the source terms of the intermittency equation.
 * The variables that are subject to change in each variation/correction have their own class. 
 */
template <class FlowIndices>
class CSourceBase_TransIntermittency : public CNumerics {
 protected:

  /*--- Residual and Jacobian ---*/
  su2double Residual, *Jacobian_i;
  su2double Jacobian_Buffer; /*!< \brief Static storage for the Jacobian (which needs to be pointer for return type). */

  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */
  const INTERMITTENCY_ParsedOptions options; /*!< \brief Struct with Intermittency options. */


  bool transition;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceBase_TransIntermittency(unsigned short nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(nDim, 1, config),
        idx(nDim, config->GetnSpecies()),
        options(config->GetINTERMITTENCYParsedOptions()) {
    /*--- Setup the Jacobian pointer, we need to return su2double** but we know
     * the Jacobian is 1x1 so we use this trick to avoid heap allocation. ---*/
    Jacobian_i = &Jacobian_Buffer;
  }


  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
    const auto& density = V_i[idx.Density()];
    const auto& laminar_viscosity = V_i[idx.LaminarViscosity()];

    AD::StartPreacc();
    AD::SetPreaccIn(density, laminar_viscosity, StrainMag_i, Volume, dist_i);
    AD::SetPreaccIn(ScalarVar_i, nVar);    
    AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
    AD::SetPreaccIn(TransVar_i, nVar);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim + idx.Velocity(), nDim);
    AD::SetPreaccIn(Vorticity_i, 3);
    AD::SetPreaccIn(StrainMag_i);
    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);
    AD::SetPreaccIn(V_i[idx.Velocity() + 2]);
    AD::SetPreaccIn(Coord_i, Coord_j);
    AD::SetPreaccIn(V_i[idx.Pressure()], V_i[idx.SoundSpeed()], V_i[idx.Temperature()]);

    /*--- Common auxiliary variables and constants of the model. ---*/
    CIntermittencyVariables var;

    for( int iMarker = 0; iMarker < config->GetnMarker_All() ; iMarker++)
      switch (config->GetMarker_All_KindBC(iMarker))
      {
      case HEAT_FLUX:
        var.Liu2022_WallType = true;
        break;

      case ISOTHERMAL:
        var.Liu2022_WallType = false;
        break;
      
      default:
        break;
      }
    
    var.vel_u = V_i[idx.Velocity()];
    var.vel_v = V_i[1 + idx.Velocity()];
    var.vel_w = (nDim == 3) ? V_i[2 + idx.Velocity()] : 0.0;
    var.Velocity_Mag = sqrt(var.vel_u * var.vel_u + var.vel_v * var.vel_v + var.vel_w * var.vel_w);
    var.density = density;
    var.laminar_viscosity = laminar_viscosity;
    var.dist = dist_i;
    var.tke = ScalarVar_i[0];
    var.omega = ScalarVar_i[1];
    var.intermittency = TransVar_i[0];
    var.cordiX = Coord_i[0];
    var.cordiY = Coord_i[1];
    var.eddyViscousity =  V_i[idx.EddyViscosity()];

    su2double Eu_i = 0.0, Eu_j = 0.0, Eu_k = 0.0, TT = 0.0;

    var.Eu = 0.5*pow( var.Velocity_Mag, 2.0);
    TT = 0.5 * (V_i[idx.Velocity()] * V_i[idx.Velocity()] + V_i[idx.Velocity()+1] * V_i[idx.Velocity()+1]);

    

    Residual = 0.0;
    Jacobian_i[0] = 0.0;

    if (dist_i > 1e-10) {

      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Eu_i += V_i[iDim +idx.Velocity()] * PrimVar_Grad_i[iDim + idx.Velocity()][0];
      Eu_j += V_i[iDim +idx.Velocity()] * PrimVar_Grad_i[iDim + idx.Velocity()][1];
      if(nDim ==3)
      Eu_k += V_i[iDim +idx.Velocity()] * PrimVar_Grad_i[iDim + idx.Velocity()][2];
      }

      var.norm_Eu = pow(Eu_i,2) + pow(Eu_j,2) + pow(Eu_k,2);
      var.norm_Eu = pow(var.norm_Eu, 0.5);

      if(var.Eu == 0 ) {
        var.Eu = 1.0e-12;
        var.norm_Eu = 1.0e-12;
      }

      var.norm_k = pow(ScalarVar_Grad_i[0][0],2) + pow(ScalarVar_Grad_i[0][1],2) ;
      if(nDim == 3) var.norm_k += pow(ScalarVar_Grad_i[0][2],2);
      var.norm_k = pow(var.norm_k, 0.5);

      const su2double VorticityMag = GeometryToolbox::Norm(3, Vorticity_i);
      var.VorticityMag = VorticityMag;
      var.StrainMag = StrainMag_i;

      /*--- Liu2022 Local value calculation ---*/
      const su2double sos = V_i[idx.SoundSpeed()];
      const su2double rho_inf = config->GetDensity_FreeStream();
      const su2double p_inf = config->GetPressure_FreeStream();
      const su2double velU_inf = config->GetVelocity_FreeStream()[0];
      const su2double velV_inf = config->GetVelocity_FreeStream()[1];
      const su2double velW_inf = (nDim ==3) ? config->GetVelocity_FreeStream()[2] : 0.0;
      const su2double velMag_inf = pow(velU_inf*velU_inf + velV_inf * velV_inf + velW_inf *velW_inf,0.5) ;
      const su2double gamma_Spec = config->GetGamma();
      const su2double p = V_i[idx.Pressure()];
      const su2double sos_inf = pow(config->GetTemperature_FreeStream() * config->GetGas_Constant() * gamma_Spec,0.5);
      const su2double temperautre_local = V_i[idx.Temperature()];

      var.rho_eL = pow(rho_inf,gamma_Spec) * p / p_inf;
      var.rho_eL = pow(var.rho_eL,1/gamma_Spec);
      var.U_eL = gamma_Spec / ( gamma_Spec - 1.0) * p_inf / rho_inf + 0.5 * velMag_inf * velMag_inf;
      var.U_eL -= gamma_Spec / ( gamma_Spec - 1.0) * p / var.rho_eL;
      var.U_eL =pow(var.U_eL * 2.0, 0.5) ;
      var.a_eL = sos_inf * sos_inf /(gamma_Spec - 1.0) + velMag_inf * velMag_inf / 2.0;
      var.a_eL -= var.U_eL * var.U_eL /2.0;
      var.a_eL = pow(var.a_eL * (gamma_Spec - 1.0), 0.5);
      var.Ma_eL = var.U_eL / var.a_eL;
      var.T_eL = var.a_eL * var.a_eL / gamma_Spec / config->GetGas_Constant();

      if(nDim == 2) {
        var.He = 0.0;
      }
      else {
        const su2double unitU = V_i[idx.Velocity()]/var.Velocity_Mag, unitV = V_i[idx.Velocity()+1]/var.Velocity_Mag, unitW = V_i[idx.Velocity()+2]/var.Velocity_Mag;
        const su2double vorticity_x = Vorticity_i[0], vorticity_y = Vorticity_i[1], vorticity_z = Vorticity_i[2];
        const su2double UVor_x = unitU * vorticity_x, VVor_y = unitV * vorticity_y, WVor_z = unitW * vorticity_z;

        var.He = pow( UVor_x * UVor_x + VVor_y * VVor_y + WVor_z * WVor_z,0.5);
      }


      /*--- Compute production, destruction, and jacobian ---*/
      su2double Production = 0.0, Destruction = 0.0;
      switch ( options.Intermit_model )
      {
      case INTERMITTENCY_MODEL::FU2013 :
        SourceTerms::Fu2013::get( var, nDim,Production, Destruction, Jacobian_i[0]);
        Residual = ( Production );
        break;
      case INTERMITTENCY_MODEL::WANG2016 :
        SourceTerms::Wang2016::get( var, nDim, Production, Destruction, Jacobian_i[0]);
        Residual = ( Production );
        break;
      case INTERMITTENCY_MODEL::ZHOU2016 :
        SourceTerms::Zhou2016::get( var, nDim, Production, Destruction, Jacobian_i[0]);
        Residual = ( Production );
        break;
      case INTERMITTENCY_MODEL::ZHAO2020 :
        SourceTerms::Zhao2020::get( var, nDim, Production, Destruction, Jacobian_i[0]);
        Residual = ( Production );
        break;
      case INTERMITTENCY_MODEL::MENTER2015 :
        SourceTerms::Menter2015::get( var, nDim, Production, Destruction, Jacobian_i[0]);
        Residual = ( Production );
        break;
      case INTERMITTENCY_MODEL::LIU2022 :
        SourceTerms::Liu2022::get( var, nDim, Production, Destruction, Jacobian_i[0]);
        Residual = ( Production - Destruction);
        break;
      }
      
      
      //SourceTerms::get( var, Production, Destruction, Jacobian_i[0]);

      Residual *= Volume;
      Jacobian_i[0] *= Volume;
    }

    AD::SetPreaccOut(Residual);
    AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
  }

  /*!
 * \brief Intermittency source terms classes: production, destruction and their derivative.
 * \ingroup SourceDiscr
 * \param[in] var: Common SA variables struct.
 * \param[out] production: Production term.
 * \param[out] destruction: Destruction term.
 * \param[out] jacobian: Derivative of the combined source term wrt nue.
 */
struct SourceTerms {

/*! \brief Fu2013 (Transition model). */
struct Fu2013 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {

    const su2double betaStar = 0.09, C_1 = 0.6, C_2 = 0.35, C_3 = 0.005, C_4 = 0.001, C_5 = 5.0;
    const su2double C_6 = 8.0e-5, C_7 = 0.07, C_8 = 1.2, C_mu = 0.09;

    su2double F_onset = 0.0, TempVar1 = 0.0;

    su2double zeta = 0.0, lengthScale_T = 0.0, lengthScale_B = 0.0, zeta_eff = 0.0;
    zeta = pow(var.dist,2) * var.VorticityMag / pow( 2*var.Eu , 0.5);
    lengthScale_T = pow(var.tke, 0.5)/(betaStar * var.omega);
    lengthScale_B = pow(var.tke, 0.5)/(C_mu * var.StrainMag);
    zeta_eff = min(min(zeta, lengthScale_T), C_1 * lengthScale_B);
    //zeta_eff = min(zeta, lengthScale_T);
    F_onset = -exp(-C_8 * zeta_eff * pow( var.tke ,0.5) * var.norm_k / ( var.laminar_viscosity / var.density * var.norm_Eu) );
    TempVar1 = F_onset;
    F_onset += 1.0;

    su2double pg = 0.0;
    if( var.intermittency == 1.0) {
        production = 0.0;
    }
    else {
       su2double templog = -1.0 * log(1.0-var.intermittency) ;
       production = (1.0-var.intermittency) * pow( templog ,0.5) * F_onset * C_6 
                    * ( 1.0 + C_7 * pow( var.tke / 2.0 / var.Eu ,0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu ;
        }
    if( var.intermittency != 1.0 || var.intermittency != 0.0) {
      su2double templog = -1.0 * log(1.0-var.intermittency) ;
      jacobian = F_onset * C_6 * ( 1.0 + C_7 * pow( var.tke / 2.0 / var.Eu ,0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu *
               (log(1.0 - var.intermittency) +0.5)/var.density/pow( templog, 0.5);
    }
    else {
      jacobian = 0.0;
    }
    
  }

  static void ComputeDestruction( const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
    destruction = 0.0;
  }

};


/*! \brief Wang2016 (Transition model). */
struct Wang2016 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {
    const su2double C_1 = 0.7, C_2 = 0.35, C_3 = 0.005, C_4 = 0.001, C_5 = 5.0;
    const su2double C_6 = 8.0e-5, C_7 = 0.07, C_8 = 1.2, C_mu = 0.09;

    su2double F_onset = 0.0, TempVar1 = 0.0;

    su2double zeta = 0.0, lengthScale_T = 0.0, zeta_eff = 0.0;
    zeta = pow(var.dist,2) * var.VorticityMag / pow( 2*var.Eu , 0.5);
    lengthScale_T = pow(var.tke, 0.5)/(var.omega);
    zeta_eff = min(zeta, lengthScale_T);
    
    F_onset = -exp(-C_8 * zeta_eff * pow( var.tke ,0.5) * var.norm_k / ( var.laminar_viscosity / var.density * var.norm_Eu) );
    TempVar1 = F_onset;
    F_onset += 1.0;

    if( var.intermittency == 1.0) {
        production = 0.0;
    }
    else {
       su2double templog = -1.0 * log(1.0-var.intermittency) ;
       production = (1.0-var.intermittency) * pow( templog ,0.5) * F_onset * C_6 
                    * ( 1.0 + C_7 * pow( var.tke / 2.0 / var.Eu, 0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu ;
    }
    if( var.intermittency != 1.0 || var.intermittency != 0.0) {
      su2double templog = -1.0 * log(1.0-var.intermittency) ;
      jacobian = F_onset * C_6 * ( 1.0 + C_7 * pow( var.tke / 2.0 / var.Eu ,0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu *
               (log(1.0 - var.intermittency) +0.5)/var.density/pow( templog, 0.5);
    }
    else {
      jacobian = 0.0;
    }
  }

  static void ComputeDestruction(const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
    destruction = 0.0;
  }

};


/*! \brief Zhou2016 (Transition model). */
struct Zhou2016 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {

    const su2double betaStar = 0.09, C_1 = 0.7, C_2 = 0.35, C_3 = 0.005, C_4 = 8.0e-5, C_5 = 0.07;
    const su2double C_6 = 1.2, C_7 = 0.05, C_8 = 1.0, C_mu = 0.09;

    su2double F_onset = 0.0, TempVar1 = 0.0, Tu = 0.0;
    Tu = 100.0 * sqrt(2.0 *var.tke / 3.0) / pow( 2*var.Eu , 0.5);

    su2double zeta = 0.0, lengthScale_T = 0.0, lengthScale_B = 0.0, zeta_eff = 0.0;
    zeta = pow(var.dist,2) * var.VorticityMag / pow( 2*var.Eu , 0.5);
    lengthScale_T = pow(var.tke, 0.5)/( var.omega);
    lengthScale_B = pow(var.tke, 0.5)/(C_mu * var.StrainMag);
    zeta_eff = min(min(zeta, lengthScale_T), lengthScale_B);
    //zeta_eff = min(zeta, lengthScale_T);
    F_onset = -exp(-C_8 * var.eddyViscousity/ var.laminar_viscosity ) ;
    TempVar1 = F_onset;
    F_onset += 1.0;

    su2double pg = 0.0;
    if( var.intermittency == 1.0) {
        production = 0.0;
    }
    else {
       su2double templog = -1.0 * log(1.0-var.intermittency) ;
       production = C_7 * var.density * pow( Tu ,7/4*0.5) * F_onset * 2.0 * var.StrainMag * 
                    var.density * var.dist * var.dist * var.StrainMag / var.laminar_viscosity *
                    pow( templog ,0.5) * (1.0-var.intermittency) ;
    }
    if( var.intermittency != 1.0 || var.intermittency != 0.0) {
      su2double templog = -1.0 * log(1.0-var.intermittency) ;
      jacobian = C_7 * var.density * pow( Tu ,7/4*0.5) * F_onset * 2.0 * var.StrainMag * 
                    var.density * var.dist * var.dist * var.StrainMag / var.laminar_viscosity *
               (log(1.0 - var.intermittency) +0.5)/var.density/pow( templog, 0.5);
    }
    else {
      jacobian = 0.0;
    }
  }

  static void ComputeDestruction( const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
    destruction = 0.0;
  }

};

/*! \brief Zhao2020 (Transition model). */
struct Zhao2020 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {
    const su2double C_1 = 700, C_2 = 0.5, C_3 = 0.5, C_4 = 8.0e-5, C_5 = 0.07, C_6 = 1.2, C_mu = 0.09;

    su2double F_onset = 0.0, TempVar1 = 0.0;

    su2double zeta = 0.0, lengthScale_T = 0.0, zeta_eff = 0.0;
    zeta = pow(var.dist,2) * var.VorticityMag / pow( 2*var.Eu , 0.5);
    lengthScale_T = pow(var.tke, 0.5)/(var.omega);
    zeta_eff = min(zeta, C_1 * lengthScale_T);
    
    F_onset = -exp(-C_6 * zeta_eff * pow( var.tke ,0.5) * var.norm_k / ( var.laminar_viscosity / var.density * var.norm_Eu) );
    TempVar1 = F_onset;
    F_onset += 1.0;

    if( var.intermittency == 1.0) {
        production = 0.0;
    }
    else {
       su2double templog = -1.0 * log(1.0-var.intermittency) ;
       production = (1.0-var.intermittency) * pow( templog ,0.5) * F_onset * C_4
                    * ( 1.0 + C_5 * pow( var.tke / 2.0 / var.Eu, 0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu ;
    }
    if( var.intermittency != 1.0 || var.intermittency != 0.0) {
      su2double templog = -1.0 * log(1.0-var.intermittency) ;
      jacobian = F_onset * C_4 * ( 1.0 + C_5 * pow( var.tke / 2.0 / var.Eu, 0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu *
               (log(1.0 - var.intermittency) +0.5)/var.density/pow( templog, 0.5);
    }
    else {
      jacobian = 0.0;
    }
  }

  static void ComputeDestruction(const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
    destruction = 0.0;
  }

};

/*! \brief Menter2015 (Transition model). */
struct Menter2015 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {
    const su2double C_1 = 300, C_2 = 0.5, C_3 = 0.5, C_4 = 8.0e-5, C_5 = 0.07, C_6 = 1.2, C_mu = 0.09;

    su2double F_onset = 0.0, TempVar1 = 0.0;

    su2double zeta = 0.0, lengthScale_T = 0.0, zeta_eff = 0.0;
    zeta = pow(var.dist,2) * var.VorticityMag / pow( 2*var.Eu , 0.5);
    lengthScale_T = pow(var.tke, 0.5)/(var.omega);
    zeta_eff = min(zeta, C_1 * lengthScale_T);
    
    F_onset = -exp(-C_6 * zeta_eff * pow( var.tke ,0.5) * var.norm_k / ( var.laminar_viscosity / var.density * var.norm_Eu) );
    TempVar1 = F_onset;
    F_onset += 1.0;

    if( var.intermittency == 1.0) {
        production = 0.0;
    }
    else {
       su2double templog = -1.0 * log(1.0-var.intermittency) ;
       production = (1.0-var.intermittency) * pow( templog ,0.5) * F_onset * C_4
                    * ( 1.0 + C_5 * pow( var.tke / 2.0 / var.Eu, 0.5) ) * var.dist * var.density / var.laminar_viscosity * var.norm_Eu ;
    }
    jacobian = 0.0;
  }

  static void ComputeDestruction(const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
    destruction = 0.0;
  }

};

/*! \brief Liu2022 (Transition model). */
struct Liu2022 {
  static void get( const CIntermittencyVariables& var, const su2double& nDim, su2double& production, su2double& destruction,
                   su2double& jacobian) {
    ComputeProduction(var, nDim, production, jacobian);
    ComputeDestruction(var, nDim, destruction, jacobian);
  }

  static void ComputeProduction(const CIntermittencyVariables& var, const su2double& nDim, su2double& production,
                                su2double& jacobian) {

    su2double fMaeLTeL = 0.0, f_a = 0.0, f_b = 0.0, f_c = 0.0;
    su2double TuL = 0.0, RetMae_c = 0.0, Re_tc = 0.0;
    su2double delH_cf = 0.0, H_cf = 0.0;

    su2double Fonset_s = 0.0, Fonset_cf = 0.0, F_onset1 = 0.0,
              F_onset2 = 1.0, F_onset3 = 1.0, F_nose = 1.0;
    const su2double C_1 = 100.0;
  
    if(var.Liu2022_WallType){
      //heat flux
      f_a = -0.003755 * exp(-0.01786 * var.T_eL) - 0.004692 * exp(1.26e-6 * var.T_eL);
      f_b = 0.1141 * exp(-0.02721 * var.T_eL) + 0.09247 * exp(2.527e-6 * var.T_eL);
      f_c = 0.1066 * exp(-1.011e-4 * var.T_eL) -0.1889 * exp(-0.03724 * var.T_eL);
      fMaeLTeL = f_a * pow(var.Ma_eL, 3) + f_b * pow(var.Ma_eL, 2) + f_c * var.Ma_eL + 2.2;
    } 
    else {
      //isothermal
      f_a = -0.1871 * pow(var.Ma_eL, 3) + 2.9  * pow(var.Ma_eL, 2) + 2.958 * var.Ma_eL + 99.41;
      f_b = 4.715e-4 * pow(var.Ma_eL, 3) - 0.006389 * pow(var.Ma_eL, 2) - 0.02174 * var.Ma_eL - 0.7565;
      f_c = -0.001551 * pow(var.Ma_eL, 3) + 0.03085 * pow(var.Ma_eL, 2) + 0.01819  * var.Ma_eL + 0.8891;
      fMaeLTeL = f_a * pow(var.T_eL,f_b) + f_c;
    }

    TuL = min(100.0 * sqrt(2.0 * var.tke /3.0) / var.omega / var.dist , 100.0);
    RetMae_c = 1034.0 * exp(-97.56 * TuL * F_nose) + 440.0 * exp(-1.96*TuL * F_nose);
    Re_tc = RetMae_c * var.Ma_eL;
    Re_tc = min( max(Re_tc, 100.0), 1500.0);

    H_cf = var.dist * var.He / var.Velocity_Mag;
    delH_cf = H_cf * (1.0 + min(var.eddyViscousity/var.laminar_viscosity, 0.3));
    const su2double Re_v = var.density * var.dist * var.dist * var.StrainMag / var.laminar_viscosity;

    Fonset_s = Re_v / fMaeLTeL / Re_tc;
    Fonset_cf = delH_cf * Re_v / fMaeLTeL / 46.0;

    const su2double R_t = var.density * var.tke / var.laminar_viscosity / var.omega;
  
    F_onset1 = max(Fonset_s, Fonset_cf);
    F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
    F_onset3 = max(1.0 - pow( R_t / 3.5, 3.0), 0.0);
    if( var.cordiX > 0.05 && var.cordiY < 0.00001){
      su2double testmpppp = 0.0;
    }

    const su2double F_onset = max(F_onset2 - F_onset3, 0.0);
    production = C_1 * var.density * var.StrainMag * var.intermittency * (1.0 - var.intermittency) * F_onset;
    jacobian += C_1 * var.StrainMag * (1.0 - 2.0 * var.intermittency) * F_onset;

  }

  static void ComputeDestruction(const CIntermittencyVariables& var, const su2double& nDim, su2double& destruction,
                                 su2double& jacobian) {
                                  
    const su2double R_t = var.density * var.tke / var.laminar_viscosity / var.omega;
    const su2double f_turb = exp(-pow(R_t / 2, 2));
    const su2double C_2 = 0.06, C_3 = 50.0;
    destruction = C_2 * var.density * var.VorticityMag * var.intermittency * f_turb * (C_3 * var.intermittency -1.0);
    jacobian -= C_2 * var.VorticityMag * f_turb * (C_3 * 2.0 * var.intermittency -1.0);

  }

};





};

};