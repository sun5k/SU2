/*!
 * \file CTransIntermittencySolver.cpp
 * \brief Main subroutines for intermittency based Transition model solver.
 * \author A. Aranake, S. Kang.
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

#include "../../include/solvers/CTransIntermittencySolver.hpp"
#include "../../include/variables/CTransIntermittencyVariable.hpp"
#include "../../include/variables/CFlowVariable.hpp"
#include "../../include/variables/CTurbSAVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

/*---  This is the implementation of the Langtry-Menter transition model.
       The main reference for this model is:Langtry, Menter, AIAA J. 47(12) 2009
       DOI: https://doi.org/10.2514/1.42362 ---*/

// Note: TransIntermittency seems to use rho*gamma as Solution variables, thus Conservative=true

CTransIntermittencySolver::CTransIntermittencySolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
    : CTurbSolver(geometry, config, true) {
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;

  bool multizone = config->GetMultizone_Problem();
  bool MenterKOG = false;

  if( config -> GetINTERMITTENCYParsedOptions().Intermit_model == INTERMITTENCY_MODEL::LIU2022 ) MenterKOG = true;
  if( config -> GetINTERMITTENCYParsedOptions().Intermit_model == INTERMITTENCY_MODEL::MENTER2015 ) MenterKOG = true;

  /*--- Dimension of the problem --> 1 Transport equations (intermittency) ---*/
  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();  

  /*--- Single grid simulation ---*/

  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/

    Residual_RMS.resize(nVar,0.0);
    Residual_Max.resize(nVar,0.0);
    Point_Max.resize(nVar,0);
    Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);

    /*--- Initialization of the structure of the whole Jacobian ---*/

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Intermittency based transition model)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
    System.SetxIsZero(true);

    if (ReducerStrategy)
      EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

    /*--- Initialize the BGS residuals in multizone problems. ---*/
    if (multizone){
      Residual_BGS.resize(nVar,0.0);
      Residual_Max_BGS.resize(nVar,0.0);
      Point_Max_BGS.resize(nVar,0);
      Point_Max_Coord_BGS.resize(nVar,nDim) = su2double(0.0);
    }

  }

  /*--- Initialize lower and upper limits---*/
  lowerlimit[0] = 1.0e-11;
  upperlimit[0] = 1.0;

  /*--- Far-field flow state quantities and initialization. ---*/  

  su2double Intermittency_Inf  = 0.02;

  if(MenterKOG) {
    Intermittency_Inf = 1.0;
    lowerlimit[0] = 0.02;
    upperlimit[0] = 1.0;
  }

  Solution_Inf[0] = Intermittency_Inf;

  /*--- Initialize the solution to the far-field state everywhere. ---*/
  nodes = new CTransIntermittencyVariable(Intermittency_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  SlidingState.resize(nMarker);
  SlidingStateNodes.resize(nMarker);

  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE) {
      SlidingState[iMarker].resize(nVertex[iMarker], nPrimVar+1) = nullptr;
      SlidingStateNodes[iMarker].resize(nVertex[iMarker],0);
    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
    due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars.resize(nMarker);
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_TurbVars[iMarker].resize(nVertex[iMarker],nVar);
    for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; ++iVertex) {
      Inlet_TurbVars[iMarker](iVertex,0) = Intermittency_Inf;
    }
  }

   const su2double CFL = config->GetCFL(MGLevel)*config->GetCFLRedCoeff_Turb();
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    nodes->SetLocalCFL(iPoint, CFL);
  }
  Min_CFL_Local = CFL;
  Max_CFL_Local = CFL;
  Avg_CFL_Local = CFL;

  /*--- Add the solver name. ---*/
  SolverName = "Intermittency model";

}

void CTransIntermittencySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
         unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)

  /*--- Upwind second order reconstruction and gradients ---*/
  CommonPreprocessing(geometry, config, Output);
}

void CTransIntermittencySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  /*--- Compute LM model gradients. ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config);
  }
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config);
  }

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  auto* turbNodes = su2staticcast_p<CTurbVariable*>(solver_container[TURB_SOL]->GetNodes());


  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    const su2double betaStar = 0.09, C_1 = 0.6, C_2 = 0.35, C_3 = 0.005, C_4 = 0.001, C_5 = 5.0;
    const su2double C_6 = 8.0e-5, C_7 = 0.07, C_8 = 1.2, C_mu = 0.09;

    const su2double Intermittency = nodes->GetSolution(iPoint,0);
    const su2double rho = flowNodes->GetDensity(iPoint);
    const su2double mu = flowNodes->GetLaminarViscosity(iPoint);
    const su2double p = flowNodes->GetPressure(iPoint);
    const su2double sos = flowNodes->GetSoundSpeed(iPoint);    
    const su2double rho_inf = config->GetDensity_FreeStream();
    const su2double p_inf = config->GetPressure_FreeStream();
    const su2double velU_inf = config->GetVelocity_FreeStream()[0];
    const su2double velV_inf = config->GetVelocity_FreeStream()[1];
    const su2double velW_inf = (nDim ==3) ? config->GetVelocity_FreeStream()[2] : 0.0;
    const su2double velMag_inf = pow(velU_inf*velU_inf + velV_inf * velV_inf + velW_inf *velW_inf,0.5) ;
    const su2double gamma_Spec = config->GetGamma();
    const su2double dist = geometry->nodes->GetWall_Distance(iPoint);
    const su2double VorticityMag = max(GeometryToolbox::Norm(3, flowNodes->GetVorticity(iPoint)), 1e-12);
    const su2double StrainMag = max(nodes->GetStrainMag(iPoint), 1e-12);
    const su2double sos_inf = pow(config->GetTemperature_FreeStream() * config->GetGas_Constant() * gamma_Spec,0.5);
    
    su2double tke = 0.0, omega = 0.0, norm_k = 0.0;
    tke = turbNodes->GetSolution(iPoint,0);
    omega = turbNodes->GetSolution(iPoint,1);
    
    norm_k = pow(turbNodes->GetGradient(iPoint, 0, 0), 2.0) + pow(turbNodes->GetGradient(iPoint, 0, 1), 2.0);
    if(nDim == 3) norm_k += pow(turbNodes->GetGradient(iPoint, 0, 2), 2.0);
    norm_k = pow(norm_k,0.5);
    


    /*--- Compute the Time and Length Scale ---*/
    su2double Eu = 0.0, Eu_i = 0.0, Eu_j = 0.0, Eu_k = 0.0, norm_Eu = 0.0;
    const su2double vel_u = flowNodes->GetVelocity(iPoint, 0);
    const su2double vel_v = flowNodes->GetVelocity(iPoint, 1);
    const su2double vel_w = (nDim ==3) ? flowNodes->GetVelocity(iPoint, 2) : 0.0;
    const su2double Velocity_Mag = sqrt(vel_u * vel_u + vel_v * vel_v + vel_w * vel_w);


    Eu = 0.5*pow(Velocity_Mag,2);
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        Eu_i += flowNodes->GetVelocity(iPoint, iDim) * flowNodes->GetGradient_Primitive(iPoint, iDim, 0);
        Eu_j += flowNodes->GetVelocity(iPoint, iDim) * flowNodes->GetGradient_Primitive(iPoint, iDim, 1);
        if(nDim ==3)
        Eu_k += flowNodes->GetVelocity(iPoint, iDim) * flowNodes->GetGradient_Primitive(iPoint, iDim, 2);
      }
    norm_Eu = pow(Eu_i,2) + pow(Eu_j,2) + pow(Eu_k,2);
    norm_Eu = pow(norm_Eu, 0.5);
    if(Eu == 0 ) {
      Eu = 1.0e-12;
      norm_Eu = 1.0e-12;
    }

    su2double zeta = 0.0, lengthScale_T = 0.0, lengthScale_B = 0.0, zeta_eff = 0.0;
    zeta = pow(dist,2) * VorticityMag / pow( 2*Eu , 0.5);


    su2double M_rel = 0.0, Cr = 0.0, rho_e = 0.0, velMag_e = 0.0;
    rho_e = pow(rho_inf,gamma_Spec)/p_inf*p;
    rho_e = pow(rho_e,1/gamma_Spec);
    velMag_e = gamma_Spec / ( gamma_Spec - 1.0) * p_inf / rho_inf + 0.5 * velMag_inf * velMag_inf;
    velMag_e -= gamma_Spec / ( gamma_Spec - 1.0) * p / rho_e;
    velMag_e =pow(velMag_e * 2.0, 0.5) ;
    Cr = 0.94 * velMag_e;    
    M_rel = (Velocity_Mag - Cr)/sos;

    su2double M_rel2 = 0.0, a_eL = 0.0;

    a_eL = sos_inf * sos_inf /(gamma_Spec - 1.0) + velMag_inf * velMag_inf / 2.0;
    a_eL -= velMag_e * velMag_e /2.0;
    a_eL = pow(a_eL * (gamma_Spec - 1.0), 0.5);
    M_rel = (Velocity_Mag - Cr)/a_eL;

    /*--- Tau_nt1 and Tau_nt2 ---*/
    su2double Tau_nt1_extra = 0.0,  sgn = 0.0 , Tau_extra = 0.0;
    Tau_nt1_extra = C_2 * rho / pow( pow( 2.0 * Eu, 0.5) * mu ,0.5);
    sgn = abs(M_rel - 1.0) / (M_rel -1.0);

    

    su2double F_onset = 0.0, TempVar1 = 0.0;
    F_onset = -exp(-1.2 * zeta_eff * pow( tke ,0.5) * norm_k / ( mu / rho * norm_Eu) );
    TempVar1 = F_onset;
    F_onset += 1.0;
    su2double test3 = 0.0, test4 = 0.0, test5 = 0.0, test6 = 0.0, test7 = 0.0; 
    su2double pg = 0.0, templog = 0.0;
    if( Intermittency == 1.0) {
        pg = 0.0;
    }
    else {
      templog = -1.0 * log(1.0-Intermittency) ;
      pg = (1.0- Intermittency) * pow( templog ,0.5) * F_onset * C_6 
                    * ( 1 + C_7 * pow( tke / 2.0 / Eu ,0.5) ) * dist * rho / mu * norm_Eu ;
      test4 = pow( templog ,0.5);
      test5 = pow( tke / 2.0 / Eu ,0.5);
      test6 = dist * rho / mu * norm_Eu;
      test7 = ( 1 + C_7 * pow( tke / 2.0 / Eu, 0.5) );

    }

    su2double He = 0;


    if(nDim == 2) {
        He = 0.0;
      }
      else {
        const su2double unitU = vel_u/Velocity_Mag, unitV = vel_v/Velocity_Mag, unitW = vel_w/Velocity_Mag;
        const su2double vorticity_x = flowNodes->GetVorticity(iPoint)[0], vorticity_y = flowNodes->GetVorticity(iPoint)[1], vorticity_z = flowNodes->GetVorticity(iPoint)[2];
        const su2double UVor_x = unitU * vorticity_x, VVor_y = unitV * vorticity_y, WVor_z = unitW * vorticity_z;

        He = pow( UVor_x * UVor_x + VVor_y * VVor_y + WVor_z * WVor_z,0.5);
      }

      su2double Ma_eL = velMag_e / a_eL;
      su2double T_eL = a_eL * a_eL / gamma_Spec/ config->GetGas_Constant();
      su2double fMaeLTeL = 0.0, f_a = 0.0, f_b = 0.0, f_c = 0.0;
      f_a = -0.1871 * pow(Ma_eL, 3) + 2.9  * pow(Ma_eL, 2) + 2.958 * Ma_eL + 99.41;
      f_b = 4.715e-4 * pow(Ma_eL, 3) - 0.006389 * pow(Ma_eL, 2) - 0.02174 * Ma_eL - 0.7565;
      f_c = -0.001551 * pow(Ma_eL, 3) + 0.03085 * pow(Ma_eL, 2) + 0.01819  * Ma_eL + 0.8891;
      fMaeLTeL = f_a * pow(T_eL,f_b) + f_c;




    nodes -> SetIntermittency(iPoint, Intermittency);

    switch ( config->GetINTERMITTENCYParsedOptions().Intermit_model )
    {
    case INTERMITTENCY_MODEL::NONE :
      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_nt1_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      break;
    case INTERMITTENCY_MODEL::FU2013 :
      lengthScale_T = pow(tke, 0.5)/(betaStar * omega);
      lengthScale_B = pow(tke, 0.5)/(C_mu * StrainMag);
      zeta_eff = min(min(zeta, lengthScale_T), C_1 * lengthScale_B);

      F_onset = -exp(-C_8 * zeta_eff * pow( tke ,0.5) * norm_k / ( mu / rho * norm_Eu) );
      TempVar1 = F_onset;
      F_onset += 1.0;

      if( Intermittency == 1.0) {
        pg = 0.0;
    }
    else {
      templog = -1.0 * log(1.0- Intermittency) ;
      pg = (1.0- Intermittency) * pow( templog ,0.5) * F_onset * C_6 
                    * ( 1 + C_7 * pow( tke / 2.0 / Eu ,0.5) ) * dist * rho / mu * norm_Eu ;
        
      test4 = pow( templog ,0.5);
      test5 = pow( tke / 2.0 / Eu ,0.5);
      test6 = dist * rho / mu * norm_Eu, 
      test7 = ( 1.0 + C_7 * pow( tke / 2.0 / Eu ,0.5) );
        
    }
      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_nt1_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, zeta, zeta_eff, F_onset, pg, test4, test7);
      break;


    case INTERMITTENCY_MODEL::WANG2016 :
      lengthScale_T = pow(tke, 0.5)/(omega);
      zeta_eff = min(zeta, lengthScale_T);      

      F_onset = -exp(-1.2 * zeta_eff * pow( tke ,0.5) * norm_k / ( mu / rho * norm_Eu) );
      TempVar1 = F_onset;
      F_onset += 1.0;

      if( Intermittency == 1.0) {
        pg = 0.0;
      }
      else {
        templog = -1.0 * log(1.0-Intermittency) ;
        pg = (1.0- Intermittency) * pow( templog ,0.5) * F_onset * C_6 
                    * ( 1 + C_7 * pow( tke / 2.0 / Eu ,0.5) ) * dist * rho / mu * norm_Eu ;
        test4 = pow( templog ,0.5);
        test5 = pow( tke / 2.0 / Eu ,0.5);
        test6 = dist * rho / mu * norm_Eu;
        test7 = ( 1.0 + C_7 * pow( tke / 2.0 / Eu, 0.5) );
      }

      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_nt1_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, zeta, zeta_eff, F_onset, pg, test4, test7);
      break;
    case INTERMITTENCY_MODEL::ZHOU2016 :
      if(M_rel > 1.0) {
        Tau_extra = 0.005 * 2.0 /Cr;
      } else {
        Tau_extra = C_2 / pow( pow( 2.0 * Eu, 0.5) * mu / rho ,0.5);
      }
      nodes -> SetIntermittency_Zhou_Func(iPoint, zeta, M_rel, Tau_extra);
      nodes -> SetIntermittency_Wonder_Func(iPoint, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      break;

    case INTERMITTENCY_MODEL::ZHAO2020 :
      lengthScale_T = pow(tke, 0.5)/(omega);
      zeta_eff = min(zeta, 700.0 * lengthScale_T);      

      F_onset = -exp(-1.2 * zeta_eff * pow( tke ,0.5) * norm_k / ( mu / rho * norm_Eu) );
      TempVar1 = F_onset;
      F_onset += 1.0;

      if( Intermittency == 1.0) {
        pg = 0.0;
      }
      else {
        templog = -1.0 * log(1.0-Intermittency) ;
        pg = (1.0- Intermittency) * pow( templog ,0.5) * F_onset * 8.0e-5 
                    * ( 1.0 + 0.07 * pow( tke / 2.0 / Eu ,0.5) ) * dist * rho / mu * norm_Eu ;
        test4 = pow( templog ,0.5);
        test5 = pow( tke / 2.0 / Eu ,0.5);
        test6 = dist * rho / mu * norm_Eu;
        test7 = ( 1.0 + 0.07 * pow( tke / 2.0 / Eu, 0.5) );
      }

      Tau_extra = 0.5 / pow( pow( 2.0 * Eu, 0.5) * mu/rho ,0.5);
      nodes -> SetIntermittency_Zhou_Func(iPoint, zeta, M_rel, Tau_extra);
      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, zeta, zeta_eff, F_onset, pg, test4, test7);
      break;

    case INTERMITTENCY_MODEL::MENTER2015 :
      lengthScale_T = pow(tke, 0.5)/(omega);
      zeta_eff = min(zeta, 300.0 * lengthScale_T);      

      F_onset = -exp(-1.2 * zeta_eff * pow( tke ,0.5) * norm_k / ( mu / rho * norm_Eu) );
      TempVar1 = F_onset;
      F_onset += 1.0;

      if( Intermittency == 1.0) {
        pg = 0.0;
      }
      else {
        templog = -1.0 * log(1.0-Intermittency) ;
        pg = (1.0- Intermittency) * pow( templog ,0.5) * F_onset * 8.0e-5 
                    * ( 1.0 + 0.07 * pow( tke / 2.0 / Eu ,0.5) ) * dist * rho / mu * norm_Eu ;
        test4 = pow( templog ,0.5);
        test5 = pow( tke / 2.0 / Eu ,0.5);
        test6 = dist * rho / mu * norm_Eu;
        test7 = ( 1.0 + 0.07 * pow( tke / 2.0 / Eu, 0.5) );
      }

      Tau_extra = 0.5 / pow( pow( 2.0 * Eu, 0.5) * mu/rho ,0.5);
      nodes -> SetIntermittency_Zhou_Func(iPoint, zeta, M_rel, Tau_extra);
      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, zeta, zeta_eff, F_onset, pg, test4, test7);
      break;

    case INTERMITTENCY_MODEL::LIU2022 :
      nodes -> SetIntermittency_Wonder_Func(iPoint, M_rel, He, 0.0, 0.0, 0.0, 0.0 );
      break;


    default:
      nodes -> SetIntermittency_Fu_Func(iPoint, zeta, sgn, Tau_nt1_extra, Cr);
      nodes -> SetIntermittency_Wonder_Func(iPoint, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
      break;
    }

  }

}


void CTransIntermittencySolver::Viscous_Residual(unsigned long iEdge, CGeometry* geometry, CSolver** solver_container,
                                     CNumerics* numerics, CConfig* config) {

  /*--- Define an object to set solver specific numerics contribution. ---*/

  auto SolverSpecificNumerics = [&](unsigned long iPoint, unsigned long jPoint) {};

  /*--- Now instantiate the generic implementation with the functor above. ---*/

  Viscous_Residual_impl(SolverSpecificNumerics, iEdge, geometry, solver_container, numerics, config);
}


void CTransIntermittencySolver::Source_Residual(CGeometry *geometry, CSolver **solver_container,
                                     CNumerics **numerics_container, CConfig *config, unsigned short iMesh) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  auto* flowNodes = su2staticcast_p<CFlowVariable*>(solver_container[FLOW_SOL]->GetNodes());
  //auto* turbNodes = su2staticcast_p<CFlowVariable*>(solver_container[TURB_SOL]->GetNodes());
  CVariable* turbNodes = solver_container[TURB_SOL]->GetNodes();

  /*--- Pick one numerics object per thread. ---*/
  auto* numerics = numerics_container[SOURCE_FIRST_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Loop over all points. ---*/

  AD::StartNoSharedReading();

  SU2_OMP_FOR_DYN(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {



    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(flowNodes->GetPrimitive(iPoint), nullptr);

    /*--- Gradient of the primitive and conservative variables ---*/

    numerics->SetPrimVarGradient(flowNodes->GetGradient_Primitive(iPoint), nullptr);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    /*--- ScalarVar & ScalarVarGradient : Turbulence model solution(k&w) ---*/

    numerics->SetScalarVar(turbNodes->GetSolution(iPoint), nullptr);
    numerics->SetScalarVarGradient(turbNodes->GetGradient(iPoint), nullptr);

    /*--- Transition variables w/o reconstruction, and its gradient ---*/

    numerics->SetTransVar(nodes->GetSolution(iPoint), nullptr);
    numerics->SetTransVarGradient(nodes->GetGradient(iPoint), nullptr);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->nodes->GetVolume(iPoint));

    /*--- Set distance to the surface ---*/

    numerics->SetDistance(geometry->nodes->GetWall_Distance(iPoint), 0.0);

    /*--- Set vorticity and strain rate magnitude ---*/

    numerics->SetVorticity(flowNodes->GetVorticity(iPoint), nullptr);

    numerics->SetStrainMag(flowNodes->GetStrainMag(iPoint), 0.0);

    /*--- Set coordinate (for debugging) ---*/
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint), nullptr);    

    /*--- Compute the source term ---*/

    auto residual = numerics->ComputeResidual(config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, residual);
    if (implicit) Jacobian.SubtractBlock2Diag(iPoint, residual.jacobian_i);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

}

void CTransIntermittencySolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                       CConfig *config, unsigned short iMesh) {
}

void CTransIntermittencySolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Allocate the value at the infinity ---*/

      auto V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at the wall, and at infinity ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Solution_Inf);

      /*--- Set Normal (it is necessary to change the sign) ---*/
      /*--- It's mean wall normal zero flux. */

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute residuals and Jacobians ---*/

      auto residual = conv_numerics->ComputeResidual(config);

      /*--- Add residuals and Jacobians ---*/

      LinSysRes.AddBlock(iPoint, residual);
      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransIntermittencySolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                        CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);

}

void CTransIntermittencySolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Loop over all the vertices on this boundary marker ---*/

  SU2_OMP_FOR_STAT(OMP_MIN_SIZE)
  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      su2double Normal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Normal[iDim] = -geometry->vertex[val_marker][iVertex]->GetNormal(iDim);
      conv_numerics->SetNormal(Normal);

      /*--- Allocate the value at the inlet ---*/

      auto V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/

      auto V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);

      /*--- Set various quantities in the solver class ---*/

      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Non-dimensionalize Inlet_TurbVars if Inlet-Files are used. ---*/
      su2double Inlet_Vars[MAXNVAR];
      Inlet_Vars[0] = Inlet_TurbVars[val_marker][iVertex][0];
      if (config->GetInlet_Profile_From_File()) {
        Inlet_Vars[0] /= pow(config->GetVelocity_Ref(), 2);
      }

      /*--- Set the LM variable states. ---*/
      /*--- Load the inlet transition LM model variables (uniform by default). ---*/

      conv_numerics->SetScalarVar(nodes->GetSolution(iPoint), Inlet_Vars);

      /*--- Set various other quantities in the solver class ---*/

      if (dynamic_grid)
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint),
                                  geometry->nodes->GetGridVel(iPoint));

      /*--- Compute the residual using an upwind scheme ---*/

      auto residual = conv_numerics->ComputeResidual(config);
      LinSysRes.AddBlock(iPoint, residual);

      /*--- Jacobian contribution for implicit integration ---*/

      if (implicit) Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
    }
  }
  END_SU2_OMP_FOR

}

void CTransIntermittencySolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_Far_Field(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
}

void CTransIntermittencySolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                                  bool val_update_geo) {

  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar() + solver[MESH_0][TURB_SOL] ->GetnVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    const bool energy = config->GetEnergy_Equation();
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();

    if (incompressible && ((!energy) && (!weakly_coupled_heat))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        const auto index = counter * Restart_Vars[1] + skipVars;
        for (auto iVar = 0u; iVar < nVar; iVar++) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);
        nodes ->SetIntermittency_Wonder_Func(iPoint_Local, Restart_Data[index + 2], Restart_Data[index + 3]
              , Restart_Data[index + 4], Restart_Data[index + 5], Restart_Data[index + 6], Restart_Data[index + 7]);

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file!\n") +
                         string("This can be caused by empty lines at the end of the file."),
                     CURRENT_FUNCTION);
    }

  }  // end SU2_OMP_MASTER, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][TRANS_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][TRANS_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  /*--- For turbulent+species simulations the solver Pre-/Postprocessing is done by the species solver. ---*/
  if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                            RUNTIME_FLOW_SYS, false);
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
    solver[MESH_0][TRANS_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {

    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][TRANS_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][TRANS_SOL]->GetNodes()->GetSolution());
    solver[iMesh][TRANS_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][TRANS_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

    if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                            false);
      solver[iMesh][TRANS_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
    }
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Delete the class memory that is used to load the restart. ---*/

    delete[] Restart_Vars;
    Restart_Vars = nullptr;
    delete[] Restart_Data;
    Restart_Data = nullptr;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

}