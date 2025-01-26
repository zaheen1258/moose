//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMixtureOfExpertsElastoViscoplasticStress.h"
#include <torch/script.h> // One-stop header.
#include <torch/torch.h>
#include <iostream>
#include <memory>

registerMooseObject("SolidMechanicsApp", ADMixtureOfExpertsElastoViscoplasticStress);

InputParameters
ADMixtureOfExpertsElastoViscoplasticStress::validParams()
{
  InputParameters params = ADComputeStressBase::validParams();
  params.addClassDescription("Computes stress after subtracting inelastic strain increment calculated with" 
                                    "Mixture of Expert ML based elasto-viscoplastic constitutive relation");
  params.addRequiredParam<MaterialPropertyName>("mu", "Shear modulus");
  params.addRequiredParam<MaterialPropertyName>("temp", "Temperature");
  params.addRequiredParam<Real>("initial_rhoc", "Initial Dislocation density at cell");
  params.addRequiredParam<Real>("initial_rhow", "Initial Dislocation density at cell wall");
  params.addRequiredParam<MaterialPropertyName>("flux", "Neutron Flux");
  return params;
}

ADMixtureOfExpertsElastoViscoplasticStress::ADMixtureOfExpertsElastoViscoplasticStress(
    const InputParameters & parameters)
  : ADComputeStressBase(parameters),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(getADMaterialPropertyByName<RankFourTensor>(_elasticity_tensor_name)),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _mechanical_strain_old(
        getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mu(getADMaterialProperty<Real>("mu")),
    _temp(getADMaterialProperty<Real>("temp")),
    _rhoc(declareADProperty<Real>("rhoc")), 
    _rhow(declareADProperty<Real>("rhow")),
    _rhoc_old(getMaterialPropertyOld<Real>("rhoc")),
    _rhow_old(getMaterialPropertyOld<Real>("rhow")),
    _flux(getADMaterialProperty<Real>("flux")),
    _total_effective_plastic_strain(declareADProperty<Real>("effective_plastic_strain")),
    _total_effective_plastic_strain_old(getMaterialPropertyOld<Real>("effective_plastic_strain")),
    _initial_rhoc(getParam<Real>("initial_rhoc")),
    _initial_rhow(getParam<Real>("initial_rhow")),
    effective_plastic_strain_increment(0.0),
    rhoc_rate(0.0),
    rhow_rate(0.0),
    _incremental_plastic_strain_tensor(declareADProperty<RankTwoTensor>("incremental_plastic_strain_tensor")),
    _incremental_plastic_strain_tensor_old(getMaterialPropertyOld<RankTwoTensor>("incremental_plastic_strain_tensor"))

{
  // Checking if the model file exists and load it during construction
  std::string model_path = "/home/zaheen/PROJECTS/moose/modules/solid_mechanics/src/materials/model_scripted.pt";

  if (!std::filesystem::exists(model_path)) {
      std::cerr << "Model file does not exist at the specified path: " << model_path << std::endl;
      throw std::runtime_error("Model file not found");
  }

  try {
      // Load the TorchScript model
      model = torch::jit::load(model_path);
      std::cout << "Model loaded successfully." << std::endl;
  } catch (const c10::Error &e) {
      std::cerr << "Error loading the model: " << e.what() << std::endl;
      throw;  // Propagate the exception
  }
}

void
ADMixtureOfExpertsElastoViscoplasticStress::initQpStatefulProperties()
{
  _rhoc[_qp] = _initial_rhoc; 
  _rhow[_qp] = _initial_rhow;
}

void
ADMixtureOfExpertsElastoViscoplasticStress::computeQpStress()
{
  // copying dislocations to a temporary file so that any change can be enforced
  auto _rhoc_old_temp = _rhoc_old[_qp];
  auto _rhow_old_temp = _rhow_old[_qp];

  // Checking bounds of the MoE model inputs and fixing them
  if (_rhoc_old_temp <= 0)
  {
    _rhoc_old_temp = _initial_rhoc; 
  }

  if (_rhow_old_temp <= 0)
  {
    _rhow_old_temp = _initial_rhow;
    
  }

   if (_rhoc_old_temp > 8461801410123.313)
  {
    _rhoc_old_temp = 8461801410123.313; 
  }

  if (_rhow_old_temp > 11999567054170.322)
  {
    _rhow_old_temp = 11999567054170.322;
    
  }

  // Getting the strain increment
  ADRankTwoTensor elastic_strain_increment =
      (_mechanical_strain[_qp] - _mechanical_strain_old[_qp]);

  // Computing trial stress asuming the whole strain increment as elastic
  RankTwoTensor trial_stress = MetaPhysicL::raw_value(_stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment);
  ADRankTwoTensor deviatoric_trial_stress = trial_stress.deviatoric();
  auto norm_dev_stress_squared = deviatoric_trial_stress.doubleContraction(deviatoric_trial_stress);
  auto norm_dev_stress = std::sqrt(norm_dev_stress_squared);
  auto effective_trial_stress = std::sqrt(1.5) * norm_dev_stress; 
  const Real tolerance = 1e-14;

  // Checking for non-zero stress

  if (effective_trial_stress > (1.0 + tolerance)) {

    // Creating an input tensor
    at::Tensor inputs = torch::tensor({{MetaPhysicL::raw_value(effective_trial_stress), MetaPhysicL::raw_value(_temp[_qp]), MetaPhysicL::raw_value(_total_effective_plastic_strain_old[_qp]), 
    MetaPhysicL::raw_value(_rhoc_old_temp), MetaPhysicL::raw_value(_rhow_old_temp), MetaPhysicL::raw_value(_flux[_qp])}}, torch::kFloat64);

    // evaluating the model
    auto model_output = this->eval(inputs);
    auto creep_rate = model_output.index({0, 0});

    // initializing useful varialbles
    Real effective_plastic_strain_increment = 0.0;
    Real residual = 0.0;
    unsigned int k = 0;
    Real tol = 1e-8;
    Real dt = _dt;
    if (dt == 0){dt = 1;}
    
    // Starting the Return-Mapping Newton-Raphson loop
    while (std::abs(residual) > tol || k == 0) {
      if (++k > 1000) {
          std::cerr << "Max iterations reached.\n";
          break;
      }
      
      residual = effective_plastic_strain_increment - creep_rate.item().toDouble() * dt;
      auto jacobian = this->eval_jacobian(inputs);
      auto jac_1st = jacobian.index({0, 0});
      Real dphi_ddeqpl = - 3.0 * MetaPhysicL::raw_value(_mu[_qp]) * jac_1st.item().toDouble();
      Real del_deqpl = (creep_rate.item().toDouble() - effective_plastic_strain_increment / dt) / ((1 / dt) - dphi_ddeqpl);
      effective_plastic_strain_increment += del_deqpl;
      at::Tensor inputs_new = torch::tensor({{MetaPhysicL::raw_value(effective_trial_stress - 3.0 * _mu[_qp] * effective_plastic_strain_increment), MetaPhysicL::raw_value(_temp[_qp]), MetaPhysicL::raw_value(_total_effective_plastic_strain_old[_qp]), 
      MetaPhysicL::raw_value(_rhoc_old_temp), MetaPhysicL::raw_value(_rhow_old_temp), MetaPhysicL::raw_value(_flux[_qp])}}, torch::kFloat64);
      inputs = inputs_new;
      model_output = this->eval(inputs);
      auto new_rate = model_output.index({0, 0});
      creep_rate = new_rate;
    }


    // Updating dislocations and total plastic strains for next MoE input
    auto rhoc_rate_temp = model_output.index({0, 1});
    auto rhow_rate_temp = model_output.index({0, 2});

    rhoc_rate = rhoc_rate_temp.item().toDouble();
    rhow_rate = rhow_rate_temp.item().toDouble();

    _total_effective_plastic_strain[_qp] = _total_effective_plastic_strain_old[_qp] + effective_plastic_strain_increment; 
    _rhoc[_qp] = _rhoc_old_temp + rhoc_rate * dt;
    _rhow[_qp] = _rhow_old_temp + rhow_rate * dt;
    
    // Getting the plastic strain tensor
    _incremental_plastic_strain_tensor[_qp] = deviatoric_trial_stress * (1.5 * effective_plastic_strain_increment / effective_trial_stress);
    
    // getting elastic strain 
    elastic_strain_increment = elastic_strain_increment - _incremental_plastic_strain_tensor[_qp];
  }
  // Computing actual stress after returning back
  _stress[_qp] = _stress_old[_qp] + _elasticity_tensor[_qp] * elastic_strain_increment;

}


torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::eval(const torch::Tensor& tensor){
  
  auto varx1 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 0}), -1.0, 1.0, 0.119132960099, 299.985430264);
  auto varx2 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 1}), -1.0, 1.0, 600.0690046472204, 1099.9901302093344);
  auto varx3 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 2}), -1.0, 1.0, -20.0, -1.3010306066481285, 1e-20);
  auto varx4 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 3}), -1.0, 1.0, 3.8772211664273257, 12.927462828664899, 0.0);
  auto varx5 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 4}), -1.0, 1.0, 12.644195880217044, 13.079165576932896, 0.0);
  auto varx6 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 5}), -1.0, 1.0, 1.018454806316354e-09, 9.998075715064236e-07);

 
  at::Tensor transformed = torch::tensor({{varx1.item<double>(), varx2.item<double>(), varx3.item<double>(), varx4.item<double>(), varx5.item<double>(), varx6.item<double>()}}, torch::kFloat64);

  at::Tensor output = model.forward({transformed}).toTensor();


  auto out1 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_inverse_transform(output.index({0, 0}), 0.0, 1.0, -13.436071343451458, 3.905343109204577, 0.0);
  auto out2 = ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_transform(output.index({0, 1}), 1.0, -16.17746470707994, 17.114161436245546, 0.46834836458280193);
  auto out3 = ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_transform(output.index({0, 2}), 1.0, -18.45887052858157, 15.864119682881029, -1.297375422850271);


  at::Tensor evaluated = torch::tensor({{out1.item<double>(), out2.item<double>(), out3.item<double>()}}, torch::kFloat64);

  return evaluated;
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::eval_jacobian(const torch::Tensor& tensor){
  auto jax1 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_forward_derivative(tensor.index({0, 0}), -1.0, 1.0, 0.119132960099, 299.985430264);
  auto jax2 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_forward_derivative(tensor.index({0, 1}), -1.0, 1.0, 600.0690046472204, 1099.9901302093344);
  auto jax3 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_forward_derivative(tensor.index({0, 2}), -1.0, 1.0, -20.0, -1.3010306066481285, 1e-20);
  auto jax4 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_forward_derivative(tensor.index({0, 3}), -1.0, 1.0, 3.8772211664273257, 12.927462828664899, 0.0);
  auto jax5 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_forward_derivative(tensor.index({0, 4}), -1.0, 1.0, 12.644195880217044, 13.079165576932896, 0.0);
  auto jax6 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_forward_derivative(tensor.index({0, 5}), -1.0, 1.0, 1.018454806316354e-09, 9.998075715064236e-07);

  at::Tensor Jin = torch::tensor({{jax1.item<double>(), jax2.item<double>(), jax3.item<double>(), jax4.item<double>(), jax5.item<double>(), jax6.item<double>()}}, torch::kFloat64).view(-1);

  auto varx1 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 0}), -1.0, 1.0, 0.119132960099, 299.985430264);
  auto varx2 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 1}), -1.0, 1.0, 600.0690046472204, 1099.9901302093344);
  auto varx3 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 2}), -1.0, 1.0, -20.0, -1.3010306066481285, 1e-20);
  auto varx4 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 3}), -1.0, 1.0, 3.8772211664273257, 12.927462828664899, 0.0);
  auto varx5 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(tensor.index({0, 4}), -1.0, 1.0, 12.644195880217044, 13.079165576932896, 0.0);
  auto varx6 = ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(tensor.index({0, 5}), -1.0, 1.0, 1.018454806316354e-09, 9.998075715064236e-07);

  at::Tensor p = torch::tensor({{varx1.item<double>(), varx2.item<double>(), varx3.item<double>(), varx4.item<double>(), varx5.item<double>(), varx6.item<double>()}}, torch::kFloat64);

  p = p.repeat({3, 1});
  p.requires_grad_(true);

  torch::Tensor  y = model.forward({p}).toTensor();

  torch::Tensor eye_matrix = torch::eye(3, torch::kDouble);

  y.backward(eye_matrix);
  torch::Tensor Jmod = p.grad();

  auto phi = y[0].unsqueeze(0);

  auto out1 = ADMixtureOfExpertsElastoViscoplasticStress::logtransform_inverse_derivative(phi.index({0, 0}), 0.0, 1.0, -13.436071343451458, 3.905343109204577, 0.0);
  auto out2 = ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_derivative(phi.index({0, 1}), 1.0, -16.17746470707994, 17.114161436245546, 0.46834836458280193);
  auto out3 = ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_derivative(phi.index({0, 2}), 1.0, -18.45887052858157, 15.864119682881029, -1.297375422850271);

  at::Tensor Jout = torch::tensor({{out1.item<double>(), out2.item<double>(), out3.item<double>()}}, torch::kFloat64).view(-1);

  torch::Tensor J = Jout.unsqueeze(1) * (Jmod * Jin);

  return J;

}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::minmax_transform(const torch::Tensor& x, double lb, double ub, double xmin, double xmax){

  auto a = ub - lb;
  auto b = lb;
  return a * (x - xmin) / (xmax - xmin) + b;
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::logtransform_transform(const torch::Tensor& x, double lb, double ub, double logxmin, double logxmax, double eps){
  auto a = ub - lb;
  auto b = lb;
  return a * (torch::log10(x + eps) - logxmin) / (logxmax - logxmin) + b;
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::logtransform_inverse_transform(const torch::Tensor& z, double lb, double ub, double logxmin, double logxmax, double eps){
  auto a = ub - lb;
  auto b = lb;
  return torch::pow(10, (logxmax - logxmin) * (z - b) / a + logxmin) - eps;
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_transform(const torch::Tensor& z, double a, double zmin, double zmax, double zbar){
  auto x = (zmax - zmin) * z / (2 * a) + zbar;
  return torch::sign(x) * (torch::pow(10, torch::abs(x)) - 1.0);
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::minmax_forward_derivative(const torch::Tensor& x, double lb, double ub, double xmin, double xmax){
  auto a = ub - lb;
  auto b = lb;
  return torch::tensor(a / (xmax - xmin));
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::logtransform_forward_derivative(const torch::Tensor& x, double lb, double ub, double logxmin, double logxmax, double eps){
  auto a = ub - lb;
  auto b = lb;
  return a  / (logxmax - logxmin) * (1.0/((x + eps) * torch::log(torch::tensor(10.0))));
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::logtransform_inverse_derivative(const torch::Tensor& z, double lb, double ub, double logxmin, double logxmax, double eps){
  auto a = ub - lb;
  auto b = lb;
  auto u = torch::pow(10, (logxmax - logxmin) * (z - b) / a + logxmin);
  return u * torch::log(torch::tensor(10.0)) * (logxmax - logxmin) / a;
}

torch::Tensor 
ADMixtureOfExpertsElastoViscoplasticStress::symlog_inverse_derivative(const torch::Tensor& z, double a, double zmin, double zmax, double zbar){
  auto x = (zmax - zmin) * z / (2 * a) + zbar;
  auto u = torch::pow(10, torch::abs(x));
  return u * torch::log(torch::tensor(10.0)) * (zmax - zmin) / (2 * a);
}


