[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [cube]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    nx = 1
    ny = 1
    nz = 1
    elem_type = HEX8
  []
  [pin]
    type = ExtraNodesetGenerator
    input = cube
    new_boundary = pin
    coord = '0 0 0'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [elasto_viscoplastic]
        # new_system = true
        add_variables = true
        strain = SMALL
        incremental = true
        generate_output = "stress_xx stress_yy stress_zz
                           stress_xy stress_xz stress_yz
                           mechanical_strain_xx mechanical_strain_yy mechanical_strain_zz
                           mechanical_strain_xy mechanical_strain_xz mechanical_strain_yz"
        # additional_generate_output = 'vonmises_cauchy_stress'
        use_automatic_differentiation = true
      []
    []
  []
[]

[Materials]
  [elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    lambda = 116762.4944714728
    shear_modulus = 60150.37593984962
  []
  [stress]
    type = ADMixtureOfExpertsElastoViscoplasticStress
    # lambda = 116762.4944714728
    mu = 60150.37593984962
    temp = 789
    initial_rhoc = 5.85e12
    initial_rhow = 8.66e12
    flux = 3.82e-8
    output_properties = 'rhoc rhow effective_plastic_strain'
    outputs = exodus
  []
  # [stress_zz]
	# 	type = ADRankTwoCartesianComponent
	# 	rank_two_tensor = stress
	# 	property_name = stress_zz
	# 	index_i = 2
	# 	index_j = 2
	# []
[]

[Postprocessors]
	[avg_stress_zz]
		type = ADElementAverageMaterialProperty
		mat_prop = stress_zz
	[]
  [rhoc]
		type = ADElementAverageMaterialProperty
		mat_prop = rhoc
    # execute_on = LINEAR
	[]
  [rhow]
		type = ADElementAverageMaterialProperty
		mat_prop = rhow
    # execute_on = LINEAR
	[]
  [plastic_strain]
		type = ADElementAverageMaterialProperty
		mat_prop = effective_plastic_strain
    # execute_on = LINEAR
	[]
[]
# [NodalKernels]
#   [force]
#     type = UserForcingFunctionNodalKernel
#     variable = disp_z
#     function = -12.6
#     boundary = 'front'
#   []
# []

[BCs]
  [fix_x]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'pin'
    # preset = true
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'pin'
    # preset = true
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'back'
    # preset = true
  []
  [Pressure]
    [front]
      boundary = 'front'
      # type = Pressure
      # variable = disp_z
      factor =  -128.7 # 130.57899616142 MPa
      # preset = true
    []
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'lu'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       nonzero              '
  # snesmf_reuse_base = false
  automatic_scaling = true
  end_time = 5
  dt = 1.9375
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 7.75
  #   optimal_iterations = 100
  #   cutback_factor = 0.5
  #   growth_factor = 2
  # []
  # nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  nl_forced_its = 4
  nl_max_its = 150
  residual_and_jacobian_together = true
  # line_search = None
[]

[Outputs]
  # print_linear_residuals = false
  file_base = MoE_Tension_SciPy_1287_test
  exodus = true
  csv = true
  perf_graph = true
  # [pgraph]
  #   type = PerfGraphOutput
  #   # execute_on = 'initial final'  # Default is "final"
  #   level = 2                    # Default is 1
  #   heaviest_branch = true        # Default is false
  #   heaviest_sections = 10        # Default is 0
  # []
[]
