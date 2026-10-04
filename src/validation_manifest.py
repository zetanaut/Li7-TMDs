"""Required evidence IDs for the implemented baseline and reserved full profile."""

# This ordered registry freezes the historical check-to-ID correspondence.
LEGACY_SYMBOLIC_LABELS = (
    'three-dimensional quadrupole trace',
    'three-dimensional octupole traces',
    'spin-density inverse map on all 16 Hermitian basis matrices',
    'octupole longitudinal eigenvalues',
    'octupole helicity-one matrix',
    'quadrupole helicity-two matrix',
    'STF momentum tensors through rank five',
    'STF of rank-two product vanishes',
    'dual STF of rank-two product vanishes',
    'dual STF of rank-three double contraction vanishes',
    'quark complex low branch m=1',
    'quark complex high branch m=1',
    'gluon complex low branch m=1',
    'gluon complex high branch m=1',
    'sine invariant m=1',
    'quark complex low branch m=2',
    'quark complex high branch m=2',
    'gluon complex low branch m=2',
    'gluon complex high branch m=2',
    'sine invariant m=2',
    'quark complex low branch m=3',
    'quark complex high branch m=3',
    'gluon complex low branch m=3',
    'gluon complex high branch m=3',
    'sine invariant m=3',
    'quark catalogue has 32 entries',
    'quark 64-by-32 tensor map has rank 32',
    'quark multipole counts 2,6,10,14',
    'gluon catalogue has 32 entries',
    'gluon 64-by-32 tensor map has rank 32',
    'gluon multipole counts 2,6,10,14',
    'quark projection Gram matrix',
    'fragmentation projection Gram matrix',
    'electromagnetic trace channel 0,0',
    'electromagnetic trace channel 0,2',
    'electromagnetic trace channel 0,3',
    'electromagnetic trace channel 1,0',
    'electromagnetic trace channel 1,2',
    'electromagnetic trace channel 1,3',
    'electromagnetic trace channel 2,0',
    'electromagnetic trace channel 2,2',
    'electromagnetic trace channel 2,3',
    'electromagnetic trace channel 3,0',
    'electromagnetic trace channel 3,2',
    'electromagnetic trace channel 3,3',
    'SIDIS scalar-helicity-Collins master contraction',
    'ideal transition vector moment',
    'ideal transition quadrupole moment',
    'ideal transition octupole moment',
 )
LEGACY_IDS = {label: f"legacy.symbolic.{i:03d}" for i, label in enumerate(LEGACY_SYMBOLIC_LABELS, 1)}
BORN_CHECKS = (
    'on_shell_and_conservation', 'photon_Ward', 'gluon_Ward',
    'Hermitian_hard_matrix', 'positive_semidefinite_hard_matrix',
)
SOFTWARE_TESTS = ('test_symbolic_report', 'test_born_report', 'test_grid')
BASELINE_REQUIRED = tuple(LEGACY_IDS.values()) + tuple(f'born.{name}' for name in BORN_CHECKS) + ('grid.complete',) + tuple(f'software.{name}' for name in SOFTWARE_TESTS)

# Suite IDs reserve the scope of the complete manuscript validation. They are
# not synthetic PASS checks. Full remains incomplete until these suites exist.
FULL_PENDING_SUITES = (
    'suite.spin_independent', 'suite.stf_independent', 'suite.gluon_dictionary',
    'suite.correlators_projectors', 'suite.sidis_complete', 'suite.dy_complete',
    'suite.universality', 'suite.gluon_response', 'suite.born_independent_scan',
    'suite.positivity_collinear_fourier', 'suite.negative_controls',
)
PROFILE_REQUIRED = {
    'baseline': BASELINE_REQUIRED,
    'full': BASELINE_REQUIRED + FULL_PENDING_SUITES,
}

# Checkpoint-2 requirements are frozen independently of the functions that
# execute them. A missing implementation therefore remains a missing ID.
SPIN_FOUNDATION_IDS = (
    'spin.ladder_commutator','spin.casimir','spin.derived_subtractions','spin.traces',
    'spin.spherical_orthonormal','spin.spherical_conjugation','spin.anchor_conversion',
    'spin.helicity_selection','spin.spherical_generators','spin.spherical_ladder','spin.cartesian_density_inverse',
    'spin.spherical_density_inverse','spin.positive_complex_states','spin.trace_metrics',
    'spin.full_octupole_metric','spin.cross_rank_orthogonality','spin.cartesian_permutation',
    'spin.metric_inverse_weights','spin.octupole_spectrum','spin.preparation_positivity',
    'spin.pure_octupole','spin.rotation_calibration','spin.rotation_coherence',
    'spin.rotation_octupole',
    *(f'spin.direction_unit_{i}' for i in range(7)),
    'spin.rotated_cartesian_tensor','spin.rotated_components','spin.transverse_llt_ttt','spin.real_octupole_basis',
    'spin.seven_direction_tomography','spin.response_recovery',
    'spin.lower_rank_cancellation','spin.population_inverse','spin.nmr_strengths',
    'spin.nmr_inverse','spin.ls_isometry','spin.ls_projection',
)
STF_FOUNDATION_IDS = (
    'stf.rank_0.independent',
    *(f'stf.rank_{n}.{kind}' for n in range(1,6)
      for kind in ('independent', 'trace', 'dual', 'legacy_mass', 'mass_conversion',
                   'rotation', 'homogeneity', 'norm', 'zero_momentum') if not (n==1 and kind=='trace')),
    *(f'stf.rank_{n}.scalar_sine' for n in range(1,4)),
    'stf.rank_two_products','stf.rank_three_product_dual',
    'stf.dyadic_null','stf.dyadic_residual','stf.dyadic_pure_trace',
)
COVARIANT_FOUNDATION_IDS = ('target.moment_coordinates',) + tuple(
    f'{species}.{kind}' for species in ('quark','gluon') for kind in (
        'semantic_catalogue','helicity_comparison','channel_counts','symbolic_covariants',
        'zero_momentum','second_exact_point','parity_bound','rotation_covariance',
        'joint_expectation','rank_certificate','analytic_projectors',
        'independent_linear_recovery'))
DICTIONARY_FOUNDATION_IDS = (
    'gluon.dictionary_complete','gluon.dictionary_mass_polarization',
    'gluon.dictionary_derived_conversion','gluon.dictionary_rank_kernel',
)
FOUNDATION_NEGATIVE_TESTS = (
    'test_certificate_entry_label_digest_and_sign',
    'test_duplicate_omit_and_rank_preserving_sign',
    'test_imaginary_density_and_multiplicity',
    'test_lt_shift_and_polarization_mapping',
    'test_parity_dual_and_epsilon',
    'test_projector_singular_domain',
    'test_reversed_rotation_phase',
    'test_stf_dyadic_and_mass',
    'test_subtraction_coefficients',
    'test_transverse_octupole_has_llt',
)
FOUNDATION_SOFTWARE_IDS = tuple('software.foundation_negative.'+name for name in FOUNDATION_NEGATIVE_TESTS) + (
    'software.foundation_evidence_integrity',
)
FOUNDATION_REQUIRED = (BASELINE_REQUIRED + SPIN_FOUNDATION_IDS + STF_FOUNDATION_IDS +
                       COVARIANT_FOUNDATION_IDS + DICTIONARY_FOUNDATION_IDS +
                       FOUNDATION_SOFTWARE_IDS)
# Later scientific suites stay required for the full profile.
LATER_PENDING_SUITES = (
    'suite.sidis_complete','suite.dy_complete','suite.universality',
    'suite.gluon_response','suite.born_independent_scan',
    'suite.positivity_collinear_fourier',
)
PROFILE_REQUIRED['foundations'] = FOUNDATION_REQUIRED
PROFILE_REQUIRED['full'] = FOUNDATION_REQUIRED + LATER_PENDING_SUITES
FULL_PENDING_SUITES = LATER_PENDING_SUITES

# Checkpoint-3 rows are independently curated data, not emitted by the
# Cartesian production generator. Every row also requires an integral result.
from response_fixtures import SIDIS_ROWS,DY_ROWS,row_id
PROCESS_FIXED_IDS=(
    'process.dirac.algebra','process.dirac.sidis_trace','process.dirac.dy_trace',
    'process.sidis.normalization','process.sidis.flavors','process.dy.current','process.dy.normalization',
    'process.dy.flavors','process.target.preparations','process.spin_differences',
    'process.convolution.momentum_sign',
    'process.reversal.density_links','process.reversal.conditional_evolution',
)
PROCESS_ROW_IDS=tuple(row_id(p,row) for p,rows in (('SIDIS',SIDIS_ROWS),('DY',DY_ROWS)) for row in rows)
PROCESS_INTEGRAL_IDS=tuple('integral.'+check_id for check_id in PROCESS_ROW_IDS)
REVERSAL_IDS=tuple(f'process.reversal.{species}.{channel}.{K}{m}.{n}'
    for species in ('quark','gluon') for K in range(4) for m in range(K+1)
    for channel,n in (
      *((('f',m),) if m or K%2==0 else ()),
      *((('g',m),) if m or K%2==1 else ()),
      *((('h',n) for n in ((1,) if species=='quark' and m==0 else
                             (2,) if species=='gluon' and m==0 else
                             (abs(m-(1 if species=='quark' else 2)),m+(1 if species=='quark' else 2))))),
    ))
PROCESS_NEGATIVE_TESTS=(
    'test_collins_sign_and_high_branch',
    'test_target_factor_and_azimuth_conversion',
    'test_momentum_sign_and_jacobian',
    'test_antiquark_beam_and_exchange',
    'test_antiunitary_and_links',
    'test_transverse_prep_and_denominator',
    'test_whole_asymmetry_sign',
)
PROCESS_SOFTWARE_IDS=tuple('software.process_negative.'+name for name in PROCESS_NEGATIVE_TESTS)+(
    'software.process_evidence_integrity',)
PROCESS_REQUIRED=(FOUNDATION_REQUIRED+PROCESS_FIXED_IDS+PROCESS_ROW_IDS+
                  PROCESS_INTEGRAL_IDS+REVERSAL_IDS+PROCESS_SOFTWARE_IDS)
PROFILE_REQUIRED['quark-processes']=PROCESS_REQUIRED
LATER_PENDING_SUITES=('suite.gluon_response','suite.born_independent_scan',
                      'suite.positivity_collinear_fourier')
PROFILE_REQUIRED['full']=PROCESS_REQUIRED+LATER_PENDING_SUITES
FULL_PENDING_SUITES=LATER_PENDING_SUITES
